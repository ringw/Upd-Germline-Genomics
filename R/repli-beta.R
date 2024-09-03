beta_prior_draws <- function(n = 200, shape = 4, rate = 0.5) {
  alpha <- rgamma(n * 2, shape, rate) %>%
    subset(. > 1) %>%
    head(n)
  beta <- rgamma(n * 2, shape, rate) %>%
    subset(. > 1) %>%
    head(n)
  x <- seq(0, 1, by = 0.02)
  y <- tibble(alpha, beta) %>%
    rowwise() %>%
    reframe(
      name = str_glue("{alpha}@{beta}"),
      x,
      y = dbeta(x, alpha, beta)
    )
  ggplot(y, aes(x, y, group = name, color = name)) +
    geom_line() +
    scale_color_hue(guide = NULL) +
    theme_cowplot() +
    coord_cartesian(c(0, 1), c(-0.01, 7.5), expand = F)
}

init_beta_dm_experiment <- function(sce) {
  # We will pretend that the SummarizedExperiment is a SingleCellExperiment
  # - having an optional colData entry named sizeFactor.
  sf <- sce$sizeFactor
  if (is.null(sf))
    sf <- colSums(assay(sce, "counts")) / 1000 / 1000
  num_fractions <- length(levels(sce$replication_value))
  sce@metadata$beta_regression_cuts <- seq(
    0, 1, length.out=1+num_fractions
  )
  # Zero is not allowed in extraDistr Dirichlet! (For an entry which is observed
  # to always be zero). Infinitesimal value works well, and is so small that it
  # does not influence the optimization process.
  beta_regression_size_factors <- matrix(
    1e-40,
    nrow = length(levels(sce$rep)),
    ncol = num_fractions,
    dimnames = list(levels(sce$rep), levels(sce$replication_value))
  )
  beta_regression_feed_y <- matrix(
    0,
    nrow = ncol(sce),
    ncol = length(beta_regression_size_factors),
    dimnames = list(colnames(sce), NULL)
  )
  for (i in seq(ncol(sce))) {
    rep <- sce$rep[i]
    replication_value <- sce$replication_value[i]
    beta_regression_size_factors[
      as.numeric(rep),
      as.numeric(replication_value)
    ] <- sf[i]
    beta_regression_feed_y[
      i,
      as.numeric(rep)
      + length(levels(sce$rep)) * (as.numeric(replication_value) - 1)
    ] <- 1
  }
  sce@metadata$beta_regression_num_fractions <- num_fractions
  sce@metadata$beta_regression_size_factors <- beta_regression_size_factors
  sce@metadata$beta_regression_feed_y <- beta_regression_feed_y
  sce
}

beta_dm_regression <- function(exper, i) {
  beta_regression_cuts <- exper@metadata[["beta_regression_cuts"]]
  beta_regression_size_factors <- exper@metadata[["beta_regression_size_factors"]]
  Y <- dot(
    exper@metadata[["beta_regression_feed_y"]],
    assay(exper, "counts")[i, ]
  ) %>%
    matrix(
      nrow = length(levels(exper$rep)),
      ncol = length(levels(exper$replication_value)),
      dimnames = list(
        levels(exper$rep),
        levels(exper$replication_value)
      )
    )
  n <- rowSums(Y)
  ll_post <- function(v, theta) {
    log_prior <- (
      sum(dgamma(v, 9*4, scale=1/12, log = TRUE))
      + dgamma(theta, 1, scale=0.5, log = TRUE)
    )
    pb <- diff(pbeta(beta_regression_cuts, v[1], v[2]))
    pb <- pmax(pb, 1e-40)
    mu_unnorm <- beta_regression_size_factors %*% diag(x = pb)
    dirparam <- (
      mu_unnorm
      /
      rowSums(mu_unnorm)
      /
      theta
    )
    log_lik <- sum(ddirmnom(Y, n, dirparam, log = TRUE))
    -log_prior - log_lik
  }
  # Start with a typical value of theta (err on the low side).
  theta <- 0.01
  # Initial parameter value is the expected value of alpha and beta in the
  # prior.
  par <- optim(c(2, 2), ll_post, theta = theta, method = "L-BFGS-B", control = list(maxit = 1000), lower = c(0.01, 0.01))$par
  theta <- optim(theta, ll_post, v = par, method = "L-BFGS-B", control = list(maxit = 1000), lower = 1e-4)$par
  par <- optim(par, ll_post, theta = theta, method = "L-BFGS-B", control = list(maxit = 1000), lower = c(0.01, 0.01))$par
  precision <- optimHess(par, ll_post, theta = theta)
  list(par = par, theta = theta, fisher_information = precision)
}

# Pulls the "Y" matrix of DM-distributed data (tuples of sample counts). As we
# are interested in a rolling view of "Y", then "i" is a numeric vector of the
# features in the SummarizedExperiment to extract.
beta_dm_pull_expr_plate <- function(exper, i)
  crossprod(
    exper@metadata[["beta_regression_feed_y"]],
    t(assay(exper, "counts")[i,, drop=F])
  ) %>%
    array(
      dim = c(
        length(levels(exper$rep)),
        length(levels(exper$replication_value)),
        length(i)
      ),
      dimnames = list(
        levels(exper$rep),
        levels(exper$replication_value),
        i
      )
    ) %>%
    aperm(c(1, 3, 2)) %>%
    matrix(
      nrow = length(levels(exper$rep)) * length(i),
      ncol = length(levels(exper$replication_value)),
      dimnames = list(
        with(
          cross_join(tibble(index=i), tibble(rep=levels(exper$rep))),
          str_glue("{rep}@{index}")
        )
      )
    )

# Regression: In SummarizedExperiment, extract a view of the sequential
# features, centered on the numeric index "center_i", according to "wts".
# A logistic link function is applied. Then the predicted log-odds value can be
# transformed, after producing a pseudotime value `(alpha-1)/(alpha+beta-2)`
# (Beta distribution). The logit transform is applied, and the pseudotime value
# is multiplied by scale and translated by center. Inverting the transformation,
# we get an updated alpha and beta. The other degree of freedom of the Beta dist
# is handled by keeping `alpha+beta` (concentration of the Beta) invariant. If
# alpha <= 1 or beta <= 1, then the pseudotime value will be extreme and we
# won't apply any translation to this predicted value.
beta_dm_regression_gaussian_plate <- function(
  exper, center_i, wts,
  xform_scale = 1, xform_center = 0
) {
  beta_regression_cuts <- exper@metadata[["beta_regression_cuts"]]
  beta_regression_size_factors <- exper@metadata[["beta_regression_size_factors"]]
  i <- seq(
    pmax(1, center_i - floor(length(wts) / 2)),
    pmin(nrow(exper), center_i + floor(length(wts) / 2))
  )
  wts <- wts[
    seq(
      1 + floor(length(wts) / 2) + min(i) - center_i,
      1 + floor(length(wts) / 2) + max(i) - center_i
    )
  ]
  Y <- beta_dm_pull_expr_plate(exper, i)
  wts <- rep(wts, each=length(levels(exper$rep)))
  n <- rowSums(Y)
  ll_post <- function(v, theta) {
    log_prior <- (
      sum(dgamma(v, 9*2, scale=1/6, log = TRUE))
      + dgamma(theta, 1, scale=0.5, log = TRUE)
    )
    if (v[1] > 1 & v[2] > 1) {
      pseudotime_value <- qlogis((v[1] - 1) / (v[1] + v[2] - 2))
      pseudotime_value <- pseudotime_value * xform_scale + xform_center
      updated_param_value <- plogis(pseudotime_value)
      v <- c(1, 1) + (v[1] + v[2] - 2) * c(updated_param_value, 1 - updated_param_value)
    }
    pb <- diff(pbeta(beta_regression_cuts, v[1], v[2]))
    pb <- pmax(pb, 1e-40)
    mu_unnorm <- beta_regression_size_factors %*% diag(x = pb)
    dirparam <- (
      mu_unnorm
      /
      rowSums(mu_unnorm)
      /
      theta
    )
    log_lik <- dot(
      ddirmnom(
        Y,
        n,
        dirparam[rep(seq(nrow(dirparam)), length(i)), ],
        log = TRUE
      ),
      wts
    )
    -log_prior - log_lik
  }
  # Start with a typical value of theta (err on the low side).
  theta <- 0.01
  # Initial parameter value is the expected value of alpha and beta in the
  # prior.
  par <- optim(c(2, 2), ll_post, theta = theta, method = "L-BFGS-B", control = list(maxit = 1000), lower = c(0.01, 0.01))$par
  theta <- optim(theta, ll_post, v = par, method = "L-BFGS-B", control = list(maxit = 1000), lower = 1e-4)$par
  optim_result <- optim(par, ll_post, theta = theta, method = "L-BFGS-B", control = list(maxit = 1000), lower = c(0.01, 0.01))
  par <- optim_result$par
  precision <- optimHess(par, ll_post, theta = theta)
  list(par = par, theta = theta, objective_value = -optim_result$value, fisher_information = precision)
}

full_quantification_gaussian_plate <- function(rpkm, wts, num_fractions=4) {
  data <- exp(rpkm$Beta[, 1:num_fractions, drop=F])
  data <- data / rowSums(data)
  if (nrow(rpkm$Beta) > 1) {
    fft_length <- 2^(1+ceil_discrete_log(nrow(rpkm$Beta)))
    data <- data %>% replace(!is.finite(.), 0)
    data <- rbind(data, matrix(0, nrow=fft_length-nrow(data), ncol=ncol(data)))
    # As the increased value of the wts does not actually increase the effect that
    # the log-likelihood term has, wts will just sum to 1 for applying filter.
    wts <- wts / sum(wts)
    filt <- c(
      wts[seq(ceil(length(wts)/2), length(wts))],
      rep(0, fft_length - length(wts)),
      wts[seq(1, floor(length(wts)/2))]
    )
    data <- data %>%
      mvfft %>%
      `*`(fft(filt)) %>%
      mvfft(inv=T) %>%
      `/`(fft_length) %>%
      abs %>%
      head(nrow(rpkm$Beta))
  }
  timing_values <- seq(-1, 1, length.out = 1 + 2 * num_fractions)
  timing_values <- timing_values[seq(length(timing_values)-1, 1, by=-2)]
  dot(t(data), timing_values) / dot(t(data), rep(1, num_fractions))
}

beta_dm_regression_calculate_prior <- function() {
  nl_prior <- \(v) -sum(dgamma(v, 4, scale=0.5, log = TRUE))
  par <- optim(c(2, 2), nl_prior)$par
  precision <- optimHess(par, nl_prior)
  list(par = par, theta = 0, fisher_information = precision)
}

laplace_approx_sliding_transform <- function(lst, bw=2) {
  sigmaMu <- sapply(lst, \(elem) solve(elem$fisher_information, elem$par))
  # Laplace posterior s.d. sliding estimate. An element-wise sliding filter on
  # the posterior covariance mats.
  covPoint <- array(
    sapply(lst, \(elem) solve(elem$fisher_information)),
    c(2, 2, length(lst))
  )
  covHat <- array(dim=c(2, 2, length(lst)))
  for (i in 1:2)
    for (j in 1:2)
      covHat[i, j, ] <- density(
        seq_along(lst),
        bw = bw,
        weights = covPoint[i, j, ],
        n = length(lst),
        from = 1,
        to = length(lst)
      )$y
  
  sigmaMuEst <- apply(
    sigmaMu,
    1,
    \(v) density(
      seq_along(lst),
      bw = bw,
      weights = v,
      n = length(lst),
      from = 1
      ,
      to = length(lst)
    )$y
  ) %>%
    t
  muHat <- mapply(
    solve,
    apply(covHat, 3, identity, simplify=FALSE),
    apply(sigmaMuEst, 2, identity, simplify=FALSE)
  )
  mapply(
    \(elem, cov, mu) list(
      par = mu,
      theta = elem$theta,
      fisher_information = solve(cov)
    ),
    lst,
    apply(covHat, 3, identity, simplify=FALSE),
    apply(muHat, 2, identity, simplify=FALSE),
    SIMPLIFY=FALSE
  )
}

laplace_approx_expectation_ratio <- function(elem) {
  par <- elem$par
  sigma <- solve(elem$fisher_information)
  denom <- pmvnorm(
    lower=c(1, 1),
    mean = elem$par,
    sigma = sigma
  )
  if (denom < 1e-10) {
    # There is not much probability mass inside the quadrant. Then because it is
    # a convex distribution, almost all of the probability mass inside the
    # quadrant (closed lower bound indices) is on the x-intercept or the
    # y-intercept.
    if (par[1] > 1) return(1) # Our Late draw from the Beta distribution.
    else if (par[2] > 1) return(0) # Our Early draw from the Beta distribution.

    mean <- par - 1
    v <- eigen(sigma, sym=T)$vectors[, 1]
    if (all(v < 0)) v <- -v
    par.update <- mean + v / max(-v / mean)
    angle.update <- atan2(par.update[1], par.update[2])
    if (angle.update > pi/2)
      # We hit the alpha = 1-intercept, while the beta parameter is < 1. So
      # there is a great amount of probability mass with the mode at alpha,
      # conditioned on alpha >= 1, beta >= 1.
      return(1)
    else if (angle.update < 0)
      return(0)
  }
  # Now we are interested in the moment of sin(a)/(sin(a)+cos(a)) in our
  # quadrant. We didn't actually check that the MAP (par) is in the >1-quadrant,
  # but we are able to compute probabilities in this quadrant regardless,
  # because the values are not infinitesimally small.
  numer <- integral(
    # Vectorize(purrr::partial(integrate_gauss2d_arc, origin=c(1,1), mean=par, sigma=sigma)),
    Vectorize(\(angle) sin(angle)/(sin(angle)+cos(angle)) * integrate_gauss2d_arc(c(1,1), angle, par, sigma)),
    0, pi/2,
    reltol = 1e-3
  )
  res <- numer / as.numeric(denom)
}

integrate_gauss2d_arc_naive <- function(origin, angle, mean, sigma) {
  integral(
    Vectorize(
      \(r) r * dmvnorm(
        origin + r * c(sin(angle), cos(angle)),
        mean,
        sigma
      )
    ),
    0,
    100
  )
}

integrate_gauss2d_arc <- function(origin, angle, mean, sigma) {
  v <- c(sin(angle), cos(angle))
  A <- dot(v, solve(sigma, v))
  B <- dot(v, solve(sigma, mean - origin))
  C <- dot(mean - origin, solve(sigma, mean - origin))
  (
    1/(2*pi*sqrt(det(sigma)))
    *
    exp(B^2/A/2 - C/2)
    *
    (
      1/A * exp(-B^2/2/A)
      +
      B * sqrt(2 * pi) / A^1.5 * pnorm(B / sqrt(A))
    )
  )
}

# Fix Bivariate Normal precision matrix using pseudo-count identity matrix. This
# fixes cases where the Fisher information is not positive definite.
beta_correct_fit_hessian_precision_matrix <- function(fit) {
  eigs <- eigen(fit$fisher_information, sym=T)
  if (eigs$values[2] < 0)
    fit$fisher_information <- fit$fisher_information + c(
      -eigs$values[2] + 0.0001,
      0,
      0,
      -eigs$values[2] + 0.0001
    )
  fit <- fit
}

# Calculate Bayes Factor of our Laplace approximation. A reciprocal of an "inner
# product", of sorts, of our first integrated rays of probability density.
# Analogous to a likelihood ratio test with a shared or independent Beta dist
# location parameter, it scores two model fits for whether the Beta dist
# location is likely to be significantly far apart.
laplace_bayes_factor_polar_integral <- function(posterior_a, posterior_b) {
  posterior_a <- posterior_a %>% beta_correct_fit_hessian_precision_matrix()
  posterior_b <- posterior_b %>% beta_correct_fit_hessian_precision_matrix()
  integral_a <- sapply(
    seq(0, pi/2, length.out=100),
    integrate_gauss2d_arc,
    origin=c(1, 1),
    mean=posterior_a$par,
    sigma=solve(posterior_a$fisher_information)
  )
  integral_b <- sapply(
    seq(0, pi/2, length.out=100),
    integrate_gauss2d_arc,
    origin=c(1, 1),
    mean=posterior_b$par,
    sigma=solve(posterior_b$fisher_information)
  )
  sum(integral_a) * sum(integral_b) / sum(integral_a * integral_b)
}
