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

beta_dm_plain_likelihood_fn <- function(
  Y, beta_regression_cuts, beta_regression_size_factors,
  wts, heads, tails, theta
) {
  n <- rowSums(Y)
  beta_values <- cross_join(
    tibble(cut = beta_regression_cuts),
    tibble(heads = as.numeric(heads), tails = as.numeric(tails))
  )
  pb <- matrix(
    pbeta(beta_values$cut, beta_values$heads, beta_values$tails),
    nrow = length(heads)
  ) %>%
    rowDiffs
  pb <- pmax(pb, 1e-40)
  probs <- purrr::reduce(
    sapply(
      seq(nrow(Y)),
      \(j) {
        y <- Y[j, ]
        sf <- beta_regression_size_factors[j, ]
        mu_unnorm <- t(t(pb) * sf)
        dirparam <- (
          mu_unnorm
          /
          rowSums(mu_unnorm)
          /
          theta
        )
        ddirmnom(
          y,
          n[j],
          dirparam
        )^wts[j]
      },
      simplify=FALSE
    ),
    `*`
  )
}

beta_dm_prior_heads_tails_fn <- function(heads, tails) {
  prior <- dnorm(
    sqrt((heads - 1)^2 + (tails - 1)^2),
    mean = 2.5,
    sd = 1
  )
  prior[!(heads >= 1 & tails >= 1)] <- 0
  prior <- prior
}

# The prior on "heads" and "tails" as a separable filter of "X" (angle) and "Y"
# (heads + tails).
beta_dm_prior_filter <- function(repli.polar.coordinates) {
  angle <- repli.polar.coordinates$X
  rho <- repli.polar.coordinates$Y
  heads <- 1 + rho * sin(angle)
  tails <- 1 + rho * cos(angle)
  M <- beta_dm_prior_heads_tails_fn(heads, tails)
  s <- svd(M)
  Xsum <- ncol(M) / (pi/2)
  Xfactor <- Xsum / sum(abs(s$v[,1]))
  list(
    X = abs(s$v[,1]) * Xfactor,
    Y = s$d[1] / Xfactor * abs(s$u[,1])
  )
}

beta_dm_regression_likelihood <- function(
  exper, center_i, wts,
  repli.polar.coordinates,
  theta = NULL,
  xform_scale = 1, xform_center = 0
) {
  angle <- repli.polar.coordinates$X
  rho <- repli.polar.coordinates$Y
  heads <- 1 + rho * sin(angle)
  tails <- 1 + rho * cos(angle)
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
  beta_regression_size_factors <- beta_regression_size_factors[
    rep(seq(nrow(beta_regression_size_factors)), length(i)),
  ]
  wts <- rep(wts, each=length(levels(exper$rep)))

  # Mutate the amount placed on heads or tails according to the xform.
  update_param <- heads > 1.001 & tails > 1.001
  pseudotime_value <- qlogis(
    ((heads - 1) / (heads + tails - 2))[update_param]
  )
  pseudotime_value <- pseudotime_value * xform_scale + xform_center
  updated_param_value <- plogis(pseudotime_value)
  heads[update_param] <- (
    1 + (heads + tails - 2)[update_param] * updated_param_value
  )
  tails[update_param] <- (
    1 + (heads + tails - 2)[update_param] * (1 - updated_param_value)
  )

  if (is.null(theta)) {
    probs <- beta_dm_plain_likelihood_fn(
      Y, beta_regression_cuts, beta_regression_size_factors,
      wts, heads, tails, theta = 0.02
    ) *
      beta_dm_prior_heads_tails_fn(heads, tails)
    param0 <- which.max(probs)
    if (is.matrix(heads)) {
      i <- param0 %% nrow(heads)
      j <- 1 + floor((param0 - 1) / nrow(heads))
      heads0 <- heads[i, j]
      tails0 <- tails[i, j]
    } else {
      heads0 <- heads[param0]
      tails0 <- tails[param0]
    }
    theta <- optimize(
      \(theta) -log(
        beta_dm_plain_likelihood_fn(
          Y, beta_regression_cuts, beta_regression_size_factors,
          wts, heads0, tails0, theta
        )
      ),
      interval = c(0.001, 20)
    )$minimum %>%
      tryCatch(error = \(e) 0.02)
  }

  probs <- beta_dm_plain_likelihood_fn(
    Y, beta_regression_cuts, beta_regression_size_factors,
    wts, heads, tails, theta
  )
  probs <- if (is.matrix(heads)) {
    probs %>%
      matrix(
        nrow = nrow(heads), ncol = ncol(heads), dimnames = dimnames(heads)
      )
  } else {
    probs
  }
  probs <- probs %>% `attr<-`("theta", value = theta)
}

beta_dm_regression_post_probability <- function(
  exper, center_i, wts,
  repli.polar.coordinates,
  prior_filter,
  theta = NULL,
  xform_scale = 1, xform_center = 0
) {
  grid <- beta_dm_regression_likelihood(
    exper, center_i, wts,
    repli.polar.coordinates,
    theta,
    xform_scale, xform_center
  )
  angle <- repli.polar.coordinates$X[1, ]
  rho <- repli.polar.coordinates$Y[, 1]
  tibble(
    value = sin(angle) / (sin(angle) + cos(angle)),
    prob = colSums(rho * prior_filter$Y * grid),
    theta = attr(grid, "theta")
  )
}

analyze_repli_experiment <- function(
  repli.experiment,
  repli.sliding.weights,
  repli.polar.coordinates,
  repli.prior.distribution,
  xform_scale,
  xform_center,
  theta = NULL
) {
  if (is.null(theta))
    theta <- numeric(0)
  # For each set of indices into the experiment (each chromosome),
  # parallel-apply the beta fitting function with smooth-weight plate.
  mapply(
    \(inds, thetas) {
      exper <- repli.experiment[inds, ]
      future_lapply(
        seq_along(inds),
        \(lookup) with(
          beta_dm_regression_post_probability(
            exper,
            lookup,
            repli.sliding.weights,
            repli.polar.coordinates,
            repli.prior.distribution,
            theta = if (length(thetas)) thetas[lookup] else NULL,
            xform_scale = xform_scale,
            xform_center = xform_center
          ),
          tibble(
            rowname = inds[lookup],
            value,
            prob,
            theta
          )
        )
      ) %>%
        do.call(rbind, .)
    },
    split(
      seq(nrow(repli.experiment)),
      rowData(repli.experiment)$seqnames
    ),
    split(
      theta,
      rowData(repli.experiment)$seqnames
    ),
    SIMPLIFY = FALSE
  ) %>%
    do.call(rbind, .)
}

quantify_repli_experiment <- function(prob, prior) {
  angle <- seq(0, pi/2, length.out=length(prob))
  timing <- sin(angle)/(sin(angle)+cos(angle))
  d_timing_d_angle <- 1/(1 + sin(2*angle))
  sum(timing * prob * prior * d_timing_d_angle) / sum(prob * prior * d_timing_d_angle)
}

bayes_factor_repli_experiment <- function(prob1, prob2, prior) {
  angle <- seq(0, pi/2, length.out=length(prob1))
  d_timing_d_angle <- 1/(1 + sin(2*angle))
  full_model <- sum(
    prob1 * d_timing_d_angle * prior
  ) * sum(
    prob2 * d_timing_d_angle * prior
  )
  reduced_model <- sum(prob1 * prob2 * d_timing_d_angle * prior)
  full_model / reduced_model
}
