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

beta_dm_regression <- function(exper, i, rolling=FALSE, solve=TRUE) {
  sf <- colSums(assay(exper, "counts"))[order(exper$name)]
  sf <- matrix(
    sf,
    ncol = 4
  )
  if (rolling && i %in% c("2L", "2R", "3L", "3R", "4", "X", "Y")) {
    final_chr_y <- -1 + findInterval(8, as.numeric(rowData(exper)$seqnames))
    blk <- assay(exper, "counts")[
      seq(max(1, i - 3), min(final_chr_y, i + 3)),
      order(exper$name)
    ]
    obs <- matrix(
      blk,
      ncol = 4
    )
    sfExtend <- sparseMatrix(
      i = seq(nrow(obs)),
      j = rep(seq(nrow(sf)), each = nrow(obs) / nrow(sf))
    )
    sf <- as.matrix(sfExtend %*% sf)
  } else {
    matrix_row <- assay(exper, "counts")[i, ]
    obs <- matrix(
      matrix_row[order(exper$name)],
      ncol = 4
    )
  }
  n <- rowSums(obs)
  obs <- obs[n > 0, ]
  n <- n[n > 0]
  ll_post <- function(v, theta) {
    pb <- diff(pbeta(c(0, 0.25, 0.5, 0.75, 1), v[1], v[2]))
    log_prior <- sum(dgamma(v, 4, scale=0.5, log = TRUE))
    expected_value_shape <- sf %*% diag(pb)
    expected_value_shape <- expected_value_shape / rowSums(expected_value_shape)
    dir_shape <- expected_value_shape / theta
    log_lik <- sum(ddirmnom(obs, n, dir_shape, log = TRUE))
    -log_prior - log_lik
  }
  # Start with a typical value of theta (err on the low side).
  theta <- 0.05
  # Initial parameter value is the expected value of alpha and beta in the
  # prior.
  par <- optim(c(2, 2), ll_post, theta = theta, method = "L-BFGS-B", control = list(maxit = 1000), lower = c(0.01, 0.01))$par
  theta <- optim(theta, ll_post, v = par, method = "L-BFGS-B", control = list(maxit = 1000), lower = 1e-4)$par
  par <- optim(par, ll_post, theta = theta, method = "L-BFGS-B", control = list(maxit = 1000), lower = c(0.01, 0.01))$par
  precision <- optimHess(par, ll_post, theta = theta)
  if (solve) {
    sigma <- solve(precision)
    moment <- \(x, y) matrix(as.numeric((x - 1) / (x + y - 2)) * dmvnorm(cbind(as.numeric(x), as.numeric(y)), par, sigma), nrow = nrow(x))
    mydens <- \(x, y) matrix(dmvnorm(cbind(as.numeric(x), as.numeric(y)), par, sigma), nrow = nrow(x))
    est <- integral2(
      moment, 1, 20, 1, 20,
      abstol = 1e-3
    )$Q / integral2(
      mydens, 1, 20, 1, 20,
      abstol = 1e-3
    )$Q
    list(est = est, par = par, theta = theta, fisher_information = precision)
  } else {
    list(par = par, theta = theta, fisher_information = precision)
  }
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

laplace_approx_expectation_ratio_old_integral2 <- function(elem) {
  par <- elem$par
  sigma <- solve(elem$fisher_information)
  moment <- \(x, y) matrix(as.numeric((x - 1) / (x + y - 2)) * dmvnorm(cbind(as.numeric(x), as.numeric(y)), par, sigma), nrow = nrow(x))
  mydens <- \(x, y) matrix(dmvnorm(cbind(as.numeric(x), as.numeric(y)), par, sigma), nrow = nrow(x))
  est <- integral2(
    moment, 1, 100, 1, 100
  )$Q / integral2(
    mydens, 1, 100, 1, 100
  )$Q
}

laplace_approx_expectation_ratio <- function(elem) {
  par <- elem$par
  sigma <- solve(elem$fisher_information)
  moment <- \(x, y) matrix(as.numeric((x - 1) / (x + y - 2)) * dmvnorm(cbind(as.numeric(x), as.numeric(y)), par, sigma), nrow = nrow(x))
  mydens <- \(x, y) matrix(dmvnorm(cbind(as.numeric(x), as.numeric(y)), par, sigma), nrow = nrow(x))
  numer <- integral(
    # Vectorize(purrr::partial(integrate_gauss2d_arc, origin=c(1,1), mean=par, sigma=sigma)),
    Vectorize(\(angle) sin(angle)/(sin(angle)+cos(angle)) * integrate_gauss2d_arc(c(1,1), angle, par, sigma)),
    0, pi/2
  )
  denom <- pmvnorm(
    lower=c(1, 1),
    mean = elem$par,
    sigma = sigma
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

analyze_beta_repli.experiment_Germline_chr <- function(
    repli.experiment_Germline_chr) {
  library(future)
  library(future.apply)
  library(extraDistr)
  library(pracma)
  plan(multicore, workers=12)

  fits <- future_sapply(
    seq(nrow(repli.experiment_Germline_chr)),
    \(i) tryCatch(beta_dm_regression(repli.experiment_Germline_chr, i, solve=F), error=\(c) NULL),
    simplify=FALSE
  )
  gr <- tar_read(repli.granges_nos_E2_chr) %>% unlist
  gr$score <- as.numeric(sapply(fits, \(l) if (is.null(l$est)) NA else (1-2*l$est))) %>%
    replace(is.na(.), 0)
  gr
}
