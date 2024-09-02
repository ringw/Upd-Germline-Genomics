dm_regression_gaussian_plate <- function(
  exper, center_i, wts
) {
  beta_regression_size_factors <- exper@metadata[["beta_regression_size_factors"]]
  beta_regression_size_factors <- beta_regression_size_factors /
    rowMeans(beta_regression_size_factors)
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
  ll_post <- function(pb) {
    theta <- 1/sum(pb)
    log_prior <- dgamma(theta, 1, scale=0.5, log = TRUE)
    dirparam <- beta_regression_size_factors %*% diag(x = pb)
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
  optim(
    rep(1, exper@metadata$beta_regression_num_fractions),
    ll_post,
    method = "L-BFGS-B",
    control = list(maxit = 1000),
    lower = rep(0.01, exper@metadata$beta_regression_num_fractions)
  )
}