devianceResiduals <- function(Y) {
  n <- colSums(Y)
  N <- matrix(n, nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE)
  fit <- glm_gp(
    Y[, n != 0],
    size_factors = n[n != 0],
    overdispersion = FALSE,
    overdispersion_shrinkage = FALSE,
    verbose = TRUE
  )
  Mu <- fit$Mu[, match(seq(ncol(Y)), which(n != 0))]
  sign(Y - Mu) *
    sqrt(
      2 * Y * log(Y / Mu)
        + 2 * (N - Y) * log((N - Y) / (N - Mu))
    ) %>%
      replace(!is.finite(.), 0)
}
