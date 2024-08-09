devianceResiduals <- function(Y) {
  n <- colSums(Y)
  N <- matrix(n, nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE)
  fit <- glm_gp(
    Y,
    size_factors = n,
    overdispersion = FALSE,
    overdispersion_shrinkage = FALSE,
    verbose = TRUE
  )
  Mu <- fit$Mu
  sign(Y - Mu) *
    sqrt(
      2 * Y * log(Y / Mu)
        + 2 * (N - Y) * log((N - Y) / (N - Mu))
    ) %>%
      replace(!is.finite(.), 0)
}
