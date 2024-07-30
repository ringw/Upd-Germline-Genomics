# Logistic Distribution - on tanh scale - link functions.
# These are the Logistic Distribution CDF. The tanh/arctanh functions would
# actually rescale the predicted values before link function as well as
# afterwards. However, the range of our responses will actually follow the
# range of tanh - scale Logistic probabilities to be on the range (-1,1), so
# that they look good in IGV.
qlogistanh <- \(q) qlogis(0.5*q + 0.5)
plogistanh <- \(x) (2*plogis(x) - 1)

repli_logistic_beta_to_tanh <- function(betas) {
  -2*betas + 1
}

repli_logistic_link_apply_scale <- function(granges, num_fractions) {
  score <- granges$score %>% qlogistanh
  cntr <- scale(score[which(is.finite(score))])
  score[which(is.finite(score))] <- as.numeric(cntr)
  score <- score * (
    qlogistanh((num_fractions - 1) / num_fractions)
    /
    qnorm(0.975)
  )
  score <- score %>% plogistanh
  granges$score <- score
  granges@metadata[["scaled:center"]] <- attr(cntr, "scaled:center")
  granges@metadata[["scaled:scale"]] <- attr(cntr, "scaled:scale")
  granges
}