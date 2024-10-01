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

dm_barplot <- function(
  repli.dm,
  chromosome_arms_diameter_1000,
  grouping = as.factor(seqnames(chromosome_arms_diameter_1000)) %>%
    factor(levels = head(levels(.), 7))
) {
  data <- sapply(
    repli.dm,
    \(lst) lst$par
  ) %>%
    t()
  data <- data / rowSums(data)
  chrs <- t(
    sapply(
      split(as.data.frame(data), grouping),
      \(df) colSums(as.matrix(df))
    )
  ) %>%
    `/`(rowSums(.))
  ggplot(
    melt(chrs) %>%
      mutate(
        chr = Var1,
        fraction = factor(Var2) %>%
          recode(V1="E", V2="EM", V3="ML", V4="L") %>%
          fct_relevel(c("L", "ML", "EM", "E"))
      ),
    aes(Var1, value, fill=fraction)
  ) +
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values = unlist(repli_level_colors)) +
    scale_x_discrete(limits = rev) +
    scale_y_continuous(labels = percent) +
    labs(
      x = "dm6 chromosome",
      y = "Percent in Fraction"
    ) +
    coord_flip(expand = FALSE) +
    theme(
      aspect.ratio = 1/3 * nrow(chrs) / 7,
      legend.position = "none"
    )
}