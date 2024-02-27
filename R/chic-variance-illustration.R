# F-test of equality of variances of the input tracks (biological replicates).
# Are there many genomic windows where, when pairing up input tracks at the
# window, there is not homogeneity of variance? In this case,
# pct_flagged_genomic_windows will far exceed the 5% of genomic windows where
# this is expected by random chance. If pct_flagged_genomic_windows passes this
# QC test, then we can estimate the standard deviation for our t-test pooling
# the input tracks together.
mean_coverage_track_var_test <- function(track1, track2, bin_size = 100) {
  mu1 <- track1 %>% subsample_unlist_rle_list(bin_size = bin_size)
  mu2 <- track2 %>% subsample_unlist_rle_list(bin_size = bin_size)
  series1 <- track1 %>% attr("standard_deviation") %>% subsample_unlist_rle_list(bin_size = bin_size)
  series2 <- track2 %>% attr("standard_deviation") %>% subsample_unlist_rle_list(bin_size = bin_size)
  square <- \(x) x^2
  pval <- pf(
    square(series1 / series2),
    attr(track1, "n") - 1,
    attr(track2, "n") - 1
  )
  pval <- ifelse(
    pval < 0.5,
    2 * pval,
    2 - 2 * pval
  ) %>%
    replace(pmin(mu1, mu2) < 0.1, NA)
  list(pvals = pval, pct_flagged_genomic_windows = mean(pval < 0.05))
}

estimate_variance_input_tracks <- function(input_tracks, bin_size=50) {
  square <- \(x) x^2
  sd <- input_tracks %>%
    sapply(
      \(tr) tr %>%
        attr("standard_deviation") %>%
        subsample_unlist_rle_list(bin_size = bin_size)
    )
  sqrt(rowMeans(square(sd)))
}

plot_chic_anova <- function(mod_track, input_tracks, bin_size=50) {
  pair_tracks <- combn(length(input_tracks), 2, simplify=FALSE)
  pvals <- pair_tracks %>%
    sapply(
      \(inds) ecdf(mean_coverage_track_var_test(input_tracks[[inds[1]]], input_tracks[[inds[2]]], bin_size=bin_size)$pvals),
      simplify=FALSE
    )
  pair_names <- sapply(
    pair_tracks,
    \(v) v %>% sapply(
        \(idx) names(input_tracks)[idx] %>%
          str_extract("H3K[0-9]+")
      ) %>%
      replace(1, paste0("Input: ", .[1])) %>%
      as.list %>%
      append(list(sep="/")) %>%
      do.call(paste, .)
  )
  plot_data <- pvals %>% mapply(
    \(n, fn) tibble(
      comparison = n,
      F = seq(0, 1, length.out = 101),
      p = fn(F)
    ),
    pair_names, .,
    SIMPLIFY = FALSE
  ) %>% do.call(rbind, .)

  # Mod/Input f-test
  sd_est <- estimate_variance_input_tracks(input_tracks, bin_size)
  sd_mod <- mod_track %>% attr("standard_deviation") %>% subsample_unlist_rle_list(bin_size = bin_size)
  square <- \(x) x^2
  pval_est <- pf(
    square(sd_mod / sd_est),
    attr(mod_track, "n") - 1,
    sum(
      sapply(input_tracks, \(tr) attr(tr, "n") - 1)
    )
  )
  plot_data <- plot_data %>% rbind(
    tibble(
      comparison = "Mod/Inputs",
      F = seq(0, 1, length.out=101),
      p = ecdf(pval_est)(F)
    )
  )

  ggplot(
    plot_data, aes(F, p)
  ) + geom_line(aes(group=comparison, color=comparison), \(df) df %>% subset(comparison != "Mod/Inputs")) + geom_line(aes(group=comparison, color=comparison), \(df) df %>% subset(comparison == "Mod/Inputs"), linetype="dotted") + scale_color_manual(
    values = c("purple", "cyan", "orange", "magenta")
  ) + geom_line(
    data=data.frame(F=c(0,1), p=c(0,1)),
    linetype="dashed"
  ) + geom_line(
    data=data.frame(
      F=c(0.025,0.025,0, 0.975,0.975,1),
      p=c(0.025,0,0.025, 1,0.975,0.975),
      group=rep(c("b", "t"), each=3)
    ),
    aes(group=group),
    color="blue"
  ) + theme_bw() + theme(
  ) + coord_cartesian(
    expand=FALSE
  )
}