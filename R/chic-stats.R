chic_stats_boxplot <- function(track) {
  tibble(
    x = droplevels(as.factor(seqnames(track))) %>%
      fct_recode(`2`="2L", `2`="2R", `3`="3L", `3`="3R") %>%
      relevel("X"),
    y = track$score
  ) %>%
    ggplot(aes(x, y)) +
    geom_boxplot(fill="#f0f1a1") +
    scale_y_continuous(limits=c(0, 5)) +
    coord_cartesian(
      c(0.5 - 0.05, 5.5 + 0.05),
      c(0, 3),
      expand=FALSE
    ) +
    theme_bw() +
    labs(x = NULL, y = "Enrichment (vs Auto Monosome Median)")
}