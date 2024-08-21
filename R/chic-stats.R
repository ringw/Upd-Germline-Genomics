chic_stats_boxplot <- function(track) {
  jitter_sd_width <- 0.12
  jitter_height <- 0.001
  tibble(
    x = droplevels(as.factor(seqnames(track))) %>%
      fct_recode(`2`="2L", `2`="2R", `3`="3L", `3`="3R") %>%
      relevel("X"),
    y = track$score
  ) %>%
    ggplot(aes(x, y)) +
    geom_boxplot(color="#3d3c3b", fill="#f0f1a1", outlier.shape = NA) +
    rasterise(
      geom_point(
        data=\(data) data %>%
          group_by(x) %>%
          slice_sample(n=5000),
        shape=".",
        position=ggforce::position_jitternormal(jitter_sd_width, jitter_height)
      ),
      dpi=200
    ) +
    scale_y_continuous(limits=c(0, 5)) +
    coord_cartesian(
      c(0.5 - 0.05, 5.5 + 0.05),
      c(0, 3),
      expand=FALSE
    ) +
    theme_bw() +
    labs(x = NULL, y = "Enrichment (vs Auto Monosome Median)")
}