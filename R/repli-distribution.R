plot_repli_timing_byfeature <- function(repli.timing.byfeature, tss_only = FALSE, box_ci = FALSE) {
  if (!tss_only) {
    repli.timing.byfeature$occur_tss <- 1
  }
  gg <- ggplot(
    # Reverse (top to bottom) feature. Reverse again celltype!
    repli.timing.byfeature %>%
      mutate(celltype = factor(celltype, c("Somatic", "Germline"))),
    aes(feature, timing, fill=celltype)
  ) +
    facet_wrap(vars(row), ncol=3, scales="free") +
    geom_violin(aes(weight=occur_tss), \(df) df %>% subset(occur_tss != 0)) +
    geom_blank(
      aes(),
      tibble(
        feature="", timing=0, celltype=c("Germline", "Somatic"), row=factor("", c("2", "3", ""))
      )
    ) +
    scale_x_discrete(limits = rev) +
    scale_y_reverse() +
    scale_fill_manual(values = cell_type_violin_colors) +
    labs(x = "Genomic Range") +
    coord_flip() +
    theme(
      aspect.ratio = 1.5,
      panel.grid.major.x = element_blank(),
      legend.position = "none"
    )
  if (box_ci) {
    data <- ggplot_build(gg)$data[[1]]
    # 75% CI information. We persist the density results in ggplot_build data, as
    # well as the final coordinate space x and y, so we need to match the quantile
    # to plot to the cumsum of % density and then get the y value for this
    # particular variable value.
    box_data <- data %>%
      group_by(group) %>%
      summarise(
        row = (levels(repli.timing.byfeature$row) %>% factor(., .))[as.numeric(PANEL[1])],
        xmin = xmin[1] + (xmax - xmin)[1] * (1 - width[1]) / 2,
        xmax = xmax[1] - (xmax - xmin)[1] * (1 - width[1]) / 2,
        ymin = y[
          findInterval(
            0.25,
            cumsum(density) / sum(density)
          )
        ],
        ymax = y[
          findInterval(
            0.75,
            cumsum(density) / sum(density)
          )
        ]
      )
    median_data <- data %>%
      group_by(group) %>%
      summarise(
        row = (levels(repli.timing.byfeature$row) %>% factor(., .))[as.numeric(PANEL[1])],
        x = xmin[1] + (xmax - xmin)[1] * (1 - width[1]) / 2,
        x = xmin[1],
        y = y[
          findInterval(
            0.5,
            cumsum(density) / sum(density)
          )
        ],
        xend = xmax[1] - (xmax - xmin)[1] * (1 - width[1]) / 2,
        yend = y
      )
    gg +
      geom_rect(
        aes(x = NULL, y = NULL, xmin = xmin, xmax = xmax, ymin = -ymin, ymax = -ymax, fill = NULL),
        box_data,
        color = "black",
        fill = "transparent"
      ) +
      geom_segment(
        aes(x, -y, xend = xend, yend = -yend, fill = NULL),
        median_data,
        color = "black",
        linewidth = 1.25
      )
  } else {
    gg
  }
}
