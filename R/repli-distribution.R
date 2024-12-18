plot_repli_timing_byfeature <- function(repli.timing.byfeature, tss_only = FALSE, box_ci = FALSE) {
  if (!tss_only) {
    repli.timing.byfeature$occur_tss <- 1
  }
  gg <- ggplot(
    # Reverse (top to bottom) feature. Reverse again celltype!
    repli.timing.byfeature %>%
      mutate(celltype = factor(celltype, c("Somatic", "Germline"))),
    aes(feature, timing, fill = celltype)
  ) +
    facet_wrap(vars(row), ncol = 3, scales = "free") +
    geom_violin(aes(weight = occur_tss), \(df) df %>% subset(occur_tss != 0)) +
    geom_blank(
      aes(),
      tibble(
        feature = "", timing = 0, celltype = c("Germline", "Somatic"), row = factor("", c("2", "3", ""))
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

flatten_diff_rep_program_chromatin <- function(chromosomewise) {
  data <- chromosomewise %>%
    rowwise() %>%
    reframe(
      timing,
      tribble(
        ~celltype, ~value, ~H3K4, ~H3K27, ~H3K9,
        factor("Germline", c("Germline", "Somatic")),
        timing_germline,
        chromatin_germline[, "H3K4"],
        chromatin_germline[, "H3K27"],
        chromatin_germline[, "H3K9"],
        factor("Somatic", c("Germline", "Somatic")),
        timing_somatic,
        chromatin_somatic[, "H3K4"],
        chromatin_somatic[, "H3K27"],
        chromatin_somatic[, "H3K9"],
      )
    )
}

plot_diff_rep_program_chromatin <- function(chromosomewise_flat_feature) {
  values <- chromosomewise_flat_feature %>%
    group_by(timing, celltype) %>%
    reframe(
      tibble(H3K4, H3K27, H3K9) %>%
        as.matrix() %>%
        melt() %>%
        as_tibble()
    )
  g <- values %>%
    # Reverse (top to bottom) feature. Reverse again celltype!
    # mutate(celltype = factor(celltype, c("Somatic", "Germline"))) %>%
    ggplot(aes(Var2, value, fill = celltype)) +
    facet_wrap(vars(timing), ncol = 1, scales = "free") +
    geom_violin(scale = "width", width = 0.75, position = position_dodge(width = 0.9)) +
    scale_y_reverse(breaks = -1:1) +
    scale_fill_manual(values = cell_type_violin_colors) +
    coord_cartesian(c(0.4, 3.6), c(-1.5, 1.5), expand=F) +
    labs(x = NULL, y = "") +
    theme(
      aspect.ratio = 1.5,
      panel.grid.minor.x = element_blank(),
      legend.position = "none"
    )
  p <- ggplotGrob(g) %>%
    gtable_add_grob(
      list(
        textGrob(
          "L2FC",
          rot = 90,
          gp = gpar(fontfamily = "Helvetica", fontsize = 10)
        )
      ) %>%
        rep(2),
      t = c(8, 13),
      l = 3
    )
}

chromosomewise_test_t <- function(
    chromosomewise_flat_feature, facet.name, evaluate.name) {
  # Grouping variables data frame:
  faceter <- chromosomewise_flat_feature[facet.name] %>%
    subset(!duplicated(pull(., 1)))
  tibble(
    cross_join(faceter, chic.mark.data$mark %>% factor(., .) %>% tibble(mark = .)),
    do.call(
      rbind,
      mapply(
        \(facet.use, mark.use) chromosomewise_flat_feature %>%
          subset(pull(., facet.name) == facet.use) %>%
          with(
            # evaluate.name gives the factor (two levels) to test.
            split(
              pull(., mark.use),
              pull(., evaluate.name)
            )
          ) %>%
          setNames(c("x", "y")) %>%
          with(
            cbind(
              with(
                t.test(x, y),
                tibble(
                  t.diameter_500 = statistic,
                  df.diameter_500 = parameter,
                  x = estimate["mean of x"],
                  y = estimate["mean of y"],
                )
              ),
              tibble(
                sd.x = sd(x, na.rm = T),
                sd.y = sd(y, na.rm = T)
              )
            )
          ),
        get(facet.name),
        factor(mark, levels = chic.mark.data$mark),
        SIMPLIFY = FALSE
      )
    ),
    effect = x - y,
    t.statistic = t.diameter_500 / sqrt(7),
    df.statistic = df.diameter_500 / 7,
    signif = (pt(-abs(t.statistic), df.statistic) * 2 * 2 * length(t.statistic)) %>%
      cut(c(-Inf, 1e-4, 1e-3, 1e-2, 5e-2, Inf)) %>%
      `levels<-`(value = c("****", "***", "**", "*", "")),
  )
}

plot_chromosomewise_test_t <- function(data, evaluate) {
  facet <- pull(data, 1)
  eval.x <- seq(-1.5, 1.5, by = 0.05)
  if (!("celltype" %in% colnames(data))) data$celltype <- NA
  bell <- data %>%
    tibble(., facet) %>%
    rowwise() %>%
    reframe(
      mark,
      facet,
      celltype,
      dist.x = rep(eval.x, 2),
      dist.y = c(
        dnorm(eval.x, x, sd.x),
        dnorm(eval.x, y, sd.y)
      ),
      group = rep(evaluate, rep(length(eval.x), 2)),
    )
  vline <- data %>%
    tibble(., facet) %>%
    rowwise() %>%
    reframe(
      mark,
      facet,
      celltype,
      x = c(x, y),
      group = evaluate,
    )
  signif <- data %>%
    tibble(., facet) %>%
    reframe(
      mark,
      facet,
      dist.x = pmin(x, y),
      xend = pmax(x, y),
      dist.y = 0.75,
      yend = 0.75,
      signif,
    ) %>%
    subset(signif != "")
  if (is.na(data$celltype[1])) {
    bell$celltype <- bell$group
    vline$celltype <- vline$group
  }
  g <- bell %>%
    ggplot(aes(dist.x, dist.y, color = group)) +
    facet_grid(vars(facet), vars(mark)) +
    geom_line() +
    geom_segment(
      aes(x, y = -Inf, xend = x, yend = Inf),
      data = vline
    ) +
    geom_segment(
      aes(xend = xend, yend = yend, color = NULL),
      data = signif
    ) +
    geom_text(
      aes(label = signif, color = NULL),
      data = signif,
      hjust = 1.12
    ) +
    scale_color_manual(
      values = c(
        setNames(
          unlist(chic_line_track_colors, use.names = F),
          c("Germline", "Somatic")
        ),
        setNames(
          unlist(repli_level_colors[c("E", "ML")], use.names = F),
          c("GermlineEarlier", "GermlineLater")
        )
      )
    ) +
    scale_x_continuous(breaks = -1:1) +
    coord_cartesian(
      c(-1.5, 1.5),
      c(0, 1),
      expand = F
    ) +
    labs(x = "L2FC", y = "P(L2FC)") +
    theme(
      aspect.ratio = 1,
      legend.position = "none",
      panel.spacing.y = unit(10, "pt"),
      plot.margin = margin(5.5, 5.5, 5.5, 72.5),
      strip.background.y = element_blank(),
      strip.text.y = element_blank(),
    )
  p <- set_panel_size(g, w = unit(1.25, "in"), h = unit(1.25, "in"))
  p <- p %>%
    gtable_add_grob(
      list(
        textGrob(
          levels(facet)[1],
          1, 0.5,
          just = "right",
          gp = gpar(fontfamily = "Helvetica", fontsize = 10)
        ),
        textGrob(
          levels(facet)[2],
          1, 0.5,
          just = "right",
          gp = gpar(fontfamily = "Helvetica", fontsize = 10)
        )
      ),
      t = c(8, 10),
      l = 1
    )
}
