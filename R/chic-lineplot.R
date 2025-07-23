chic_line_track_colors <- list(
  germline = "#5DC663",
  somatic = "#D845E8"
)

chic_lineplot_limit_data_TSS <- tribble(
  ~mark_name, ~value,
  "H3K4me3", 0.67,
  "H3K4me3", 4.88,
  "H3K27me3", 0.57,
  "H3K27me3", 3.73,
  "H3K9me3", 0.5,
  "H3K9me3", 1.35,
) %>%
  tibble(mark = mark_name %>% factor(., unique(.)))

quant_quartile_factor <- function(v, q1_threshold) {
  chop_lower = (v < q1_threshold) %>% replace_na(TRUE) %>% which
  vals <- v[-chop_lower]
  thirds <- cut(vals, c(0, quantile(vals, c(1/3, 2/3)), Inf))
  levels(thirds) <- c("Q2", "Q3", "Q4")
  fct <- factor(rep("", length(v)), c("Q1", "Q2", "Q3", "Q4"))
  fct[chop_lower] <- "Q1"
  fct[-chop_lower] <- as.character(thirds)
  setNames(fct, names(v))
}

chic_heatmap_facet_genes <- function(
  enrichment_mat,
  facet_genes
) {
  if (!length(enrichment_mat))
    enrichment_mat <- matrix(
      nrow = length(facet_genes$gene),
      ncol = 0,
      dimnames = list(facet_genes$gene, character(0))
    )
  df <- facet_genes %>%
    group_by(subset(facet_genes, select=-gene)) %>%
    reframe(
      n = length(gene),
      enrichment_mat[gene,, drop=F] %>%
        colMeans(na.rm=T) %>%
        enframe("pos")
    )
  if (length(attr(enrichment_mat, "x")))
    as_tibble(cbind(df, x = attr(enrichment_mat, "x")))
  else
    df %>% mutate(pos = factor(pos, unique(pos)))
}

# Analysis with TSS profile as x-axis, with ChIC tracks.
chic_average_profile_limits <- c(0.25, 5)
chic_average_breaks <- c(1/2, 1, 2, 3, 4, 8)
chic_average_minor_breaks <- c(1/sqrt(2), sqrt(2), 4*sqrt(2))

# TSS plot. Rows come from var named "facet". Cols come from "mark".
# Grouping within the panel is always called "genes" because it normally
# consists of a group-by then mean (using disjoint sets).
# Line plot is called "l2FC" although it is actually Fold Change and log
# transform is applied to the y-axis in ggplot.
chic_plot_average_profiles_facet_grid <- function(
  facet_data, legend_title, quartile_colors, linewidth=c(0.33, 0.6, 0.75, 1),
  faceter = facet_grid(rows = vars(mark), cols = vars(facet)),
  x_intercept_red = NA,
  chic_average_profile_limits = get("chic_average_profile_limits", envir=globalenv()),
  chic_average_breaks = get("chic_average_breaks", envir=globalenv())
) {
  break_labels <- tibble(
    pos = seq(head(levels(facet_data$pos), 1), tail(levels(facet_data$pos), 1)),
    label = levels(facet_data$pos)
  )
  # Get the facet factor levels to use from a global. However, if these are not
  # the factor levels being used, then just convert the data to factor and use
  # alphabetical order instead.
  mark_names <- str_glue("{chic.mark.data$mark}me3")
  if (any(mark_names %in% facet_data$mark))
    mark <- factor(mark_names, mark_names, ordered=TRUE)
  else
    mark <- factor(unique(facet_data$mark))
  facet_data$pos.continuous <- break_labels$pos[match(facet_data$pos, break_labels$label)]
  facet_data %>% ggplot(
    aes(x=pos.continuous, y=value, color=genes, linewidth=genes, group=genes)
  ) + geom_line(
    data = tribble(~pos.continuous, ~value, -Inf, x_intercept_red, Inf, x_intercept_red) %>%
      cross_join(
        tibble(
          genes = NA,
          mark = mark
        )
      ) %>%
      cross_join(
        tibble(
          facet = factor(levels(facet_data$facet), levels(facet_data$facet), ordered=TRUE)
        )
      ),
    color = "darkred",
    linewidth = 0.25
  ) + geom_line() + faceter + scale_color_manual(
    values = quartile_colors,
    guide = guide_legend(title = legend_title)
  ) + scale_linewidth_manual(
    values = linewidth,
    guide = guide_legend(title = legend_title)
  ) + scale_x_continuous(
    labels = \(n) break_labels$label[match(n, break_labels$pos)],
    expand = c(0, 0)
  ) + scale_y_continuous(
    trans = "log",
    labels = \(v) round(log(v) / log(2), 1),
    breaks = chic_average_breaks,
    minor_breaks = chic_average_minor_breaks
  ) + labs(
    x = "bp (from TSS)", y = "log2(mean(mark/input))"
  ) + theme(
    aspect.ratio = 1
  )
}

chic_plot_average_profiles_facet_grid2 <- function(
  facet_data, legend_title, quartile_colors, linewidth=c(0.33, 0.6, 0.75, 1),
  faceter = facet_grid(rows = vars(mark), cols = vars(facet)),
  x_intercept_red = NA,
  chic_average_profile_limits = get("chic_average_profile_limits", envir=globalenv()),
  chic_average_breaks = get("chic_average_breaks", envir=globalenv())
) {
  break_labels <- tibble(
    pos = seq(head(levels(facet_data$pos), 1), tail(levels(facet_data$pos), 1)),
    label = levels(facet_data$pos)
  )
  # Get the facet factor levels to use from a global. However, if these are not
  # the factor levels being used, then just convert the data to factor and use
  # alphabetical order instead.
  mark_names <- str_glue("{chic.mark.data$mark}me3")
  if (any(mark_names %in% facet_data$mark))
    mark <- factor(mark_names, mark_names, ordered=TRUE)
  else
    mark <- factor(unique(facet_data$mark))
  facet_data$pos.continuous <- break_labels$pos[match(facet_data$pos, break_labels$label)]
  facet_data %>% ggplot(
    aes(x=pos.continuous, y=value, color=genes, linewidth=genes, group=genes)
  ) +
  geom_point(
    data = tibble(
      pos.continuous = 0,
      mark = rep(mark, each=2),
      value = c(1/2, 4, 1/2, 4, 1/2, 2.0001),
      genes = levels(facet_data$genes)[1],
    ),
    color = "transparent",
    fill = "transparent"
  ) +
  geom_line(
    data = tribble(~pos.continuous, ~value, -Inf, x_intercept_red, Inf, x_intercept_red) %>%
      cross_join(
        tibble(
          genes = NA,
          mark = mark
        )
      ) %>%
      cross_join(
        tibble(
          facet = factor(levels(facet_data$facet), levels(facet_data$facet), ordered=TRUE)
        )
      ),
    color = "darkred",
    linewidth = 0.25
  ) + geom_line() + faceter + scale_color_manual(
    values = quartile_colors,
    guide = guide_legend(title = legend_title)
  ) + scale_linewidth_manual(
    values = linewidth,
    guide = guide_legend(title = legend_title)
  ) + scale_x_continuous(
    labels = \(n) break_labels$label[match(n, break_labels$pos)],
    expand = c(0, 0)
  ) + scale_y_continuous(
    trans = "log",
    labels = \(v) round(log(v) / log(2), 1),
    breaks = chic_average_breaks,
    minor_breaks = chic_average_minor_breaks,
    expand = c(0, 0)
  ) + labs(
    x = "bp (from TSS)", y = "log2(mean(mark/input))"
  ) + theme(
    aspect.ratio = 1
  )
}

# H3 plot. The absolute panel size here is an arbitrary number, from the RNAseq
# Quartile plot being 4 inches high, because we didn't set more specific
# absolute units on the plot.
chic_plot_h3_enrichment <- function(facet_data, limits = c(0, 5.25), ...) {
  y_axis_h3 <- scale_y_continuous(
    "Enrichment (vs Auto Monosome Median)",
    limits = limits,
    expand = c(0, 0)
  )
  gg <- chic_plot_average_profiles_facet_grid(facet_data, ...) + y_axis_h3
  gr <- as_grob(gg)
  fig3PlotSize <- unit(2.471, "in")
  gr$heights[
    grid::unitType(gr$heights) == "null" &
      as.numeric(gr$heights) == 1
  ] <- fig3PlotSize
  gr$widths[
    grid::unitType(gr$widths) == "null" &
      as.numeric(gr$widths) == 1
  ] <- fig3PlotSize
  plot_grid(gr)
}

# Now an average profiles facet grid with "H3K4" and "H3K27" factors.
plot_valency_k4_k27 <- function(facet_data, plot_color) {
  chic_valency_profile_limits <- c(0.44, 8.5) %>%
    log %>%
    `+`(c(-0.02, 0.02) * diff(.)) %>%
    exp
  break_labels <- tibble(
    pos = seq(head(levels(facet_data$pos), 1), tail(levels(facet_data$pos), 1)),
    label = levels(facet_data$pos)
  )
  facet_data$pos.continuous <- break_labels$pos[match(facet_data$pos, break_labels$label)]
  facet_data$genes <- interaction(facet_data$H3K4, facet_data$H3K27) %>%
    `levels<-`(value = c("Neither", "H3K4me3 H3K27~", "H3K4~ H3K27me3", "Bivalent")) %>%
    fct_relevel(c("H3K4me3 H3K27~", "Bivalent", "Neither", "H3K4~ H3K27me3"))
  facet_data %>% ggplot(
    aes(x=pos.continuous, y=value, group=mark, linetype=mark)
  ) + facet_wrap(vars(genes), scales="free") +
    geom_line(color=plot_color) +
    scale_linetype_manual(
      values = c("solid", "dashed")
    ) + scale_x_continuous(
      labels = \(n) break_labels$label[match(n, break_labels$pos)]
    ) + scale_y_continuous(
      trans = "log",
      labels = \(v) round(log(v) / log(2), 1),
      limits = chic_valency_profile_limits,
      breaks = chic_average_breaks,
      minor_breaks = chic_average_minor_breaks,
      expand = c(0, 0)
    ) + coord_cartesian(expand=FALSE) + labs(
      x = "bp (from TSS)", y = "log2(mean(mark/input))"
    ) + theme(
      aspect.ratio = 1
    )
}

chic_plot_paneled_profiles_facet_grid <- function(
  facet_data, legend_title, quartile_colors, linewidth=c(0.33, 0.6, 0.75, 1),
  faceter = facet_grid(rows = vars(mark), cols = vars(facet)),
  x_intercept_red = 1,
  chic_average_profile_limits = get("chic_average_profile_limits", envir=globalenv()),
  chic_average_breaks = get("chic_average_breaks", envir=globalenv())
) {
  tss_left <- min(facet_data$x)
  # Infer flanking of TSS and TES.
  bp_flanking <- head(which(facet_data$pos == "TSS") - 1, 1)
  tss_right <- facet_data$x[which(facet_data$pos == bp_flanking)[1]]
  tes_left <- facet_data$x[which(facet_data$pos == -bp_flanking)[2]]
  tes_right <- max(facet_data$x)
  label_data <- tribble(
    ~x, ~value,
    tss_left, as.character(-bp_flanking),
    mean(c(tss_left, tss_right)), "TSS",
    tss_right, as.character(bp_flanking),
    mean(c(tss_right, tes_left)), "50%",
    tes_left, as.character(-bp_flanking),
    mean(c(tes_left, tes_right)), "TES",
    tes_right, as.character(bp_flanking)
  )
  minor_breaks <- unique(facet_data$x[facet_data$pos %in% c("25%", "75%")])
  facet_data %>% ggplot(
      aes(x, y = value, color=genes, linewidth=genes, group=genes)
    ) + geom_tile(
      # Colored tile (TSS and TES)
      data = tribble(
        ~x, ~value, ~width,
        mean(c(tss_left, tss_right)),
        1,
        2 * bp_flanking + 1,
        mean(c(tes_left, tes_right)),
        1,
        2 * bp_flanking + 1
      ) %>%
        tibble(
          height = Inf,
          genes = NA
        ),
      aes(width=width, height=height, color=NULL, linewidth=NULL, genes=NULL),
      color = "transparent",
      fill = "#fffccb"
    ) + faceter + scale_color_manual(
      values = quartile_colors,
      guide = guide_legend(title = legend_title, override.aes = list(fill = "transparent"))
    ) + scale_linewidth_manual(
      values = linewidth,
      guide = guide_legend(title = legend_title)
    ) + coord_cartesian(
      c(tss_left, tes_right),
      chic_average_profile_limits,
      expand=FALSE
    ) + scale_y_continuous(
      trans = "log",
      labels = \(v) round(log(v) / log(2), 1),
      limits = chic_average_profile_limits,
      breaks = chic_average_breaks,
      minor_breaks = chic_average_minor_breaks,
      expand = c(0, 0)
    ) + theme(
      panel.background = element_rect(fill = NA),
      panel.ontop = TRUE
    ) + geom_line(
      # Dark red line at the origin (enrichment of 1)
      data = cross_join(
        tribble(~x, ~value, -Inf, 1.01, Inf, 1.01),
        tibble(
          genes = NA,
          chic_mark = factor(chic.mark.data$mark, chic.mark.data$mark, ordered=TRUE)
        )
      ),
      color = "darkred",
      linewidth = 0.33
    ) + geom_line() + labs(
      x = "base pairs", y = "log2(mean(mark/input))"
    ) + scale_x_continuous(
      breaks = label_data$x,
      labels = label_data$value,
      minor_breaks = minor_breaks
    ) + theme(
      panel.margin = unit(25, "pt"),
      aspect.ratio = 1
    )
}

chic_panel_create_grob <- function(
  facet_data,
  facet_names,
  column_width,
  row_height
) {
  facet1 <- tribble(
    ~level1, ~row,
    1, 3,
    2, 5
  )
  facet2 <- tribble(
    ~level2, ~column,
    1, 2,
    2, 4
  )
  facets <- bind_rows(
    setNames(
      list(facet1, facet2),
      names(facet_names)
    ),
    .id = "name"
  )
  gtable(
    unit(c(0.1, column_width, column_width), "in"),
    unit(c(0.25, 0.25, row_height, 0.25, row_height), "in")
  ) %>%
    gtable_add_grob(
      textGrob(
        str_glue(names(facet_names)[1], levels(pull(facet_data, facet_names[1]))[2]),
        rot=90
      ),
      3, 1
    ) %>%
    gtable_add_grob(
      textGrob(
        str_glue(names(facet_names)[1], levels(pull(facet_data, facet_names[1]))[1]),
        rot=90
      ),
      5, 1
    ) %>%
    gtable_add_grob(
      textGrob(
        str_glue(names(facet_names)[2], levels(pull(facet_data, facet_names[2]))[1])
      ),
      1, 2
    ) %>%
    gtable_add_grob(
      textGrob(
        str_glue(names(facet_names)[2], levels(pull(facet_data, facet_names[2]))[2])
      ),
      1, 3
    ) %>%
    gtable_add_grob(
      textGrob(
        apply(
          filter(
            facet_data,
            pos == "TSS",
            mark == "H3K4me3",
            as.numeric(pull(facet_data, facet_names[1])) == 2,
            as.numeric(pull(facet_data, facet_names[2])) == 1
          ) %>%
            arrange(genes),
          1,
          \(v) str_glue("{v['genes']}: {v['n']}")
        ) %>%
          as.list %>%
          do.call(paste, .)
      ),
      2, 2
    ) %>%
    gtable_add_grob(
      textGrob(
        apply(
          filter(
            facet_data,
            pos == "TSS",
            mark == "H3K4me3",
            as.numeric(pull(facet_data, facet_names[1])) == 2,
            as.numeric(pull(facet_data, facet_names[2])) == 2
          ) %>%
            arrange(genes),
          1,
          \(v) str_glue("{v['genes']}: {v['n']}")
        ) %>%
          as.list %>%
          do.call(paste, .)
      ),
      2, 3
    ) %>%
    gtable_add_grob(
      textGrob(
        apply(
          filter(
            facet_data,
            pos == "TSS",
            mark == "H3K4me3",
            as.numeric(pull(facet_data, facet_names[1])) == 1,
            as.numeric(pull(facet_data, facet_names[2])) == 1
          ) %>%
            arrange(genes),
          1,
          \(v) str_glue("{v['genes']}: {v['n']}")
        ) %>%
          as.list %>%
          do.call(paste, .)
      ),
      4, 2
    ) %>%
    gtable_add_grob(
      textGrob(
        apply(
          filter(
            facet_data,
            pos == "TSS",
            mark == "H3K4me3",
            as.numeric(pull(facet_data, facet_names[1])) == 1,
            as.numeric(pull(facet_data, facet_names[2])) == 2
          ) %>%
            arrange(genes),
          1,
          \(v) str_glue("{v['genes']}: {v['n']}")
        ) %>%
          as.list %>%
          do.call(paste, .)
      ),
      4, 3
    )
}

chic_panel_gtable_binary <- function(
  facet_data,
  facet_fn,
  facet_names,
  column_width,
  row_height,
  quartile_color,
  linewidth,
  chic_average_profile_limits
) {
  facet_data <- facet_data %>% subset(activity == "active")
  chic_panel_create_grob(facet_data, facet_names, column_width, row_height) %>%
    gtable_add_grob(
      as_grob(
        facet_fn(
          facet_data %>%
            filter(as.numeric(pull(., facet_names[1])) == 2, as.numeric(pull(., facet_names[2])) == 1),
          faceter=facet_wrap(vars(mark), nrow=1),
          quartile_color=quartile_color,
          linewidth=linewidth,
          legend_title=NULL,
          chic_average_profile_limits = chic_average_profile_limits
        ) +
          theme(legend.position = "none")
      ),
      3,
      2
    ) %>%
    gtable_add_grob(
      as_grob(
        facet_fn(
          facet_data %>%
            filter(as.numeric(pull(., facet_names[1])) == 2, as.numeric(pull(., facet_names[2])) == 2),
          faceter=facet_wrap(vars(mark), nrow=1),
          quartile_color=quartile_color,
          linewidth=linewidth,
          legend_title=NULL,
          chic_average_profile_limits = chic_average_profile_limits
        ) +
          theme(legend.position = "none")
      ),
      3,
      3
    ) %>%
    gtable_add_grob(
      as_grob(
        facet_fn(
          facet_data %>%
            filter(as.numeric(pull(., facet_names[1])) == 1, as.numeric(pull(., facet_names[2])) == 1),
          faceter=facet_wrap(vars(mark), nrow=1),
          quartile_color=quartile_color,
          linewidth=linewidth,
          legend_title=NULL,
          chic_average_profile_limits = chic_average_profile_limits
        ) +
          theme(legend.position = "none")
      ),
      5,
      2
    ) %>%
    gtable_add_grob(
      as_grob(
        facet_fn(
          facet_data %>%
            filter(as.numeric(pull(., facet_names[1])) == 1, as.numeric(pull(., facet_names[2])) == 2),
          faceter=facet_wrap(vars(mark), nrow=1),
          quartile_color=quartile_color,
          linewidth=linewidth,
          legend_title=NULL,
          chic_average_profile_limits = chic_average_profile_limits
        ) +
          theme(legend.position = "none")
      ),
      5,
      3
    )
}

facet_diff_replication_program <- function(data) {
  p <- data %>%
    tibble(
      genes = interaction(activity, celltype),
      row = relevel(timing, "GermlineEarlier"),
      newfacet = interaction(mark, row) %>%
        `levels<-`(
          levels(.) %>%
            replace(1:3, c("H3K4me3", "H3K27me3", "H3K9me3"))
        )
    ) %>%
    chic_plot_average_profiles_facet_grid(
      NULL,
      c(muted(classification_colors_fig4$germline, l=70), classification_colors_fig4$germline, muted(classification_colors_fig4$somatic, l=70), classification_colors_fig4$somatic),
      linewidth = c(0.33, 0.66, 0.33, 0.66),
      faceter = facet_wrap(vars(newfacet), scales="free"),
      x_intercept = NA
    ) +
    geom_blank(
      aes(x=0, y=value, color=NULL, linewidth=NULL, group=NULL),
      \(data) tibble(
        dplyr::slice(
          chic_lineplot_limit_data_TSS,
          rep(1:6, 3)
        ),
        newfacet = levels(data$newfacet) %>%
          rep(each = 2) %>%
          factor(., unique(.)),
      )
    ) +
    theme(
      legend.position = "none",
      panel.spacing.y = unit(29.5, "pt"),
      plot.margin = margin(24, 5.5, 5.5, 5.5, "pt"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank(),
    )
  g <- ggplotGrob(p)
  # Shrink away the facet wraps that we are not interested in
  g$heights[12] <- g$heights[17] <- unit(0, "cm")
  g %>%
    gtable_add_grob(
      list(
        textGrob(
          "Germline Earlier",
          0, 0,
          just = "left",
          vjust = -0.5,
          gp = gpar(fontfamily = "Helvetica", fontsize = 12)
        ),
        textGrob(
          "Germline Later",
          0, 0,
          just = "left",
          vjust = -0.5,
          gp = gpar(fontfamily = "Helvetica", fontsize = 12)
        ),
        textGrob(
          "n.s.",
          0, 0,
          just = "left",
          vjust = -0.5,
          gp = gpar(fontfamily = "Helvetica", fontsize = 12)
        )
      ),
      t = c(1, 10, 15),
      l = 5
    )
}
