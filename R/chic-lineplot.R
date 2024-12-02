chic_line_track_colors <- list(
  germline = "#5DC663",
  somatic = "#D845E8"
)

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
    labels = \(n) break_labels$label[match(n, break_labels$pos)]
  ) + scale_y_continuous(
    trans = "log",
    labels = \(v) round(log(v) / log(2), 1),
    limits = chic_average_profile_limits,
    breaks = chic_average_breaks,
    minor_breaks = chic_average_minor_breaks,
    expand = c(0, 0)
  ) + coord_cartesian(expand=FALSE) + labs(
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
    expand = FALSE
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

chic_panel_fpkm_third_density <- function(
  peakcalling, Upd_cpm, celltype
) {
  data <- log(Upd_cpm) / log(10)
  inter_cutoffs <- fpkm_quartile_factor_make_cutoffs(data)
  facet_data_cols <- paste0("H3K", c(4, 27), "_", str_to_title(celltype))
  pull_data <- chic.gene.enrichment %>%
    subset(select = facet_data_cols) %>%
    apply(
      2,
      \(v) factor(!is.na(v) & v < 0.001) %>%
        `levels<-`(value = c("~", "+")),
      simplify=F
    ) %>%
    as_tibble() %>%
    tibble(
      symbol = chic.gene.enrichment$symbol,
      .
    ) %>%
    subset(Upd_cpm[, celltype] >= 5) %>%
    group_by(subset(., select=c(2, 3)))
  facet_data <- tally(pull_data) %>%
    tibble(pos = "TSS", mark = "H3K4", genes = str_to_title(celltype))
  facet_names <- setNames(colnames(facet_data)[1:2], c("K4me3", "K27me3"))
  celltype_x_coordinate <- c(germline=1, somatic=2)
  custom_fpkm_third_density <- function(data) fpkm_third_density(
    data,
    inter_cutoffs = inter_cutoffs
  ) +
    coord_cartesian(celltype_x_coordinate[celltype] + c(0.5, 1.5), c(0, 4.25))
  chic_panel_create_grob(facet_data, facet_names, 4, 2) %>%
    gtable_add_grob(
      as_grob(
        data[
          pull_data$symbol[group_data(pull_data)$`.rows`[[3]]],
        ] %>% custom_fpkm_third_density()
      ),
      3,
      2
    ) %>%
    gtable_add_grob(
      as_grob(
        data[
          pull_data$symbol[group_data(pull_data)$`.rows`[[4]]],
        ] %>% custom_fpkm_third_density()
      ),
      3,
      3
    ) %>%
    gtable_add_grob(
      as_grob(
        data[
          pull_data$symbol[group_data(pull_data)$`.rows`[[1]]],
        ] %>% custom_fpkm_third_density()
      ),
      5,
      2
    ) %>%
    gtable_add_grob(
      as_grob(
        data[
          pull_data$symbol[group_data(pull_data)$`.rows`[[2]]],
        ] %>% custom_fpkm_third_density()
      ),
      5,
      3
    )
}

chic_average_gene_list_profiles <- function(
  gene_list,
  bg_gene_list,
  chic_path,
  metafeatures_path,
  chic_driver,
  track_color = "goldenrod"
) {
  metafeatures <- read.csv(metafeatures_path, row.names = 1)
  metafeatures_on <- metafeatures[gene_list, ] %>%
    subset(chr %in% names(chr.lengths))
  metafeatures_off <- metafeatures[bg_gene_list, ] %>%
    subset(chr %in% names(chr.lengths))

  before <- 500
  after <- 1500
  facet_data <- sapply(
    list(on=metafeatures_on, off=metafeatures_off),
    \(metafeatures) {
      mark_tracks <- sapply(
        chic.mark.data$mark,
        \(mark) flybase_big_matrix(
          metafeatures %>% subset(!is.na(chr) & !is.na(start) & !is.na(end)) %>% subset(!duplicated(flybase)),
          paste0(chic_path, "/", chic_driver, "_", mark, ".new.FE.bw"),
          before = before,
          after = after
        ) %>%
          subset(rowAlls(. > 0.001)) %>%
          # log %>%
          # `/`(log(2)) %>%
          colMeans,
        simplify=FALSE
      )
      facet_data <- mark_tracks %>%
        sapply(
          \(v) data.frame(pos=seq(-before, after-1), value=v),
          simplify=F) %>%
        bind_rows(.id = "marks") %>%
        mutate(marks = marks %>% factor(chic.mark.data$mark))
    },
    simplify=FALSE
  ) %>%
    bind_rows(.id = "group")
  facet_data$group <- facet_data$group %>% factor(c("on", "off"))
 
  facet_data %>% ggplot(
    aes(x=pos, y=l2FC)
  ) + geom_line(
    data = tribble(~pos, ~l2FC, -Inf, 1, Inf, 1),
    color = "darkred",
    linewidth = 0.5
  ) + geom_line(
    # Implement the actual ChIC profile track.
    aes(group = group, color = group, linewidth = group)
  ) + facet_wrap(
    vars(marks)
  ) + scale_color_manual(
    values = c(track_color, muted(track_color))
  ) + scale_linewidth_manual(
    values = c(1, 0.25)
  ) + scale_y_continuous(
    trans = "log",
    labels = \(v) round(log(v) / log(2), 1),
    limits = chic_average_profile_limits,
    breaks = chic_average_breaks,
    minor_breaks = chic_average_minor_breaks,
    expand = c(0, 0)
  ) + coord_cartesian(expand = FALSE) + labs(
    x = "bp (from TSS)", y = "log2(mean(mark/input))"
  ) + theme(
    aspect.ratio = 1
  )
}

chic_multi_genes_tss_profile <- function(
  facet_data,
  color,
  linewidth
) {
  before <- 500
  after <- 1500
 
  facet_data %>% ggplot(
    aes(x=pos, y=l2FC)
  ) + geom_line(
    data = tribble(~pos, ~l2FC, -Inf, 1, Inf, 1),
    color = "darkred",
    linewidth = 0.5
  ) + geom_line(
    # Implement the actual ChIC profile track.
    aes(group = genes, color = genes, linewidth = genes)
  ) + facet_wrap(
    vars(marks)
  ) + scale_color_manual(
    values = color
  ) + scale_linewidth_manual(
    values = linewidth
  ) + scale_y_continuous(
    trans = "log",
    labels = \(v) round(log(v) / log(2), 1),
    limits = chic_average_profile_limits,
    breaks = chic_average_breaks,
    minor_breaks = chic_average_minor_breaks,
    expand = c(0, 0)
  ) + coord_cartesian(expand = FALSE) + labs(
    x = "bp (from TSS)", y = "log2(mean(mark/input))"
  ) + theme(
    aspect.ratio = 1
  )
}

chic_average_gene_list_gene_body <- function(
  gene_list,
  chic_bw,
  metafeatures_path,
  after_tss = 1000,
  before_tes = 1000,
  num_points_to_sample = 100
) {
  metafeatures <- read.csv(metafeatures_path, row.names = 1)
  metafeatures <- metafeatures[gene_list, ] %>%
    subset(chr %in% names(chr.lengths))
  coverage <- import(chic_bw, "bigwig") %>% coverage(weight="score")

  metafeatures <- metafeatures %>%
    rownames_to_column %>%
    mutate(
      display_start = ifelse(
        strand == "+",
        start + after_tss,
        end - after_tss
      ),
      display_end = ifelse(
        strand == "+",
        end - before_tes,
        start + before_tes
      )
    )
  split_features <- metafeatures %>% split(.$transcript.length > (after_tss + before_tes))

  # Genes to downsample a variable-sized window from the chic track.
  display_genes <- split_features$`TRUE` %>%
    rowwise %>%
    mutate(
      values = seq(
        display_start, display_end, length.out = num_points_to_sample
      ) %>%
        round %>%
        list
    )
  display_gene_tracks <- display_genes %>%
    split(.$chr) %>%
    mapply(
      \(df, track) df$values %>%
        do.call(c, .) %>%
        `[`(track, .) %>%
        matrix(nrow = nrow(df), byrow = TRUE, dimnames = list(df$rowname, NULL)),
      .,
      coverage[names(.)],
      SIMPLIFY=FALSE
    ) %>%
    do.call(rbind, .)

  average_gene_tracks <- split_features$`FALSE` %>%
    split(.$chr) %>%
    mapply(
      \(df, track) c(df$display_start, df$display_end) %>%
        `[`(track, .) %>%
        matrix(nrow = nrow(df), ncol = 2) %>%
        rowMeans %>%
        rep(times = num_points_to_sample) %>%
        matrix(nrow = nrow(df), ncol = num_points_to_sample, dimnames = list(df$rowname, NULL)),
      .,
      coverage[names(.)],
      SIMPLIFY=FALSE
    ) %>%
    do.call(rbind, .)

  rbind(
    display_gene_tracks,
    average_gene_tracks
  )[metafeatures$rowname, ]
}

pull_chic_average_gene_list_paneled_data <- function(
  gene_list,
  chic_bw,
  metafeatures_path,
  num_tss,
  num_tes
) {
  metafeatures <- read.csv(metafeatures_path, row.names = 1)
  metafeatures <- metafeatures[gene_list, ] %>%
    subset(chr %in% names(chr.lengths))

  tss_data <- flybase_big_matrix(
    metafeatures,
    chic_bw,
    "TSS",
    before = num_tss,
    after = num_tss
  )
  tes_data <- flybase_big_matrix(
    metafeatures,
    chic_bw,
    "TES",
    before = num_tes,
    after = num_tes
  )
  inter_data = chic_average_gene_list_gene_body(
    gene_list,
    chic_bw,
    metafeatures_path,
    num_points_to_sample = 99
  )
  result <- tibble(
    labels = c(
      (-num_tss):(-1),
      "TSS",
      paste0("+", 1:(num_tss - 1)),
      paste0(1:99, "%"),
      (-num_tes):(-1),
      "TES",
      paste0("+", 1:(num_tes - 1))
    ),
    track_value = c(
        colMeans(tss_data),
        colMeans(inter_data),
      colMeans(tes_data)
    )
  )
  colnames(result)[2] <- basename(chic_bw) %>% str_replace("[.].*", "")
  result
}

pull_chic_average_gene_list_paneled_profiles_data <- function(
  gene_lists,
  chic_path,
  metafeatures_path,
  chic_driver,
  tss_size,
  inter_size,
  tes_size,
  num_tss = 500,
  num_tes = 500
) {
  num_inter <- 99
  panel_data <- tibble(
    chic_driver = chic_driver,
    chic_mark = chic.mark.data$mark,
    filename = paste0(chic_path, "/", chic_driver, "_", chic_mark, ".new.FE.bw")
  ) %>%
    cross_join(
      tibble(
        group = names(gene_lists) %>%
          coalesce(as.character(seq_along(gene_lists))) %>%
          factor(., ., ordered=TRUE),
        gene_list = gene_lists
      )
    ) %>%
    rowwise %>%
    mutate(
      track_data = list(
        pull_chic_average_gene_list_paneled_data(
          gene_list,
          filename,
          metafeatures_path,
          num_tss=num_tss,
          num_tes=num_tes
        )
      )
    )
  panel_data %>%
    rowwise %>%
    reframe(
      chic_mark,
      group,
      gene_list = list(gene_list),
      x = c(
        seq(0, tss_size, length.out = num_tss * 2),
        seq(tss_size, tss_size + inter_size, length.out = num_inter),
        seq(tss_size + inter_size, tss_size + inter_size + tes_size, length.out = num_tes * 2)
      ),
      x_label = track_data$labels,
      y = track_data %>% pull(2)
    ) %>%
    mutate(chic_mark = chic_mark %>% factor(chic.mark.data$mark, ordered=TRUE))
}

chic_quartile_gene_list_paneled_profiles <- function(
  quartile.factor,
  chic_path,
  metafeatures_path,
  chic_driver,
  legend_title,
  quartile_colors,
  heatmap_before_after = 500,
  # 1 kb - before and after
  tss_size = 1,
  # 2 kb long - by convention
  inter_size = 2,
  tes_size = 1
) {
  gene_lists <- names(quartile.factor) %>% split(quartile.factor)
  data <- pull_chic_average_gene_list_paneled_profiles_data(
    gene_lists,
    chic_path,
    metafeatures_path,
    chic_driver,
    tss_size,
    inter_size,
    tes_size,
    num_tss=heatmap_before_after,
    num_tes=heatmap_before_after
  ) %>%
    subset(select = -gene_list)

  label_data <- data %>%
    subset(
      as.numeric(chic_mark) == 1 & as.numeric(group) == 1,
      select = c(x, x_label)
    )
  label_data_show <- label_data %>%
    subset(
      x_label %in% c(
        paste0("-", heatmap_before_after),
        paste0("+", heatmap_before_after - 1),
        "50%",
        "TSS",
        "TES"
      )
    )
  minor_breaks_show <- label_data %>% subset(x_label %in% c("25%", "75%"))

  custom_plot <- data %>% ggplot(
    aes(x, y, color=group, linewidth=group, group=group)
  ) + geom_tile(
    # Colored tile (TSS and TES)
    data = tribble(
      ~x, ~y, ~width,
      # Full-size background color
      1,
      1,
      Inf,
      tss_size / 2,
      1,
      tss_size,
      tss_size + inter_size + (tes_size / 2),
      1,
      tes_size
    ) %>%
      cross_join(
        tibble(
          height = Inf
        )
      ),
    aes(width=width, height=height, color=NULL, linewidth=NULL, group=NULL),
    color = "transparent",
    fill = rep(
      # Use cream color from magma scale
      c("#eeeeee", viridis(2, option = "magma", end = 0.99)[2])[
        c(1, 2, 2)
      ],
      length(chic.mark.data$mark)
    )
  ) + facet_wrap(
    vars(chic_mark)
  ) + scale_color_manual(
    values = quartile_colors,
    guide = guide_legend(title = legend_title, reverse = TRUE, override.aes = list(fill = "transparent"))
  ) + scale_linewidth_manual(
    values = c(0.25, 0.5, 0.7, 1),
    guide = guide_legend(title = legend_title, reverse = TRUE)
  ) + coord_cartesian(
    c(0, tss_size + inter_size + tes_size),
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
  )
  custom_plot + geom_line(
    # Dark red line at the origin (enrichment of 1)
    data = cross_join(
      tribble(~x, ~y, -Inf, 1.01, Inf, 1.01),
      tibble(
        group = NA,
        chic_mark = factor(chic.mark.data$mark, chic.mark.data$mark, ordered=TRUE)
      )
    ),
    color = "darkred",
    linewidth = 0.33
  ) + geom_line() + labs(
    x = "base pairs", y = "log2(mean(mark/input))"
  ) + scale_x_continuous(
    breaks = label_data_show$x,
    labels = label_data_show$x_label,
    minor_breaks = minor_breaks_show$x
  ) + theme(
    panel.margin = unit(25, "pt"),
    aspect.ratio = 1
  )
}

chic_custom_gene_list_paneled_profile <- function(
  gene_list,
  bg_gene_list,
  chic_path,
  metafeatures_path,
  chic_driver,
  heatmap_before_after = 500,
  # 1 kb - before and after
  tss_size = 1,
  # 2 kb long - by convention
  inter_size = 2,
  tes_size = 1,
  track_color = "goldenrod"
) {
  data <- sapply(
    list(on = gene_list, off = bg_gene_list),
    \(gene_list) pull_chic_average_gene_list_paneled_profiles_data(
      list(gene_list),
      chic_path,
      metafeatures_path,
      chic_driver,
      tss_size,
      inter_size,
      tes_size,
      num_tss=heatmap_before_after,
      num_tes=heatmap_before_after
    ) %>%
      subset(select = -gene_list),
    simplify=FALSE
  ) %>%
    bind_rows(.id = "group")
  data$group <- data$group %>% factor(c("on", "off"))

  label_data <- data %>%
    subset(
      as.numeric(chic_mark) == 1 & as.numeric(group) == 1,
      select = c(x, x_label)
    )
  label_data_show <- label_data %>%
    subset(
      x_label %in% c(
        paste0("-", heatmap_before_after),
        paste0("+", heatmap_before_after - 1),
        "50%",
        "TSS",
        "TES"
      )
    )
  minor_breaks_show <- label_data %>% subset(x_label %in% c("25%", "75%"))

  custom_plot <- data %>% ggplot(
    aes(x, y)
  ) + geom_tile(
    # Colored tile (TSS and TES)
    data = tribble(
      ~x, ~y, ~width,
      # Full-size background color
      1,
      1,
      Inf,
      tss_size / 2,
      1,
      tss_size,
      tss_size + inter_size + (tes_size / 2),
      1,
      tes_size
    ) %>%
      cross_join(
        tibble(
          height = Inf
        )
      ),
    aes(width=width, height=height),
    color = "transparent",
    fill = rep(
      # Use cream color from magma scale
      c("#eeeeee", viridis(2, option = "magma", end = 0.99)[2])[
        c(1, 2, 2)
      ],
      length(chic.mark.data$mark)
    )
  ) + facet_wrap(
    vars(chic_mark)
  ) + coord_cartesian(
    c(0, tss_size + inter_size + tes_size),
    chic_average_profile_limits,
    expand=FALSE
  ) + theme(
    panel.background = element_rect(fill = NA),
    panel.ontop = TRUE
  )
  custom_plot + geom_line(
    # Dark red line at the origin (enrichment of 1)
    data = cross_join(
      tribble(~x, ~y, -Inf, 1, Inf, 1),
      tibble(
        group = NA,
        chic_mark = factor(chic.mark.data$mark, chic.mark.data$mark, ordered=TRUE)
      )
    ),
    color = "darkred",
    linewidth = 0.5
  ) + geom_line(
    # Implement the actual ChIC profile track.
    aes(group = group, color = group, linewidth = group)
  ) + labs(
    x = "base pairs", y = "log2(mean(mark/input))"
  ) + scale_x_continuous(
    breaks = label_data_show$x,
    labels = label_data_show$x_label,
    minor_breaks = minor_breaks_show$x
  ) + scale_y_continuous(
    trans = "log",
    labels = \(v) round(log(v) / log(2), 1),
    limits = chic_average_profile_limits,
    breaks = chic_average_breaks,
    minor_breaks = chic_average_minor_breaks,
    expand = c(0, 0)
  ) + scale_color_manual(
    values = c(track_color, muted(track_color, l=70))
  ) + scale_linewidth_manual(
    values = c(1, 0.25)
  ) + theme(
    panel.margin = unit(25, "pt"),
    aspect.ratio = 1
  )
}

# Draws any lines grouped by "genes", along K4/K27/K9 facet wrap and yellow
# TSS/TES tiles.
chic_multi_genes_paneled_profile <- function(
  custom_data,
  color,
  linewidth,
  heatmap_before_after = 500,
  # 1 kb - before and after
  tss_size = 1,
  # 2 kb long - by convention
  inter_size = 2,
  tes_size = 1
) {
  label_data <- custom_data %>%
    subset(
      as.numeric(chic_mark) == 1 & as.numeric(genes) == 1,
      select = c(x, x_label)
    )
  label_data_show <- label_data %>%
    subset(
      x_label %in% c(
        paste0("-", heatmap_before_after),
        paste0("+", heatmap_before_after - 1),
        "50%",
        "TSS",
        "TES"
      )
    )
  minor_breaks_show <- label_data %>% subset(x_label %in% c("25%", "75%"))

  custom_plot <- custom_data %>% ggplot(
    aes(x, y)
  ) + geom_tile(
    # Colored tile (TSS and TES)
    data = tribble(
      ~x, ~y, ~width,
      # Full-size background color
      1,
      1,
      Inf,
      tss_size / 2,
      1,
      tss_size,
      tss_size + inter_size + (tes_size / 2),
      1,
      tes_size
    ) %>%
      cross_join(
        tibble(
          height = Inf
        )
      ),
    aes(width=width, height=height),
    color = "transparent",
    fill = rep(
      # Use cream color from magma scale
      c("#eeeeee", viridis(2, option = "magma", end = 0.99)[2])[
        c(1, 2, 2)
      ],
      length(chic.mark.data$mark)
    )
  ) + facet_wrap(
    vars(chic_mark)
  ) + coord_cartesian(
    c(0, tss_size + inter_size + tes_size),
    chic_average_profile_limits,
    expand=FALSE
  ) + theme(
    panel.background = element_rect(fill = NA),
    panel.ontop = TRUE
  )
  custom_plot + geom_line(
    # Dark red line at the origin (enrichment of 1)
    data = cross_join(
      tribble(~x, ~y, -Inf, 1, Inf, 1),
      tibble(
        genes = NA,
        chic_mark = factor(chic.mark.data$mark, chic.mark.data$mark, ordered=TRUE)
      )
    ),
    color = "darkred",
    linewidth = 0.5
  ) + geom_line(
    # Implement the actual ChIC profile track.
    aes(group = genes, color = genes, linewidth = genes)
  ) + labs(
    x = "base pairs", y = "log2(mean(mark/input))"
  ) + scale_x_continuous(
    breaks = label_data_show$x,
    labels = label_data_show$x_label,
    minor_breaks = minor_breaks_show$x
  ) + scale_y_continuous(
    trans = "log",
    labels = \(v) round(log(v) / log(2), 1),
    limits = chic_average_profile_limits,
    breaks = chic_average_breaks,
    minor_breaks = chic_average_minor_breaks,
    expand = c(0, 0)
  ) + scale_color_manual(
    values = color
  ) + scale_linewidth_manual(
    values = linewidth
  ) + theme(
    panel.margin = unit(25, "pt"),
    aspect.ratio = 1
  )
}
