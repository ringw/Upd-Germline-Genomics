approx_track <- function(score_track, fine_track) {
  stopifnot(all.equal(as.character(unique(seqnames(fine_track))), intersect(seqlevels(fine_track), unique(seqnames(fine_track)))))
  score_track <- score_track[score_track@seqnames %in% levels(fine_track@seqnames)]
  seqlevels(score_track) <- levels(seqnames(fine_track))

  # The two tracks can now be split according to the same levels (names).
  score_track <- split(score_track, seqnames(score_track))
  # First, initialize the score track from the first entry in each coarse
  # (score_track) seq's first value. Note that some of the score_track seqs have
  # more than one range in them, while the short seqs have a single window
  # covering the entire seq and benefit from this optimization.
  fine_track$score <- rep(
    sapply(score_track, \(gr) gr$score[1]) %>%
      replace(. == 0, NA),
    unlist(table(seqnames(fine_track)))
  )
  fine_track <- split(fine_track, seqnames(fine_track))
  # Now loop over the seqs that are actually affected by the finer windows that
  # are found in the fine_track.
  fine_seqs <- names(score_track)[sapply(score_track, length) > 1]
  for (n in fine_seqs) {
    x <- score_track[[n]]@ranges@start + round(score_track[[n]]@ranges@width / 2)
    y <- score_track[[n]]$score
    # rtracklayer forces us to replace NA with 0.
    x <- x[y != 0]
    y <- y[y != 0]
    if (length(y) == 0)
      fine_track[[n]]$score <- NA
    else if (length(y) == 1)
      fine_track[[n]]$score <- y
    else
      fine_track[[n]]$score <- approx(
        x,
        y,
        xout = fine_track[[n]]@ranges@start + round(fine_track[[n]]@ranges@width / 2),
        rule = 2
      )$y
  }
  unlist(fine_track %>% setNames(NULL))
}

plot_repli_chic_bin2d <- function(repli, chic, chic_step_size = 20) {
  repli.timing <- approx_track(repli, GRanges(seqnames(chic), ranges(chic), seqinfo=seqinfo(chic)))$score
  plt <- tibble(`Replication Timing`=repli.timing, Enrichment=chic$score) %>%
    filter(between(Enrichment, -2, 2)) %>%
    ggplot(aes(`Replication Timing`, Enrichment)) +
    geom_bin2d(bins=c(120, 80)) +
    scale_x_reverse(
      breaks = c(-0.75, -0.375, 0, 0.375, 0.75),
      labels = purrr::partial(format, drop0 = TRUE)
    ) +
    coord_cartesian(
      c(1, -1),
      c(-2, 2),
      expand=F
    ) +
    theme_bw() +
    theme(aspect.ratio = 2/3, panel.background = element_rect(fill=NA), panel.ontop = TRUE)
  printable_plt <- ggplot_build(plt)
  max_data <- min(tail(sort(printable_plt$data[[1]]$count), 2))
  max_data <- max(printable_plt$data[[1]]$count)
  scl <- scale_fill_gradientn(
    name = "kb covered",
    labels = \(n) n * chic_step_size / 1000,
    colors = c(
      "white",
      rev(magma(101))
    ),
    values = c(
      0,
      seq(
        median(printable_plt$data[[1]]$count) / max_data,
        1,
        length.out = 101
      )
    ),
    limits = c(0, max_data),
    na.value = "#CCCCCC"
  )
  plt + scl
}

violin_plot_repli_chic <- function(df) {
  df <- df %>%
    mutate(
      timing = timing %>% fct_relabel(\(names) sapply(names, \(n) str_glue(n, " (", scales::percent(mean(df$timing == n)), ")")))
    )
  df %>%
    ggplot(aes(timing, enrichment, fill=timing)) +
    annotate("line", c(-Inf,Inf), 0) +
    geom_violin() +
    geom_boxplot(outlier.shape = NA, fill = "transparent") +
    scale_x_discrete(labels = names(repli_level_colors)) +
    scale_fill_manual(values = unlist(repli_level_colors, use.names=F)) +
    coord_cartesian(NULL, c(-4, 4)) +
    theme_cowplot() +
    theme(aspect.ratio = 1.25)
}

bin_track_by_rank <- function(track, nbins = 1000) {
  # newnames <- with(
  #   attributes(track[!is.na(track$rank)]),
  #   str_glue("{seqnames}_{ranges@start}")
  # )
  newperm <- order(track$rank[!is.na(track$rank)])
  newlookup <- order(track$rank) %>% head(sum(!is.na(track$rank)))
  newtrack <- track@ranges@width[newlookup] %>%
    IRanges(
      start = cumsum(c(1, .[-length(.)])),
      width = . # ,
      # names = newnames[newperm]
    )
  bin_loc <- seq(1, sum(newtrack@width), length.out=1+nbins)[-(1+nbins)]
  bin_fct <- rep(
    seq(nbins),
    diff(c(0, findInterval(bin_loc[-1], newtrack@start), length(newtrack)))
  )
  bin <- rep(as.integer(NA), length(newtrack))
  bin[newlookup] <- bin_fct
  bin_sizes_iranges <- newtrack@width %>% split(bin_fct) %>% sapply(sum)
  size_of_bin_bp <- rep(as.integer(NA), length(newtrack))
  size_of_bin_bp[newlookup] <- bin_sizes_iranges[bin_fct]
  GRanges(track, bin = bin, size_of_bin_bp = size_of_bin_bp)
}

cut_track <- function(track, breaks) {
  result <- cut(track$score, breaks)
  splt <- width(track) %>%
    split(result) %>%
    sapply(sum)
  size_of_bin_bp <- setNames(splt[result], NULL)
  GRanges(track, score = track$score, bin = as.numeric(result), size_of_bin_bp = size_of_bin_bp)
}

granges_bin_to_projection <- function(granges) {
  proj_nrow <- max(granges$bin, na.rm=T)
  granges %>%
    split(seqnames(.)) %>%
    sapply(
      \(gr) sparseMatrix(
        i = gr$bin[!is.na(gr$bin)],
        j = which(!is.na(gr$bin)),
        x = (gr@ranges@width / gr$size_of_bin_bp)[!is.na(gr$bin)],
        dims = c(proj_nrow, length(gr))
      ),
      simplify=F
    )
}

apply_granges_projection <- function(projs, granges) {
  rowSums(
    mapply(
      \(mat, v) (mat %*% v)[,1],
      projs,
      split(granges$score, factor(seqnames(granges), names(projs)))
    )
  )
}

repli_chromosome_arms_factor_scaffolds <- function(chromosome.arms) {
  fct.to.replace <- seqnames(chromosome.arms) %in%
    c("4", "X", "Y")
  fct.to.scaffold <- !(
    seqnames(chromosome.arms) %in% names(chr.lengths)
  )
  fct <- chromosome.arms$group %>%
    factor(c(levels(.), "4", "X", "Y", "scaffolds")) %>%
    replace(
      which(fct.to.replace),
      as.character(seqnames(chromosome.arms)[fct.to.replace])
    ) %>%
    replace(which(fct.to.scaffold), "scaffolds")
}

# Plotting functions for Repli ChIC heatmap ----
weight_repli_mode_prop_table <- function(
  repli.experiment,
  chromosome.arms
) {
  Y <- repli.experiment %>%
    assay() %>%
    as.data.frame() %>%
    group_by(rowname = chromosome.arms) %>%
    summarize_all(sum) %>%
    column_to_rownames() %>%
    as.matrix()
  g <- glm_gp(Y, verbose = TRUE)
  Yhat <- g %>% predict(newdata = matrix(1), offset = 0, type = "response") %>%
    apply(1, identity)
  other_margin <- 1000 * 1000 * Yhat / sum(Yhat)
  1 / other_margin
}

my_mix_color_correct_overplotting <- function(mix_columns) {
  # sums <- rowSums(mix_columns)
  # mix_columns <- mix_columns *
  #   ifelse(sums > 1, 1/sums, 1)
  cumulative_color <- 0
  for (i in seq(ncol(mix_columns))) {
    mix_columns[, i] <- mix_columns[, i] %>%
      pmin(1 - cumulative_color)
    cumulative_color <- cumulative_color + mix_columns[, i]
  }
  mix_columns <- mix_columns
}

my_mix_color <- function(mix_columns, mix_colors) {
  cumulative_paint <- 1 - rowSums(mix_columns, na.rm=T)
  result <- hex2RGB("#ffffff")
  for (i in seq(ncol(mix_columns))) {
    result <- mixcolor(
      mix_columns[, i] / (mix_columns[, i] + cumulative_paint),
      result,
      mix_colors[i],
      "LAB"
    )
    cumulative_paint <- cumulative_paint + mix_columns[, i]
  }
  result <- result
}

plot_repli_track_raster <- function(
  data, log2_limits = c(-0.5, 0.52), repli.mode.chr.weights = 1
) {
  # All columns except for the "repli" column may need to be NA out.
  data[, -1] <- data[, -1] +
    ifelse(
      is.na(data[, "sample_size_bp"]) |
        data[, "sample_size_bp"] < 500,
      NA,
      0
    )
  if (any(grepl("TSS", colnames(data)))) {
    sc_track <- tibble(
      series = "TSS Density",
      x = data[, "repli"],
      fill = my_mix_color(
        data[, c("TSS_high", "TSS_medium", "TSS_low", "TSS_off")] %>%
          my_mix_color_correct_overplotting(),
        hex2RGB(rev(sc_quartile_annotations))
      ) %>%
        hex()
    )
  } else {
    sc_track <- NULL
  }
  chr_track <- tibble(
    series = "Mode Chr.",
    x = data[, "repli"],
    fill = data[, c("chr2L", "chr2LC", "chr2RC", "chr2R", "chr3L", "chr3LC", "chr3RC", "chr3R", "chr4", "chrX", "chrY", "chrscaffolds")] %>%
      `*`(rep(repli.mode.chr.weights, each = nrow(.))) %>%
      apply(1, \(v) if (any(is.finite(v))) which.max(v) else NA) %>%
      `[`(c(arm.colors, "#ffffff"), value = .)
  )
  sample_size_mb <- round(sum(data[-match(0, data[, "repli"]), "sample_size_bp"], na.rm=T) / 1000 / 1000, 1)
  data[, grep("^H", colnames(data))] <- (
    log(data[, grep("^H", colnames(data))]) / log(2)
  )
  data <- melt(data) %>%
    mutate(series = series %>% fct_recode(`Timing Est.`="repli")) %>%
    group_by(timing) %>%
    mutate(
      series = series %>%
        factor(c(levels(.), "Quartile")) %>%
        fct_relevel("Timing Est.", "Quartile"),
      x = value[which.max(series == "ranking")]
    )
  x_values <- unique(data$value[data$series == "ranking"])
  range_x <- range(data$x[grepl("^H", data$series) & !is.na(data$value)])
  rnd_x <- round(range_x, 2)
  gg <- (
    ggplot(data, aes(x, series))
    + geom_raster(aes(fill=value), subset(data, series == "Timing Est."))
    + scale_fill_gradientn(
      colors = c(
        hcl(30, 40, 10),
        repli_level_colors$L,
        repli_level_colors$ML,
        repli_level_colors$EM,
        repli_level_colors$E,
        hcl(170, 55, 95)
      ),
      values = c(0, 0.125, 0.375, 0.625, 0.875, 1),
      guide = guide_colorbar(title = "timing"),
      limits = c(-1, 1),
      breaks = c(rnd_x[1], -0.5, 0, 0.5, rnd_x[2]),
      labels = c(rnd_x[1], "-0.5", "0", "0.5", rnd_x[2]),
      na.value = "transparent"
    )
  )
  if (!is.null(sc_track)) {
    gg <- gg +
      new_scale_fill() +
      geom_raster(aes(fill=fill), sc_track) +
      scale_fill_identity()
  }
  gg <- gg +
    new_scale_fill() +
    geom_raster(aes(fill=fill), chr_track) +
    scale_fill_identity() +
    new_scale_fill()
  raster_height <- 1.01
  for (n in grep("^H", levels(data$series), val = TRUE)) {
    gg <- gg + geom_blank(data = data[head(match(n, data$series), 1),, drop=FALSE])
    gg <- gg +
      geom_raster(
        aes(fill=value2),
        mutate(
          cross_join(
            subset(data, series == n),
            tibble(heightadjust = c(-0.25, 0.25) * raster_height)
          ),
          value2 = value,
          series = c(
            H3K4=3,
            H3K27=2,
            H3K9=1
          )[as.character(series)] +
            heightadjust
        )
      )
  }
  gg <- gg +
    create_direction_invert_tss_tile_matrix_gradient(
      invert_tss_limits = if (is.null(log2_limits))
        range(data$value[grepl("^H", data$series)], na.rm=T)
      else
        log2_limits,
      breaks = pretty_breaks(3),
      limits = log2_limits, oob = squish, na.value = "transparent"
    ) +
    scale_x_reverse(
      breaks = c(range_x[1], -.5, 0, .5, range_x[2]),
      labels = c(
        round(range_x[1], 2),
        "-0.5",
        "0",
        "0.5",
        round(range_x[2], 2)
      ),
      name = str_glue("Sorted Timing Est ({sample_size_mb} Mb)")
    ) +
    scale_y_discrete(limits=rev) +
    coord_cartesian(rev(range_x), NULL, expand=F) +
    guides(fill = guide_colorbar(title = "log2(mark/input)")) +
    theme(
      aspect.ratio = 2/3,
      panel.border = element_blank(),
      panel.grid.major.y = element_blank()
    )
}

# Plotting with chromosome arms legend wrapped along the bottom edge ----
plot_tile_example_arm_colors <- function() {
  ggplot(
    tibble(
      name = names(arm.colors) %>% factor(., .),
      x = seq_along(name)
    ),
    aes(x, 0, fill=name)
  ) +
    geom_tile() +
    scale_fill_manual(
      NULL,
      labels = \(text) paste0(
        text,
        rep(
          c(
            strrep(" ", 3),
            ""
          ),
          c(length(text) - 1, 1)
        )
      ),
      values = arm.colors,
      guide = guide_legend(nrow = 1)
    ) +
    theme(
      legend.key.size = unit(0.12, "in"),
      legend.position = "bottom",
    )
}

grob_add_arm_colors_legend <- function(p, width, height) {
  arm.colors.legend.height <- unit(0.5, "in")
  gr <- gtable(
    width,
    unit.c(height, arm.colors.legend.height)
  ) %>%
    gtable_add_grob(
      list(ggplotGrob(p), get_legend(plot_tile_example_arm_colors())),
      t = 1:2,
      l = 1
    )
}
