# Line-plotting of a feature on chr X/2/3/4/Y. ----
plot_track <- function(
    track,
    positive_color = "#e0eff9",
    negative_color = "#f7f9e7",
    name = "Timing",
    limits = c(-1, 1),
    breaks = c(-1, 0, 1)) {
  df <- tibble(
    chr = as.factor(seqnames(track)),
    pos = mid(track),
    value = track$score
  ) %>%
    subset(chr %in% names(chr.lengths))
  chr_lookup <- as.numeric(df$chr)
  chr_levels <- list(`2` = 1, `3` = 2, `4` = 3, `X` = 4, `Y` = 5)
  df <- df %>%
    dplyr::slice(
      df %>%
        group_by(chr) %>%
        reframe(
          keep = seq(length(pos)) %in% round(seq(1, length(pos), length.out=1000)) |
            chr == "4"
        ) %>%
        pull(keep) %>%
        which()
    )
  plot_chr <- function(chr.) ggplot(subset(df, chr == chr.), aes(pos, value)) +
    annotate(
      "rect",
      xmin=-Inf, xmax=Inf,
      ymin=-Inf, ymax=mean(limits),
      fill=negative_color
    ) +
    annotate(
      "rect",
      xmin=-Inf, xmax=Inf,
      ymin=Inf, ymax=mean(limits),
      fill=positive_color
    ) +
    geom_line(linewidth = 0.25) +
    scale_x_continuous(
      name = NULL,
      breaks = c(1, 1000000 * seq(2, 100, by=2)),
      minor_breaks = 1000000 * seq(1, 101, by=2),
      labels = NULL
    ) +
    scale_y_continuous(
      name = name,
      limits = limits,
      breaks = breaks,
      oob = scales::squish
    ) +
    coord_cartesian(
      NULL, limits, expand=F
    ) +
    theme_bw() +
    theme(
      panel.background = element_rect(fill = NA),
      panel.border = element_blank(),
      panel.grid = element_blank()
    )
  ggplot_build_panel_absolute <- function(gg, height, width, margin_right = unit(5.5, "points")) {
    gr <- as_grob(gg)
    gr$heights[7] <- height
    gr$widths[5] <- width
    gr$widths[9] <- margin_right
    gr <- gr
  }
  # Contents width will be 5.25 in. Two-column (chromosome arms) layout will be
  # handled by the two arms having precisely 2 * default margins (5.5 pt) in
  # between the columns, and no other middle content (axis title, axis text).
  chrX <- (
    plot_chr("X") +
      labs(title = "Chromosome X")
  ) %>%
    ggplot_build_panel_absolute(
      unit(0.35, "in"), unit(4.5, "in"), unit(0.55, "in") + unit(11, "points")
    )
  chr2 <- cbind(
    (plot_chr("2L") + labs(title = "Chromosome 2")) %>%
      ggplot_build_panel_absolute(unit(0.35, "in"), unit(2.5, "in")),
    (
      plot_chr("2R") +
        theme(
          axis.title = element_blank(),
          axis.text = element_blank()
        )
    ) %>%
      ggplot_build_panel_absolute(
        unit(0.35, "in"), unit(2.5, "in"), unit(0.05, "in")
      )
  )
  chr3 <- cbind(
    (plot_chr("3L") + labs(title = "Chromosome 3")) %>%
      ggplot_build_panel_absolute(unit(0.35, "in"), unit(2.5, "in")),
    (
      plot_chr("3R") +
        theme(
          axis.title = element_blank(),
          axis.text = element_blank()
        )
    ) %>%
      ggplot_build_panel_absolute(
        unit(0.35, "in"), unit(2.5, "in"), unit(0.05, "in")
      )
  )
  chr4 <- (
    plot_chr("4") +
      labs(title = "Chromosome 4")
  ) %>%
    ggplot_build_panel_absolute(
      unit(0.35, "in"), unit(2, "in"), unit(3.05, "in") + unit(11, "points")
    )
  chrY <- (
    plot_chr("Y") +
      labs(title = "Chromosome Y")
  ) %>%
    ggplot_build_panel_absolute(
      unit(0.35, "in"), unit(2, "in"), unit(3.05, "in") + unit(11, "points")
    )
  mylayout <- gtable(
    widths = unit(1, 'null'), heights = unit(rep(0.75, 5), 'in')
  ) %>%
    gtable_add_grob(
      list(
        chrX, chr2, chr3, chr4, chrY
      ),
      1:5,
      1
    )
  plot_grid(mylayout)
}

# Plot two tracks, col names score_1 and score_2, overlaid.
plot_track_2score <- function(
    track,
    positive_color = "#e0eff9",
    negative_color = "#f7f9e7",
    background_color = NA,
    line_color = NA,
    line_width = 0.5,
    score_1 = chic_line_track_colors$germline,
    score_2 = chic_line_track_colors$somatic,
    roi = data.frame(
      chr = character(0),
      xmin = numeric(0),
      xmax = numeric(0),
      ymin = numeric(0),
      ymax = numeric(0),
      fill = character(0)
    ),
    impute_score = T,
    name = "Timing",
    limits = c(-1, 1),
    breaks = c(-1, 0, 1)) {
  df <- apply(
    elementMetadata(track)[, c("score_1", "score_2")],
    2,
    \(v) tibble(
      chr = as.factor(seqnames(track)),
      pos = mid(track),
      value = v
    ) %>%
      subset(chr %in% names(chr.lengths)),
    simplify = FALSE
  ) %>%
    setNames(c("1", "2")) %>%
    bind_rows(.id = "group")
  chr_lookup <- as.numeric(df$chr)
  chr_levels <- list(`2` = 1, `3` = 2, `4` = 3, `X` = 4, `Y` = 5)
  # Impute geom_line: Remove all NA values before further subsampling the track
  # if needed. Gaps in the track will only be evident at the start and end of
  # the sequence. Otherwise, it is a geom_line where we draw a line segment
  # connecting the values which are non-NA.
  if (impute_score)
    df <- df %>% subset(!is.na(value))
  df <- df %>%
    dplyr::slice(
      df %>%
        group_by(group, chr) %>%
        reframe(
          keep = seq(length(pos)) %in% round(seq(1, length(pos), length.out=1000)) |
            chr == "4"
        ) %>%
        pull(keep) %>%
        which()
    )
  plot_chr <- function(chr.) ggplot(subset(df, chr == chr.), aes(pos, value, color=group, group=group)) +
    annotate(
      "rect",
      xmin=-Inf, xmax=Inf,
      ymin=-Inf, ymax=mean(limits),
      fill=negative_color
    ) +
    annotate(
      "rect",
      xmin=-Inf, xmax=Inf,
      ymin=Inf, ymax=mean(limits),
      fill=positive_color
    ) +
    geom_rect(
      aes(x = NULL, y = NULL, color = NULL, group = NULL, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
      data = subset(roi, chr == chr.)
    ) +
    geom_line(linewidth = line_width * 25.4 / 72) +
    scale_x_continuous(
      name = NULL,
      breaks = c(1, 1000000 * seq(2, 100, by=2)),
      minor_breaks = 1000000 * seq(1, 101, by=2),
      labels = NULL
    ) +
    scale_y_continuous(
      name = name,
      breaks = breaks
    ) +
    scale_color_manual(values = c(`1`=score_1, `2`=score_2)) +
    scale_fill_identity() +
    coord_cartesian(
      NULL, limits, expand=F
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = background_color),
      axis.ticks = if (is.na(line_color))
          waiver()
        else
          element_line(color = line_color),
      panel.border = if (is.na(line_color))
          element_blank()
        else
          element_rect(color = line_color),
      panel.grid = element_blank(),
    )
  ggplot_build_panel_absolute <- function(gg, height, width, margin_right = unit(5.5, "points")) {
    gr <- as_grob(gg)
    gr$heights[7] <- height
    gr$widths[5] <- width
    gr$widths[9] <- margin_right
    gr <- gr
  }
  # Contents width will be 5.25 in. Two-column (chromosome arms) layout will be
  # handled by the two arms having precisely 2 * default margins (5.5 pt) in
  # between the columns, and no other middle content (axis title, axis text).
  chrX <- (
    plot_chr("X") +
      labs(title = "Chromosome X")
  ) %>%
    ggplot_build_panel_absolute(
      unit(0.35, "in"), unit(4.5, "in"), unit(0.55, "in") + unit(11, "points")
    )
  chr2 <- cbind(
    (plot_chr("2L") + labs(title = "Chromosome 2")) %>%
      ggplot_build_panel_absolute(unit(0.35, "in"), unit(2.5, "in")),
    (
      plot_chr("2R") +
        theme(
          axis.title = element_blank(),
          axis.text = element_blank()
        )
    ) %>%
      ggplot_build_panel_absolute(
        unit(0.35, "in"), unit(2.5, "in"), unit(0.05, "in")
      )
  )
  chr3 <- cbind(
    (plot_chr("3L") + labs(title = "Chromosome 3")) %>%
      ggplot_build_panel_absolute(unit(0.35, "in"), unit(2.5, "in")),
    (
      plot_chr("3R") +
        theme(
          axis.title = element_blank(),
          axis.text = element_blank()
        )
    ) %>%
      ggplot_build_panel_absolute(
        unit(0.35, "in"), unit(2.5, "in"), unit(0.05, "in")
      )
  )
  chr4 <- (
    plot_chr("4") +
      labs(title = "Chromosome 4")
  ) %>%
    ggplot_build_panel_absolute(
      unit(0.35, "in"), unit(2, "in"), unit(3.05, "in") + unit(11, "points")
    )
  chrY <- (
    plot_chr("Y") +
      labs(title = "Chromosome Y")
  ) %>%
    ggplot_build_panel_absolute(
      unit(0.35, "in"), unit(2, "in"), unit(3.05, "in") + unit(11, "points")
    )
  mylayout <- gtable(
    widths = unit(1, 'null'), heights = unit(rep(0.75, 5), 'in')
  ) %>%
    gtable_add_grob(
      list(
        chrX, chr2, chr3, chr4, chrY
      ),
      1:5,
      1
    )
  # plot_grid(mylayout)
}

plot_genomic_ranges_score <- function(
    track,
    positive_color = "#e0eff9",
    negative_color = "#f7f9e7",
    name = "Timing",
    limits = c(-1, 1),
    breaks = c(-1, 0, 1)) {
  df <- tibble(
    chr = as.factor(seqnames(track)),
    pos = mid(track),
    value = track$score
  )
  df$group <- cumsum(
    c(
      1,
      end(track)[-length(track)] + 1 != start(track)[-1] |
        diff(as.numeric(seqnames(track))) != 0
    )
  )
  df <- df %>% subset(chr %in% names(chr.lengths))
  df <- df %>%
    group_by(group) %>%
    dplyr::slice(
      unique(round(seq(1, length(pos), length.out=pmax(3, round(diff(range(pos)) / 10000)))))
    )
  chr_lookup <- as.numeric(df$chr)
  chr_levels <- list(`2` = 1, `3` = 2, `4` = 3, `X` = 4, `Y` = 5)
  plot_chr <- function(chr.) ggplot(subset(df, chr == chr.), aes(pos, value, group=group)) +
    annotate(
      "rect",
      xmin=-Inf, xmax=Inf,
      ymin=-Inf, ymax=0,
      fill=negative_color
    ) +
    annotate(
      "rect",
      xmin=-Inf, xmax=Inf,
      ymin=Inf, ymax=0,
      fill=positive_color
    ) +
    geom_line(linewidth = 0.25) +
    scale_x_continuous(
      name = NULL,
      breaks = c(1, 1000000 * seq(2, 100, by=2)),
      minor_breaks = 1000000 * seq(1, 101, by=2),
      labels = NULL
    ) +
    scale_y_continuous(
      name = name,
      limits = limits,
      breaks = breaks,
      oob = scales::squish
    ) +
    coord_cartesian(
      NULL, limits, expand=F
    ) +
    theme_bw() +
    theme(
      panel.background = element_rect(fill = NA),
      panel.border = element_blank(),
      panel.grid = element_blank()
    )
  ggplot_build_panel_absolute <- function(gg, height, width, margin_right = unit(5.5, "points")) {
    gr <- as_grob(gg)
    gr$heights[7] <- height
    gr$widths[5] <- width
    gr$widths[9] <- margin_right
    gr <- gr
  }
  # Contents width will be 5.25 in. Two-column (chromosome arms) layout will be
  # handled by the two arms having precisely 2 * default margins (5.5 pt) in
  # between the columns, and no other middle content (axis title, axis text).
  chrX <- (
    plot_chr("X") +
      labs(title = "Chromosome X")
  ) %>%
    ggplot_build_panel_absolute(
      unit(0.35, "in"), unit(4.5, "in"), unit(0.55, "in") + unit(11, "points")
    )
  chr2 <- cbind(
    (plot_chr("2L") + labs(title = "Chromosome 2")) %>%
      ggplot_build_panel_absolute(unit(0.35, "in"), unit(2.5, "in")),
    (
      plot_chr("2R") +
        theme(
          axis.title = element_blank(),
          axis.text = element_blank()
        )
    ) %>%
      ggplot_build_panel_absolute(
        unit(0.35, "in"), unit(2.5, "in"), unit(0.05, "in")
      )
  )
  chr3 <- cbind(
    (plot_chr("3L") + labs(title = "Chromosome 3")) %>%
      ggplot_build_panel_absolute(unit(0.35, "in"), unit(2.5, "in")),
    (
      plot_chr("3R") +
        theme(
          axis.title = element_blank(),
          axis.text = element_blank()
        )
    ) %>%
      ggplot_build_panel_absolute(
        unit(0.35, "in"), unit(2.5, "in"), unit(0.05, "in")
      )
  )
  chr4 <- (
    plot_chr("4") +
      labs(title = "Chromosome 4")
  ) %>%
    ggplot_build_panel_absolute(
      unit(0.35, "in"), unit(2, "in"), unit(3.05, "in") + unit(11, "points")
    )
  chrY <- (
    plot_chr("Y") +
      labs(title = "Chromosome Y")
  ) %>%
    ggplot_build_panel_absolute(
      unit(0.35, "in"), unit(2, "in"), unit(3.05, "in") + unit(11, "points")
    )
  mylayout <- gtable(
    widths = unit(1, 'null'), heights = unit(rep(0.75, 5), 'in')
  ) %>%
    gtable_add_grob(
      list(
        chrX, chr2, chr3, chr4, chrY
      ),
      1:5,
      1
    )
  plot_grid(mylayout)
}

plot_points <- function(
    track,
    positive_color = "#caff84",
    negative_color = "#ee9999",
    name = "Timing",
    limits = c(-1, 1),
    breaks = c(-1, 0, 1)) {
  seq_offsets <- seqnames(track) %>%
    as.factor() %>%
    fct_relabel(
      \(names) names %>%
        sapply(\(n) if (n == "2R") chr.lengths["2L"] else if (n == "3R") chr.lengths["3L"] else 0) %>%
        as.character()
    )
  df <- tibble(
    chr = fct_recode(as.factor(seqnames(track)), `2` = "2L", `2` = "2R", `3` = "3L", `3` = "3R"),
    pos = mid(track) + as.numeric(as.character(seq_offsets)),
    value = track$score
  ) %>%
    subset(chr %in% c("2", "3", "4", "X", "Y"))
  chr_lookup <- as.numeric(df$chr)
  chr_levels <- list(`2` = 1, `3` = 2, `4` = 3, `X` = 4, `Y` = 5)
  plot_chr <- function(chr.) ggplot(subset(df, chr == chr.), aes(pos, value)) +
    annotate(
      "rect",
      xmin=-Inf, xmax=Inf,
      ymin=-Inf, ymax=0,
      fill=negative_color
    ) +
    annotate(
      "rect",
      xmin=-Inf, xmax=Inf,
      ymin=Inf, ymax=0,
      fill=positive_color
    ) +
    geom_point(size=0.75, stroke=NA) +
    scale_x_continuous(
      name = NULL,
      breaks = c(1, 1000000 * seq(2, 100, by=2)),
      minor_breaks = 1000000 * seq(1, 101, by=2),
      labels = NULL
    ) +
    scale_y_continuous(name = name, breaks = breaks) +
    coord_cartesian(
      NULL, limits, expand=F
    ) +
    theme_bw() +
    theme(
      panel.background = element_rect(fill = NA),
      panel.ontop = TRUE
    )
  chrX <- as_grob(
    plot_chr("X") +
      labs(title = "X")
  )
  chrX$heights[7] <- unit(0.35, 'in')
  chrX$widths[5] <- unit(4, 'in')
  chrX$widths[9] <- unit(5.5, 'points') + unit(1, 'in')
  chr2 <- as_grob(
    plot_chr("2") +
      annotate("line", rep(chr.lengths["2L"], 2), c(-Inf, Inf)) +
      labs(title = "2")
  )
  chr2$heights[7] <- unit(0.35, 'in')
  chr2$widths[5] <- unit(5, 'in')
  chr3 <- as_grob(
    plot_chr("3") +
      annotate("line", rep(chr.lengths["3L"], 2), c(-Inf, Inf)) +
      labs(title = "3")
  )
  chr3$heights[7] <- unit(0.35, 'in')
  chr3$widths[5] <- unit(5, 'in')
  chr4 <- as_grob(
    plot_chr("4") +
      labs(title = "4")
  )
  chr4$heights[7] <- unit(0.35, 'in')
  chr4$widths[5] <- unit(2, 'in')
  chr4$widths[9] <- unit(5.5, 'points') + unit(3, 'in')
  chrY <- as_grob(
    plot_chr("Y") +
      labs(title = "Y")
  )
  chrY$heights[7] <- unit(0.35, 'in')
  chrY$widths[5] <- unit(2, 'in')
  chrY$widths[9] <- unit(5.5, 'points') + unit(3, 'in')
  mylayout <- gtable(
    widths = unit(1, 'null'), heights = unit(rep(0.75, 5), 'in')
  ) %>%
    gtable_add_grob(
      list(
        chrX, chr2, chr3, chr4, chrY
      ),
      1:5,
      1
    )
  plot_grid(mylayout)
}
