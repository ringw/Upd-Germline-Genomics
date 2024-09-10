plot_track <- function(
    track,
    positive_color = "#e0eff9",
    negative_color = "#f7f9e7",
    name = "Timing",
    limits = c(-1, 1),
    breaks = c(-1, 0, 1),
    subsample_draw = TRUE) {
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
  if (subsample_draw)
    df <- df %>%
      dplyr::slice(
        c(
          round(seq(1, -1 + match(chr_levels$`3`, chr_lookup), length.out = 1000)),
          round(seq(match(chr_levels$`3`, chr_lookup), -1 + match(chr_levels$`4`, chr_lookup), length.out = 1000)),
          round(seq(match(chr_levels$`4`, chr_lookup), -1 + match(chr_levels$X, chr_lookup))),
          round(seq(match(chr_levels$X, chr_lookup), -1 + match(chr_levels$Y, chr_lookup), length.out = 500)),
          round(seq(match(chr_levels$Y, chr_lookup), nrow(df), length.out = 1000))
        )
      )
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
  chrX <- (
    plot_chr("X") +
      labs(title = "X")
  ) %>%
    ggplot_build_panel_absolute(
      unit(0.35, "in"), unit(4, "in"), unit(5.5, "points") + unit(1, "in")
    )
  chr2 <- (
    plot_chr("2") +
      annotate("line", rep(chr.lengths["2L"], 2), c(-Inf, Inf)) +
      labs(title = "2")
  ) %>%
    ggplot_build_panel_absolute(unit(0.35, 'in'), unit(5, 'in'))
  chr3 <- (
    plot_chr("3") +
      annotate("line", rep(chr.lengths["3L"], 2), c(-Inf, Inf)) +
      labs(title = "3")
  ) %>%
    ggplot_build_panel_absolute(unit(0.35, 'in'), unit(5, 'in'))
  chr4 <- (
    plot_chr("4") +
      labs(title = "4")
  ) %>%
    ggplot_build_panel_absolute(
      unit(0.35, 'in'), unit(2, 'in'), unit(5.5, 'points') + unit(3, 'in')
    )
  chrY <- (
    plot_chr("Y") +
      labs(title = "Y")
  ) %>%
    ggplot_build_panel_absolute(
      unit(0.35, "in"), unit(2, "in"), unit(5.5, "points") + unit(3, "in")
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
