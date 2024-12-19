source("_targets.R")
sapply(tar_option_get("packages"), \(n) require(n, ch=T))
library(grid)
library(grobblR)
library(gtable)
library(purrr)

tar_load(matches("repli.timing_") | repli.peaks_chr)

tracks <- list(
  repli.timing_Germline_chr,
  repli.timing_Somatic_chr,
  repli.timing_Kc167_chr,
  repli.timing_S2_chr
)
line_cmap <- c("#6D9965", "#B377B0", "#C27A00", "#4087F4")
bg_colors <- c("#255529", "#7200ab", "#767323")
repli_peaks_lookup <- with(
  repli.peaks_chr,
  list(
    "1" = list(
      "2" = Germline_Somatic,
      "3" = Germline_Kc167,
      "4" = Germline_S2
    ),
    "2" = list(
      "3" = Somatic_Kc167,
      "4" = Somatic_S2
    ),
    "3" = list(
      "4" = Kc167_S2
    )
  )
)
repli_titles_lookup <- with(
  repli.peaks_chr,
  list(
    "1" = list(
      "2" = "Germline - Somatic",
      "3" = "Germline - Kc167",
      "4" = "Germline - S2"
    ),
    "2" = list(
      "3" = "Somatic - Kc167",
      "4" = "Somatic - S2"
    ),
    "3" = list(
      "4" = "Kc167 - S2"
    )
  )
)

plot_track_custom <- function(i, j) {
  peaks <- repli_peaks_lookup[[as.character(i)]][[as.character(j)]]
  title <- repli_titles_lookup[[as.character(i)]][[as.character(j)]]
  width <- unit(6, "in")
  gr <- plot_track_2score(
    GRanges(tracks[[i]], score_1=tracks[[i]]$score, score_2=tracks[[j]]$score),
    "transparent",
    "transparent",
    score_1 = line_cmap[i],
    score_2 = line_cmap[j],
    roi = tibble(
      as_tibble(peaks),
      sgn = grepl("Earlier", names(peaks)),
      chr = seqnames,
      xmin = start,
      xmax = end,
      ymin = -Inf,
      ymax = Inf,
      fill = c(
        "FALSE" = "#461B56",
        "TRUE" = bg_colors[i]
      )[
        as.character(sgn)
      ]
    ) %>%
      arrange(fill)
  ) %>%
    rbind(
      gtable(width, unit(0.5, "in"), respect = FALSE) %>%
        gtable_add_grob(
          textGrob(title, x = unit(12, "pt"), just = "left", gp = gpar(fontfamily = "Helvetica", fontsize = 24)),
          1,
          1
        ),
      .
    )
  gr <- gr %>%
    grob_col() %>%
    grob_row() %>%
    grob_layout(
      height = as.numeric(convertUnit(reduce(gr$heights, `+`), "mm")),
      width = as.numeric(convertUnit(width, "mm"))
    )
}

grob_to_pdf(
  plot_track_custom(1, 2),
  plot_track_custom(1, 3),
  plot_track_custom(1, 4),
  plot_track_custom(2, 3),
  plot_track_custom(2, 4),
  plot_track_custom(3, 4),
  cbind(
    cowplot::get_legend(
      ggplot(
        tibble(name = c("Germline", "Somatic", "Kc167", "S2") %>% factor(., .)),
        aes(x = 0, y = 0, xend = 1, yend = 1, color = name)
      ) +
        geom_segment() +
        scale_color_manual("Cell Type", values = line_cmap)
    ),
    cowplot::get_legend(
      ggplot(
        tibble(name = c("Germline", "Somatic", "Kc167") %>% factor(., .)),
        aes(x = 0, y = 0, fill = name)
      ) +
        geom_tile() +
        scale_fill_manual("Earlier Cell Type Test", values = bg_colors)
    )
  ) %>%
    grob_col() %>%
    grob_row() %>%
    grob_layout(),
  file_name = "figure/Repli-Timing-Report.pdf"
)
