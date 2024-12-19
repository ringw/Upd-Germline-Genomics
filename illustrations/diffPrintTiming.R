source("_targets.R")
sapply(tar_option_get("packages"), \(n) require(n, ch=T))
library(grid)
library(grobblR)
library(gtable)
library(purrr)

tar_load(matches("repli.timing_"))

tracks <- list(
  repli.timing_Germline_chr,
  repli.timing_Somatic_chr,
  repli.timing_Kc167_chr,
  repli.timing_S2_chr
)
line_cmap <- c("#6D9965", "#B377B0", "#C27A00", "#4087F4")

plot_track_custom <- function(i, title) {
  width <- unit(6, "in")
  gr <- gtable(
    unit(1, "null"),
    unit(c(0.5, 4.5), "in")
  ) %>%
    gtable_add_grob(
      list(
        textGrob(title, x = unit(12, "pt"), just = "left", gp = gpar(fontfamily = "Helvetica", fontsize = 24)),
        plot_track(
          tracks[[i]],
          "transparent",
          "transparent",
          plot_grid = FALSE
        )
      ),
      t = 1:2,
      l = 1
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
  plot_track_custom(1, "Germline"),
  plot_track_custom(2, "Somatic"),
  plot_track_custom(3, "Kc167"),
  plot_track_custom(4, "S2"),
  file_name = "figure/Repli-Timing-Each-Track.pdf"
)
