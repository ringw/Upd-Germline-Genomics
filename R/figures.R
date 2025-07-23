# ggsave wrapper takes data frame of list of ggplot objects to write ----
save_figures <- function(dir, extension, named_figures, dpi = 120) {
  dir.create(dir, recursive = TRUE, showW = FALSE)
  filenames <- NULL
  for (i in seq(nrow(named_figures))) {
    filename <- paste0(
      dir,
      "/",
      if (type_sum(named_figures %>% pull(1)) == "chr") named_figures[i, 1]
      else rownames(named_figures)[i],
      extension
    )
    ggsave(
      filename,
      named_figures[i, ]$figure[[1]],
      width = named_figures[i, ]$width,
      height = named_figures[i, ]$height,
      dpi = dpi,
      bg = "white",
      device = ifelse(extension == ".pdf", cairo_pdf, png)
    )
    filenames = c(filenames, filename)
  }
  filenames
}

save_png <- function(filename, pl, width, height) {
  png(filename, width=width, height=height)
  print(pl)
  dev.off()
  filename
}