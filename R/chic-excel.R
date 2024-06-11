publish_chic_fragments <- function(
  frag_size_tables_publish,
  mapq_tables_publish,
  output_path
) {
  wb <- createWorkbook()

  title <- "H3 ChIC"
  addWorksheet(wb, title)
  start_row <- 1
  writeData(wb, title, "Fragment Sizes", colNames=F, startRow = start_row, startCol = 1)
  start_row <- start_row + 1
  for (n in names(frag_size_tables_publish)) {
    writeData(wb, title, n, colNames=F, startRow = start_row, startCol = 1)
    frag_size_tables_publish[[n]]$fragment_sizes <- frag_size_tables_publish[[n]]$fragment_sizes %>%
      factor(unique(.))
    tbl <- dcast(frag_size_tables_publish[[n]], chr ~ fragment_sizes, fun=sum, value.var="n")
    writeDataTable(
      wb,
      title,
      tbl,
      startRow = start_row + 1, startCol = 1
    )
    for (iloc in seq(nrow(tbl)))
      tbl[iloc, -1] <- tbl[iloc, -1] %>% `/`(sum(.)) %>% round(3)
    writeDataTable(
      wb,
      title,
      tbl,
      startRow = start_row + 1, startCol = 1 + ncol(tbl) + 1
    )
    addStyle(
      wb, title,
      createStyle(numFmt = "0.0%"),
      rows = start_row + 1 + seq(nrow(tbl)),
      cols = seq(1 + ncol(tbl) + 1 + 1, 1 + ncol(tbl) + ncol(tbl)),
      gridExpand=TRUE)
    start_row <- start_row + 1 + 1 + nrow(tbl) + 1
  }

  writeData(wb, title, "MAPQ Values", colNames=F, startRow = start_row, startCol = 1)
  start_row <- start_row + 1
  for (n in names(mapq_tables_publish)) {
    writeData(wb, title, n, colNames=F, startRow = start_row, startCol = 1)
    tbl <- dcast(mapq_tables_publish[[n]], chr ~ mapq_range, fun=sum, value.var="n")
    writeDataTable(
      wb,
      title,
      tbl,
      startRow = start_row + 1, startCol = 1
    )
    for (iloc in seq(nrow(tbl)))
      tbl[iloc, -1] <- tbl[iloc, -1] %>% `/`(sum(.)) %>% round(3)
    writeDataTable(
      wb,
      title,
      tbl,
      startRow = start_row + 1, startCol = 1 + ncol(tbl) + 1
    )
    addStyle(
      wb, title,
      createStyle(numFmt = "0.0%"),
      rows = start_row + 1 + seq(nrow(tbl)),
      cols = seq(1 + ncol(tbl) + 1 + 1, 1 + ncol(tbl) + ncol(tbl)),
      gridExpand=TRUE)
    start_row <- start_row + 1 + 1 + nrow(tbl) + 1
  }

  dir.create(dirname(output_path), showW=F, rec=F)
  saveWorkbook(wb, output_path, overwrite=T)
  output_path <- output_path
}