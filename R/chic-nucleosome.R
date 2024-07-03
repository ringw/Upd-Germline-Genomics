nucleosome_normal_to_bed <- function(named_normal_files, output_path) {
  data <- tibble(
    bind_rows(
      sapply(
        named_normal_files,
        \(n) read.table(n, comment="<"),
        simplify=F
      ),
      .id = "chr"
    ),
    start = round(V1 - qnorm(0.975) * V2) - 1,
    end = round(V1 + qnorm(0.975) * V2)
  )
  with_options(
    list(scipen=100),
    write.table(
      data[c("chr", "start", "end")],
      output_path,
      quote=F,
      sep="\t",
      row.names=F,
      col.names=F
    )
  )
  output_path
}

bed_to_mark_bigwig <- function(bed_file, bw_path) {
  bed_file <- bed_file %>% read.table %>% split(.$V1)
  granges_list <- list()
  for (data in bed_file) {
    if (nrow(data) == 0) next()
    ir_orig <- IRanges::reduce(IRanges(1 + data$V2, data$V3))
    loc <- c(1, as.numeric(rbind(start(ir_orig), end(ir_orig))))
    ir <- IRanges(
      start = loc[-length(loc)],
      end = loc[-1] - 1
    )
    gr <- GRanges(
      data$V1[1],
      ir,
      score = seq(0, length(ir) - 1) %% 2,
      seqlengths = setNames(loc[length(loc)], data$V1[1])
    )
    granges_list[[data$V1[1]]] <- gr
  }
  export(unlist(GRangesList(granges_list)), BigWigFile(bw_path))
}