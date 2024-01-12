sample_bams <- function(bam_list, sample_rate=0.01) {
  sample_arg <- with_options(list(scipen=100), as.character(sample_rate))
}

read_pileup_df <- function(bam_file, bin_size=10, incl_flags=NULL) {
  bin_index_from_zero <- floor(bin_size / 2)
  df <- NULL
  for (chr in names(chr.lengths)) {
    p <- processx::process$new(
      "sh",
      c(
        "-c",
        paste(
          "samtools",
          "mpileup",
          if (!is.null(incl_flags)) "--incl-flags" else "",
          if (!is.null(incl_flags)) incl_flags else "",
          "-r",
          chr,
          "--output-extra",
          "MAPQ",
          bam_file,
          "|",
          "perl",
          "-ne",
          "'",
          paste0("@f=split; print if ($f[1] % ", bin_size, ") == ", bin_index_from_zero),
          "'"
        )
      ),
      stdout = "|"
    )
    chr_result <- read.table(
      textConnection(p$read_all_output()),
      sep = "\t",
      quote = "",
      comment = "",
      col.names = c("chr", "pos", "ref", "count", "seq", "baseq", "mapq")
    )
    stopifnot(p$get_exit_status() == 0)
    df <- rbind(df, chr_result)
  }
  df
}

mapq_table_strand_chars <- list(
  any = ".ACGTN,acgtn",
  A = ".ACGTN", # forward
  a = ",acgtn" # reverse
)

make_coverage_rowData <- function(resolution) {
  mapply(
    \(name, len) data.frame(
      chr = name,
      pos = seq(floor(resolution/2), len, by = resolution)
    ),
    names(chr.lengths),
    chr.lengths,
    USE.NAMES = FALSE,
    SIMPLIFY = FALSE
  ) %>%
    do.call(rbind, .)
}

mapq_table_strand <- function(pileup, strand = c("any", "A", "a")) {
  strand <- match.arg(strand)
  match_chars <- utf8ToInt(mapq_table_strand_chars[[strand]])
  all_chars <- utf8ToInt(mapq_table_strand_chars$any)

  pileup_res <- pileup %>%
    split(pileup$chr) %>%
    sapply(\(df) min(diff(df$pos))) %>%
    min

  new_rowData <- make_coverage_rowData(pileup_res)

  mapq_max = 42
  mapq_ascii = (seq(0, mapq_max) + 33) %>% sapply(intToUtf8)
  result <- sparseMatrix(
    i = NULL,
    p = 0,
    dims = c(0, length(mapq_ascii)),
    dimnames = list(NULL, mapq_ascii),
    repr = "R"
  )
  
  for (rows in seq(nrow(new_rowData)) %>% split(cut(., seq(min(.)-1, max(.), length.out = 20)))) {
    pileup_batch <- new_rowData[rows, ] %>% left_join(
      pileup,
      by = c("chr", "pos")
    )
    seq_logical <- pileup_batch$seq %>%
      sapply(
        \(seq) seq %>%
          utf8ToInt %>%
          subset(. %in% all_chars) %>%
          `%in%`(match_chars),
        simplify = FALSE,
        USE.NAMES = FALSE
      )
    mapq_vec <- mapply(
      \(seq, mask) seq %>%
        utf8ToInt %>%
        subset(mask) %>%
        `-`(33) %>%
        sort %>%
        as("Rle"),
      pileup_batch$mapq,
      seq_logical,
      SIMPLIFY = FALSE,
      USE.NAMES = FALSE
    )
    # p comes from the number of nonzero bin entries in each row.
    mapq_p <- c(
      0,
      mapq_vec %>%
        sapply(\(rle) length(rle@values)) %>%
        cumsum
    )
    mat_batch <- sparseMatrix(
      p = mapq_p,
      j = do.call(
        c,
        mapq_vec %>% sapply(\(rle) 1 + rle@values, simplify = FALSE)
      ),
      # Now perform a bin count by storing the length (count of value).
      x = do.call(
        c,
        mapq_vec %>% sapply(\(rle) rle@lengths, simplify = FALSE)
      ),
      repr = "R",
      dims = c(length(rows), length(mapq_ascii)),
      dimnames = list(NULL, mapq_ascii)
    )
    result <- result %>% rbind(mat_batch)
  }

  result
}
