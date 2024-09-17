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

nucleosome_repeat_length_analysis <- function(fragments) {
  # Our NRL could operate on either the set of Watson strand-mapped reads or the
  # set of Crick strand-mapped reads! For now, we will use the starts of the
  # IRanges, which are Watson-mapped.
  chr_length <- seqlengths(fragments)[
    as.character(seqnames(fragments)[1])
  ]
  stopifnot(is.finite(chr_length))
  lst <- fragments %>%
    split(
      cut(
        start(fragments),
        union(
          seq(0, chr_length, by=2000),
          chr_length
        )
      )
    )
  starts <- sapply(lst, start, simplify=FALSE)
  dist_within <- sapply(
    starts,
    \(v) abs(
      rep(v, length(v)) - rep(v, each = length(v))
    ) %>%
      tabulate(2000)
  )
  dist_neighbor <- mapply(
    \(v1, v2) abs(
      rep(v1, length(v2)) - rep(v2, each = length(v1))
    ) %>%
      tabulate(2000),
    starts[-1],
    starts[-length(starts)]
  )
  dist_within_rolling <- dist_within %>% rowCumsums()
  dist_neighbor_rolling <- dist_neighbor %>% rowCumsums()
  est_gr <- slidingWindows(
    GRanges(
      "*",
      1:length(lst)
    ),
    width = 1000
  )[[1]]
  est_gr$score <- sapply(
    seq(length(est_gr)),
    \(i) 149 + which.max(
      dist_within_rolling[-(1:149), end(est_gr[i])] -
        if (start(est_gr)[i] > 1) {
          dist_within_rolling[-(1:149), start(est_gr[i]) - 1]
        } else {
          0
        } +
      dist_neighbor_rolling[-(1:149), end(est_gr[i]) - 1] -
        if (start(est_gr)[i] > 1) {
          dist_neighbor_rolling[-(1:149), start(est_gr[i]) - 1]
        } else {
          0
        }
    )
  )
  GRanges(
    seqnames(fragments)[1],
    IRanges(
      mid(est_gr) * 2000,
      width = 1
    ),
    score = est_gr$score %>%
      replace(. < 155, NA),
    seqinfo = Seqinfo(
      seqnames = as.character(seqnames(fragments)[1]),
      seqlengths = seqlengths(fragments)[
        as.character(seqnames(fragments)[1])
      ]
    )
  )
}

granges_mean_score <- function(ga, gb, gc) {
  data <- cbind(
    A = ga$score,
    B = gb$score,
    C = gc$score
  )
  GRanges(ga, score = colMeans(data, na.rm=T))
}
