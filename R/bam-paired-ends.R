make_frag_end_filter <- function(fft_length, filter_size, bw = 25) {
  time_domain <- c(
    dnorm(seq(0, filter_size), sd = bw),
    rep(0, fft_length - (2 * filter_size + 1)),
    dnorm(seq(-filter_size, -1), sd = bw)
  )
  fft(time_domain)
}

cigar_ref_length <- function(cigars) {
  cigars %>%
    sapply(
      \(cigar) cigar %>%
        gregexpr("[0-9]+[MNDI]", .) %>%
        `[[`(1) %>%
        mapply(
          # Insertion is the only case where we are not counting reference base pairs.
          \(pos, length) if(substring(cigar, pos + length - 1, pos + length - 1) != "I")
            as.numeric(
              substring(cigar, pos, pos + length - 2)
            )
            else 0,
          .,
          attr(., "match.length")
        ) %>%
        sum
    )
}

bam_paired_fragment_ends <- function(bam_file, chr.lengths) {
  read_origin_table <- names(chr.lengths) %>%
    sapply(
      \(n) run(
        "bash",
        c(
          "-c",
          paste0(
            "samtools view -f 0x2 ",
            bam_file,
            " ",
            n,
            " | cut -f 4,6,9",
            " | perl -MList::Util=sum -ne",
            # Cigar: M,N,D are a ref base pair--I is a query base pair insertion
            " '@m=/([0-9]+)[MND]/g;",
            " chomp; $s=0+sum(@m); @p=split;",
            " print \"$p[0]\t$s\t$p[2]\n\"'"
          )
        )
      )$stdout %>%
        textConnection %>%
        read.table(header = F, quote = "", col.names = c("pos", "cigar_length", "tlen")) %>%
        as.list %>%
        with(
          ifelse(
            tlen < 0,
            pos + cigar_length - 1,
            pos
          )
        ) %>%
        sort %>%
        rle %>%
        with(
          data.frame(
            loc = values,
            count = lengths
          )
        ),
      simplify = F
    ) %>%
    bind_rows(.id = "chr")
  # Reads mapped to the chromosomes did not pass beyond the end of the
  # chromosome, so this check is useful for validation there. On the
  # transposable elements, the 3' end in alignment coordinates might
  # have extended beyond the end.
  stopifnot(
    all(
      between(read_origin_table$loc, 1, chr.lengths[read_origin_table$chr])
      | str_starts(read_origin_table$chr, "FBte")
    )
  )
  # Zero-based offset for each chr.
  chr.starts <- setNames(
    c(0, cumsum(chr.lengths[-length(chr.lengths)])),
    names(chr.lengths)
  )
  bam_track <- sparseVector(
    i = chr.starts[read_origin_table$chr] + read_origin_table$loc,
    x = read_origin_table$count,
    length = sum(chr.lengths)
  )
  bam_track
}

fpkm_aggregate <- function(vecs) {
  # FPKM calculation. Use a 1 bp window.
  vecs <- vecs %>% sapply(
    \(v) v * (1000 * 1000 * 1000 / sum(v)),
    simplify = F
  )
  purrr::reduce(vecs, `+`) / length(vecs)
}

fpkm_cbind <- function(vecs) {
  # FPKM calculation. Use a 1 bp window.
  vecs <- vecs %>% sapply(
    \(v) as(
      v * (1000 * 1000 * 1000 / sum(v)),
      "sparseMatrix"
    ),
    simplify = F
  )
  for (v in names(vecs)) colnames(vecs[[v]]) <- v
  purrr::reduce(vecs, cbind)
}

cbind_bulk_data <- function(data) {
  # colData in first n-1 columns. Final column contains sparseVector.
  m <- do.call(
    cbind,
    data$track %>%
      sapply(\(v) as(v, "sparseMatrix"), simplify = FALSE)
  )
  SingleCellExperiment(
    list(fpkm=m),
    colData = data %>% subset(select = -track) %>% as.data.frame
  )
}

smooth_chr_bp <- function(bp_vec, chic.kde.filter, chr.lengths, bin_size = 10) {
  do_smooth <- function(v) {
    len <- length(v)
    fft_to_use <- len %>% findInterval(c(0, chic.kde.filter$limit_length))
    fft_length <- pull(chic.kde.filter, "fft_length")[fft_to_use]
    filter_fft <- pull(chic.kde.filter, "filter")[fft_to_use][[1]]

    v <- c(v %>% as.numeric, rep(0, fft_length - len))
    v <- v %>% fft %>%
      `*`(filter_fft) %>%
      fft(inverse = TRUE) %>%
      `/`(fft_length)
    v <- with_options(list(warn=-1), as.numeric(v))
    v <- v[1:len]

    v <- c(
      # Nice sampling of the density at each base pair, giving us the estimate at
      # the midpoint of the range in our Rle.
      v[seq(bin_size/2, length(v), by=bin_size)],
      # If the final window is less than half-full, then we need to sample an
      # extra value at the final window.
      if (between(length(v) %% bin_size, 1, floor(bin_size/2) - 1))
        v[length(v) - (length(v) %% bin_size) + 1]
      else NULL
    )
    lengths <- rep(
      c(bin_size, len %% bin_size),
      c(
        floor(len / bin_size),
        if (len %% bin_size == 0) 0 else 1
      )
    )
    res = Rle(
      with_options(list(warn=-1), log(v) / log(2)) %>%
        # FPKM ~0.008 seems quite small
        replace(!is.finite(.) | . < -7, -7) %>%
        round(digits = 2) %>%
        `*`(log(2)) %>%
        exp,
      lengths
    )
    attr(res, "standard_deviation") <- "fake"
    res
  }
  # chr_vec <- Rle(names(chr.lengths), chr.lengths)
  bp_vec %>%
    split(Rle(names(chr.lengths), chr.lengths)) %>%
    sapply(do_smooth, simplify = FALSE) %>%
    RleList
}

smooth_average_cols <- function(
  rep_matrix, chic.kde.filter, chr.lengths, bin_size = 10
) {
  do_smooth_rowMean <- function(m) {
    len <- nrow(m)
    fft_to_use <- len %>% findInterval(c(0, chic.kde.filter$limit_length))
    fft_length <- pull(chic.kde.filter, "fft_length")[fft_to_use]
    filter_fft <- pull(chic.kde.filter, "filter")[fft_to_use][[1]]

    m <- rbind(
      m %>% as.matrix, matrix(0, nrow = fft_length - len, ncol = ncol(m))
    )
    m <- m %>% mvfft %>%
      `*`(filter_fft) %>%
      mvfft(inverse = TRUE) %>%
      `/`(fft_length)
    m <- Re(m)
    m <- m[1:len, ]

    sample_loc <- c(
      # Nice sampling of the density at each base pair, giving us the estimate at
      # the midpoint of the range in our Rle.
      seq(bin_size/2, len, by=bin_size),
      # If the final window is less than half-full, then we need to sample an
      # extra value at the final window.
      if (between(len %% bin_size, 1, floor(bin_size/2) - 1))
        len - (len %% bin_size) + 1
      else NULL
    )
    m <- m[sample_loc, ]
    lengths <- rep(
      c(bin_size, len %% bin_size),
      c(
        floor(len / bin_size),
        if (len %% bin_size == 0) 0 else 1
      )
    )
    res = Rle(
      with_options(list(warn=-1), log(rowMeans(m)) / log(2)) %>%
        # FPKM ~0.008 seems quite small
        replace(!is.finite(.) | . < -7, -7) %>%
        round(digits = 2) %>%
        `*`(log(2)) %>%
        exp,
      lengths
    )
    attr(res, "standard_deviation") <- Rle(
      (log(rowSds(m)) / log(2)) %>%
        round(digits = 2) %>%
        `*`(log(2)) %>%
        exp,
      lengths
    )
    res
  }
  row_data <- Rle(names(chr.lengths), chr.lengths)
  reps_blocks <- sapply(
    names(chr.lengths),
    # c("Y","FBte0000015"),
    # names(chr.lengths)[-(1:7)],
    \(n) rep_matrix[as.logical(row_data == n), ],
    simplify = FALSE
  )
  smooth_rles <- reps_blocks %>%
    sapply(do_smooth_rowMean, simplify = FALSE)
  res <- RleList(smooth_rles)
  attr(res, "standard_deviation") <- mapply(
    \(n, obj) attr(obj, "standard_deviation"),
    names(smooth_rles),
    smooth_rles,
    SIMPLIFY = FALSE
  ) %>%
    RleList
  attr(res, "n") <- ncol(rep_matrix)
  res
}