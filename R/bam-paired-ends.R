chr_fft_length <- bitwShiftL(1, 25)
chr_filter_size <- floor((chr_fft_length - max(chr.lengths) - 1) / 2)

make_frag_end_filter <- function(bw = 25) {
  time_domain <- c(
    dnorm(seq(0, chr_filter_size), sd = bw),
    rep(0, chr_fft_length - (2 * chr_filter_size + 1)),
    dnorm(seq(-chr_filter_size, -1), sd = bw)
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

bam_paired_fragment_ends <- function(bam_file) {
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
  stopifnot(all(between(read_origin_table$loc, 1, chr.lengths[read_origin_table$chr])))
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

smooth_chr_bp <- function(bp_vec, filter_fft, bin_size = 10) {
  do_smooth <- function(v) {
    len <- length(v)
    v <- c(v %>% as.numeric, rep(0, chr_fft_length - len))
    v <- v %>% fft %>%
      `*`(filter_fft) %>%
      fft(inverse = TRUE) %>%
      `/`(chr_fft_length)
    v <- with_options(list(warn=-1), as.numeric(v))
    v <- v[1:len]

    # Use even steps when taking log. FPKM ~0.008 seems quite small.
    v <- c(
      v[seq(bin_size/2, length(v), by=bin_size)],
      if (length(v) %% bin_size >= floor(bin_size/2))
        NULL
      else v[length(v) - (length(v) %% bin_size)]
    )
    lengths <- rep(
      c(bin_size, len %% bin_size),
      c(
        floor(len / bin_size),
        if (len %% bin_size == 0) 0 else 1
      )
    )
    Rle(
      with_options(list(warn=-1), log(v) / log(2)) %>%
        # FPKM ~0.008 seems quite small
        replace(!is.finite(.) | . < -7, -7) %>%
        round(digits = 2) %>%
        `*`(log(2)) %>%
        exp,
      lengths
    )
  }
  # chr_vec <- Rle(names(chr.lengths), chr.lengths)
  bp_vec %>%
    split(Rle(names(chr.lengths), chr.lengths)) %>%
    sapply(do_smooth, simplify = FALSE) %>%
    RleList
}