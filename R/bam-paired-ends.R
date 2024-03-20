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