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

# For the sample (bam file), plot fragment ends (FE), an estimate of the
# genomic alignment of each end of each molecule in the sample.
# feature.lengths must be named with the reference sequences that exist in the
# bam file. If a length of 1 is substituted for the actual reference length,
# then all of the fragment ends for the reference (an even number, assuming that
# we applied the filter -f 0x2) will be summed up in a scalar. The starts of
# each track or scalar value in the resulting sparseVector are given by:
#   c(1, 1 + cumsum(feature.lengths))
bam_paired_fragment_ends <- function(bam_file, feature.lengths) {
  read_origin_table <- names(feature.lengths) %>%
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
          if (feature.lengths[n] != 1)
            data.frame(
              loc = values,
              count = lengths
            )
          else
            data.frame(
              loc = 1,
              count = sum(values * lengths)
            )
        ),
      simplify = F
    ) %>%
    bind_rows(.id = "chr")
  # Zero-based offset for each chr.
  feature.starts <- setNames(
    c(0, cumsum(feature.lengths[-length(feature.lengths)])),
    names(feature.lengths)
  )
  bam_track <- sparseVector(
    i = feature.starts[read_origin_table$chr] + read_origin_table$loc,
    x = read_origin_table$count,
    length = sum(feature.lengths)
  )
  bam_track
}