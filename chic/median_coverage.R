library(callr)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(stringr)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)
input_id <- args[1]
output_filename <- args[2]
dir.create(dirname(output_filename), rec=T, showW=F)

chr.lengths <- c(
  `2L`=23513712, `2R`=25286936, `3L`=28110227, `3R`=32079331, `4`=1348131, X=23542271, Y=3667352
)

read_chr_granges <- function(input_id, chr) {
  library(magrittr)
  chr.lengths <- c(
    `2L`=23513712, `2R`=25286936, `3L`=28110227, `3R`=32079331, `4`=1348131, X=23542271, Y=3667352
  )
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

  bam_reference_coords <- function(df) {
    df %>%
      dplyr::mutate(
        cigar = factor(cigar),
        width = cigar %>%
          forcats::fct_relabel(\(v) as.character(cigar_ref_length(v))) %>%
          as.character %>%
          as.numeric,
        pos_crick = pos + width - 1
      )
  }

  paired_end_reads_to_fragment_lengths <- function(df) {
    stopifnot(all(df$strand[seq(2, nrow(df), by=2)] == "-"))
    df$fragment_end_crick <- rep(
      bam_reference_coords(df[seq(2, nrow(df), by=2),, drop=F])$pos_crick,
      each=2
    )
    df <- df[seq(1, nrow(df), by=2),, drop=F]
    df$length <- df$fragment_end_crick - df$pos + 1
    df
  }

  paired_end_reads_to_granges <- function(df, min_length = -Inf, max_length = Inf, ...) {
    df <- df %>% paired_end_reads_to_fragment_lengths
    df <- df %>% subset(dplyr::between(length, min_length, max_length))
    GenomicRanges::GRanges(df$rname, IRanges::IRanges(start = df$pos, width = df$length), ...)
  }
  arrow::read_parquet(stringr::str_glue("_targets/objects/bulk_reads_{chr}_chic.bam_{input_id}_chr")) %>%
    # subset(mapq >= 20) %>%
    # paired_end_reads_to_granges(seqlengths = chr.lengths, min_length = 150, max_length = 200)
    paired_end_reads_to_granges(seqlengths = chr.lengths)
}

frag <- sapply(
  c("2L", "2R", "3L", "3R", "4", "X", "Y"),
  \(chr) r(\(read_chr_granges, input_id, chr) read_chr_granges(input_id, chr), args = list(read_chr_granges, input_id, chr))
) %>%
  GRangesList
wnd_coverage <- slidingWindows(
  with(enframe(chr.lengths[1:5]), GRanges(str_glue("{name}:1-{value}"))),
  1000L,
  100L
) %>%
  mapply(
    \(wnd, frag) GRanges(wnd, score = countOverlaps(wnd, frag), seqlengths = seqlengths(frag)),
    .,
    frag[1:5]
  ) %>%
  GRangesList %>%
  unlist

step_size <- 50L
track <- slidingWindows(
  with(enframe(chr.lengths), GRanges(str_glue("{name}:1-{value}"))),
  step_size,
  step_size
) %>%
  mapply(
    \(wnd, frag) GRanges(wnd, score = countOverlaps(wnd, frag), seqlengths = seqlengths(frag)),
    .,
    frag
  ) %>%
  GRangesList %>%
  unlist
track$score <- track$score / (median(wnd_coverage$score) * step_size / 1000)
export(track, BigWigFile(output_filename))
