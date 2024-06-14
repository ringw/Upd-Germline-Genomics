library(dplyr)
library(future)
library(processx)
library(rtracklayer)
library(stringr)

plan(multicore, workers=2)

args <- commandArgs(trailingOnly = TRUE)
control_path <- args[1]
treatment_path <- args[2]
output_filename <- args[3]

bam_coverage <- function(bam_file) {
  bw_file <- bam_file %>%
    str_replace("(chr|masked)", "deepTools_\\1") %>%
    str_replace("\\.bam", ".bw")
  run(
    "bamCoverage",
    c(
      "--binSize",
      "20",
      "--smoothLength",
      "200",
      "--minFragmentLength",
      "100",
      "--maxFragmentLength",
      "200",
      "--extendReads",
      "200",
      "-b",
      bam_file,
      "-o",
      bw_file
    )
  )
  bw_file
}

control_bw %<-% (control_path %>% bam_coverage)
treatment_bw %<-% (treatment_path %>% bam_coverage)

control_track <- import(BigWigFile(control_bw))
treatment_track <- import(BigWigFile(treatment_bw))
tile_track <- slidingWindows(
  GRanges(
    names(seqlengths(control_track)),
    IRanges(1, width=seqlengths(control_track)),
    seqlengths = seqlengths(control_track)
  ),
  20L,
  20L
) %>%
  unlist
control_track <- GRanges(
  tile_track,
  score = findOverlaps(tile_track, control_track) %>%
    as("List") %>%
    unlist %>%
    `[`(control_track$score, .)
)
treatment_track <- GRanges(
  tile_track,
  score = findOverlaps(tile_track, treatment_track) %>%
    as("List") %>%
    unlist %>%
    `[`(treatment_track$score, .)
)

fc_track <- GRanges(
  tile_track,
  score = (log(treatment_track$score) - log(control_track$score)) %>%
    `*`(1/log(2)) %>%
    replace(!is.finite(.), NA)
)
fc_track$score <- fc_track$score %>%
  `-`(
    median(
      .[as.numeric(seqnames(control_track)) %>% between(1, 4)],
      na.rm=T
    )
  ) %>%
  replace(is.na(.), 0)
export(fc_track, BigWigFile(output_filename))
