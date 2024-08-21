library(GenomicRanges)
library(magrittr)
library(MatrixGenerics)
library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)
control_path <- args[1]
treatment_path <- args[2]
output_filename <- args[3]
bandwidth <- 100
dir.create(dirname(output_filename), rec=T, showW=F)

control <- import(BigWigFile(control_path))
treatment <- import(BigWigFile(treatment_path))
ceil_discrete_log <- function(n) {
  for (i in 1:31) {
    if (bitops::bitShiftL(1, i) >= n) return(i)
  }
  return(32)
}
smooth_score <- function(gr) {
  if (length(gr) == 1) return(gr)
  fftlength <- 2^ceil_discrete_log(length(gr))
  stopifnot(gr@ranges@start[3] - gr@ranges@start[2] == gr@ranges@start[4] - gr@ranges@start[3])
  bw_to_use <- bandwidth / (gr@ranges@start[3] - gr@ranges@start[2])
  elementMetadata(gr) <- as.data.frame(
    elementMetadata(gr) %>%
      as.matrix %>%
      `*`(as.numeric(rowAlls(. != 0)))
  )
  smooth_result <- apply(
    elementMetadata(gr),
    2,
    \(score) {
      weights <- as(score, "sparseVector")
      dens <- suppressWarnings(
        density(
          weights@i,
          weights = weights@x,
          bw = bw_to_use,
          from = 1,
          to = fftlength,
          n = fftlength
        )
      )
      head(dens$y, length(gr))
    }
  )
  elementMetadata(gr) <- as.data.frame(smooth_result)
  gr
}
smooth_seqs <- function(gr) {
  gr %>%
    split(seqnames(gr)) %>%
    sapply(smooth_score) %>%
    GRangesList %>%
    unlist
}
stopifnot(all.equal(control@ranges, treatment@ranges))
input_data <- GRanges(control@seqnames, control@ranges, seqlengths = seqlengths(control), control = control$score, treatment = treatment$score)
output_data <- smooth_seqs(input_data)
# output_data <- input_data
score <- (log(output_data$treatment) - log(output_data$control)) / log(2)
score <- score - median(score[as.logical(seqnames(input_data) %in% c("2L", "2R", "3L", "3R", "4"))], na.rm=T)
# score <- score %>% replace(
#   output_data$treatment < 0.5 | output_data$control < 0.5,
#   NA
# )
# score <- approx(
#   seq_along(score),
#   score,
#   seq_along(score),
#   rule = 2
# )$y
export(
  GRanges(
    control@seqnames,
    control@ranges,
    seqlengths = seqlengths(control),
    score = score %>% replace(!is.finite(score), 0)
  ),
  BigWigFile(output_filename)
)
