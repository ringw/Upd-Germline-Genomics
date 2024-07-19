ceil_discrete_log <- function(n) {
  for (i in 1:31) {
    if (bitops::bitShiftL(1, i) >= n) return(i)
  }
  return(32)
}

ksmooth_sliding_windows <- function(granges, bw = 25) {
  stopifnot(length(granges@seqnames@values) == 1)

  # No smoothing case.
  if (length(granges) == 1) return(granges)

  # Step size should be uniform. We have added a custom 5' and 3' end window
  # because we later need to create a bigwig track where the intervals partition
  # the entire sequence. However, because we used
  # slidingWindows/resize/restrict, don't check the very first windows (on the
  # first chromosome, which is expected to be very long), but check a later
  # window well past the resized window size/step size that we are using.
  stopifnot(granges@ranges@start[10] - granges@ranges@start[9] == granges@ranges@start[11] - granges@ranges@start[10])
  granges_bw <- bw / (granges@ranges@start[11] - granges@ranges@start[10])
  fftlength <- 2^ceil_discrete_log(length(granges))
  smooth_result <- apply(
    elementMetadata(granges) %>%
      as.data.frame %>%
      list(
        .,
        # Do not smooth "R" factor (replicate) as log-link batch effect can be
        # a -Infinity or Infinity.
        select = grep("^score(.[^R]|$)", colnames(.), val=T)
      ) %>%
      do.call(subset, .),
    2,
    \(score) {
      weights <- as(score, "sparseVector")
      dens <- suppressWarnings(
        density(
          weights@i,
          weights = weights@x,
          bw = granges_bw,
          from = 1,
          to = fftlength,
          n = fftlength
        )
      )
      head(dens$y, length(granges))
    }
  )
  elementMetadata(granges) <- as.data.frame(smooth_result)
  granges
}