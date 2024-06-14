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
  # the entire sequence.
  stopifnot(granges@ranges@start[3] - granges@ranges@start[2] == granges@ranges@start[4] - granges@ranges@start[3])
  granges_bw <- bw / (granges@ranges@start[3] - granges@ranges@start[2])
  fftlength <- 2^ceil_discrete_log(length(granges))
  smooth_result <- apply(
    elementMetadata(granges) %>%
      as.data.frame %>%
      list(
        .,
        # Do not smooth "R" factor (replicate) as log-link batch effect can be
        # a -Infinity or Infinity.
        select = grep("^score.R", colnames(.), val=T, invert=T)
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

density_est_granges_exact <- function(granges, obs_vector, bw, sample_rate) {
  est_library_size <- granges@metadata$est_library_size
  stopifnot(is.finite(est_library_size))
  weight_mult <- 1000^3 / est_library_size
  dens <- suppressWarnings(
    density(
      obs_vector@i,
      weights = obs_vector@x * weight_mult,
      bw = bw,
      from = 1,
      to = 2^ceil_discrete_log(length(obs_vector) / sample_rate),
      n = 2^(1 + ceil_discrete_log(length(obs_vector) / sample_rate))
    )
  )
  x_start <- seq(1, length(obs_vector), by=sample_rate)
  x_width <- rep(sample_rate, length(x_start))
  x_width[length(x_width)] <- length(obs_vector) - x_start[length(x_width)] + 1
  y <- dens$y[head(seq(2, length(dens$y), by=2), length(x_start))]
  gr <- GRanges(
    granges@seqnames,
    ranges = IRanges(start = x_start, width = x_width),
    score = y,
    seqlengths = seqlengths(granges)
  )
  gr@metadata <- granges@metadata
  gr
}

density_est_granges_approx <- function(granges, obs_vector, bw, sample_rate) {
  est_library_size <- granges@metadata$est_library_size
  stopifnot(is.finite(est_library_size))
  weight_mult <- 1000^3 / est_library_size
  dens <- suppressWarnings(
    density(
      obs_vector@i,
      weights = obs_vector@x * weight_mult,
      bw = bw,
      from = 1,
      to = length(obs_vector) + 1,
      n = 2^(ceil_discrete_log(length(obs_vector) / bw))
    )
  )
  x_start <- seq(1, length(obs_vector), by=sample_rate)
  x_width <- rep(sample_rate, length(x_start))
  x_width[length(x_width)] <- length(obs_vector) - x_start[length(x_width)] + 1

  x_out <- x_start + 1/2*x_width
  y <- approx(dens$x, dens$y, xout = x_out)$y
  gr <- GRanges(
    granges@seqnames,
    ranges = IRanges(start = x_start, width = x_width),
    score = y,
    seqlengths = seqlengths(granges)
  )
  gr@metadata <- granges@metadata
  gr
}