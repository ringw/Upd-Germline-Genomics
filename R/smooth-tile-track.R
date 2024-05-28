ceil_discrete_log <- function(n) {
  for (i in 1:31) {
    if (bitops::bitShiftL(1, i) >= n) return(i)
  }
  return(32)
}

ksmooth_sliding_windows <- function(granges, bw = 25) {
  stopifnot(length(granges@seqnames@values) == 1)
  granges_widths <- as(granges@ranges@width, "Rle")
  stopifnot(length(granges_widths@lengths) == 2)
  stopifnot(granges_widths@lengths[2] == 1)
  granges_bw <- bw / granges_widths@values[1]
  weights <- as(granges$score, "sparseVector")
  seqlength <- seqlengths(granges)[as.character(granges@seqnames[1])] / granges_widths@values[1]
  dens <- suppressWarnings(
    density(
      weights@i,
      weights = weights@x,
      bw = granges_bw,
      from = 1,
      to = 2^ceil_discrete_log(seqlength),
      n = 2^ceil_discrete_log(seqlength)
    )
  )
  GRanges(granges, score = head(dens$y, length(granges)))
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