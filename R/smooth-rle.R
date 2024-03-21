# Use R's "density" for a rolling (rectangular) function. The rectangular filter
# (which follows the shape of the uniform probability distribution) has the s.d.
# 1/sqrt(12), and the distribution may be scaled by a diameter other than 1.
rolling_filter_to_bw <- function(rolling_filter_diameter) {
  rolling_filter_diameter * 1 / sqrt(12)
}

# Smooth fragment ends using density estimation.
smooth_sparse_vector_to_rle_list <- function(
  counts,
  chrs,
  size_factor = 1000 * 1000 * 1000 / sum(counts),
  kernel = "gaussian",
  bw = 25,
  sample_size = 10
) {
  if (!is.factor(chrs@values))
    chrs <- Rle(factor(chrs@values, chrs@values), chrs@lengths)
  levels(chrs@values) %>%
    sapply(
      \(n) {
        v <- counts[as.logical(chrs == n)]
        # Never apply density() to a single entry. We already have a discrete
        # scale (base pairs) so the density estimate in the single bin is the
        # original value in the single bin. Size factor is included as a
        # multiplier for the scale that each counted event is on (e.g. FPKM).
        if (length(v) == 1)
          return(Rle(as.numeric(v[1]) * size_factor, lengths = 1))
        which_v <- rep(v@i, times = v@x)
        length_log_2 <- log(length(v) / sample_size) / log(2)
        density_n <- round(exp(ceiling(length_log_2) * log(2)))
        # Default weight is 1/length. We do not want to create infinitesimally
        # small values over the length of the chromosome which is very long.
        # However, density() is always going to warn about the fact that we are
        # not introducing the 1/n weights which sum to 1.
        dens <- with_options(
          list(warn=-1),
          density(
            which_v,
            from=1,
            to=length(v),
            n=density_n,
            kernel = kernel,
            bw = bw,
            weights = rep(size_factor, length(which_v))
          )
        )
        density_sampling = tibble(
          begin = seq(0, length(v), by = sample_size),
          end = pmin(begin + sample_size, length(v)),
          center = (begin + end) / 2,
          dens = approx(dens$x, dens$y, xout=pmax(1, center))$y
        )
        with(density_sampling, Rle(dens, end - begin))
      }
    ) %>%
    RleList
}

smooth_sparse_vector_macs <- function(
  counts,
  chrs,
  diameter,
  sample_size,
  normalize_tag_count = TRUE
) {
  if (!is.factor(chrs@values))
    chrs <- Rle(factor(chrs@values, chrs@values), chrs@lengths)

  # Suppose that the diameter is nearly the size of the transposon, then we will
  # pretend that we used a sliding sum with the given diameter, but we actually
  # used the chr-wide mean. Multiply by diameter to emulate the sliding sum.
  macs_chr_background <- counts %>%
    split(chrs) %>%
    sapply(mean) %>%
    `*`(diameter) %>%
    floor
  macs_chr_background <- macs_chr_background %>%
    mapply(
      \(macs_chr_background, len) Rle(macs_chr_background, len),
      .,
      chrs@lengths
    ) %>%
    RleList
  smooth_chr <- smooth_sparse_vector_to_rle_list(
    counts,
    chrs,
    size_factor = diameter,
    kernel = "rectangular",
    bw = rolling_filter_to_bw(diameter),
    sample_size = sample_size
  )
  smooth_chr <- RleList(
    mapply(
      \(len, bg_track, smooth_track)
        if (diameter * 4 < len) smooth_track
        else bg_track,
      setNames(chrs@lengths, chrs@values),
      macs_chr_background,
      smooth_chr
    )
  )
  if (normalize_tag_count)
    smooth_chr <- smooth_chr * 1000 * 1000 * 1000 / sum(counts) / diameter
  attr(smooth_chr, "num_tags") <- sum(counts)
  attr(smooth_chr, "diameter") <- diameter
  smooth_chr
}