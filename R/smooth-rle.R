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
  size_factor = 1000 * 1000 * 1000 / est_library_size,
  kernel = "gaussian",
  bw = 25,
  sample_size = 10,
  # Heterochromatin ref seqs are smaller in length (6 Mb) compared to the
  # euchromatin chromosome arm seqs (138 Mb), but we did apply different SAM
  # filtering for heterochromatin and for TEs (MAPQ 0 vs 20, -k 10 vs none).
  # The former (MAPQ) entails a decrease in alignment confidence of 100X
  # (likelihood of incorrect alignment goes from 100% to 1%), and for the
  # latter, hardly any reads have a multi-mapping. We sought to roughly equalize
  # mean 2L/2R/3L/3R/X euchromatin FPKM with mean heterochromatin FPKM on the H3
  # input samples by applying a 0.1 multiplier to the heterochromatin counts.
  heterochromatin_mapq_factor = 0.1
) {
  num_counts <- sum(counts)

  if (!is.factor(chrs@values))
    chrs <- Rle(factor(chrs@values, chrs@values), chrs@lengths)

  # First compute the sum features (length 1).
  sum_feature_mask <- Rle(chrs@lengths == 1, chrs@lengths)
  # Now that we can mask the heterochromatin features (sum features, being
  # aggregated into a length-1 Rle), then we can compute the est library size.
  # This estimates the number of reads that we will find for a given molecule
  # in the heterochromatin filtering params (multi-mapping 10 alternate
  # alignments, any MAPQ) vs euchromatin (one alignment, MAPQ 20).
  est_library_size <- sum(counts) - (1 - heterochromatin_mapq_factor) * sum(counts[which(sum_feature_mask)])
  sum_features <- counts[which(sum_feature_mask)] %>%
    as.numeric %>%
    `*`(size_factor) %>%
    setNames(chrs[sum_feature_mask]) %>%
    sapply(\(value) Rle(value, lengths = 1))

  # Second compute the smooth features (euchromatin features) which are not of
  # length 1.
  eu_features <- levels(chrs@values) %>%
    setdiff(names(sum_features)) %>%
    sapply(
      \(n) {
        v <- counts[as.logical(chrs == n)]
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
    )

  all_feature_names <- levels(chrs@values)
  append(eu_features, sum_features)[all_feature_names] %>%
    RleList %>%
    set_attr("est_library_size", est_library_size)
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