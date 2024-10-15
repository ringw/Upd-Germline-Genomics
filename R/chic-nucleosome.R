nucleosome_normal_to_bed <- function(named_normal_files, output_path) {
  data <- tibble(
    bind_rows(
      sapply(
        named_normal_files,
        \(n) read.table(n, comment="<"),
        simplify=F
      ),
      .id = "chr"
    ),
    start = round(V1 - qnorm(0.975) * V2) - 1,
    end = round(V1 + qnorm(0.975) * V2)
  )
  with_options(
    list(scipen=100),
    write.table(
      data[c("chr", "start", "end")],
      output_path,
      quote=F,
      sep="\t",
      row.names=F,
      col.names=F
    )
  )
  output_path
}

bed_to_mark_bigwig <- function(bed_file, bw_path) {
  bed_file <- bed_file %>% read.table %>% split(.$V1)
  granges_list <- list()
  for (data in bed_file) {
    if (nrow(data) == 0) next()
    ir_orig <- IRanges::reduce(IRanges(1 + data$V2, data$V3))
    loc <- c(1, as.numeric(rbind(start(ir_orig), end(ir_orig))))
    ir <- IRanges(
      start = loc[-length(loc)],
      end = loc[-1] - 1
    )
    gr <- GRanges(
      data$V1[1],
      ir,
      score = seq(0, length(ir) - 1) %% 2,
      seqlengths = setNames(loc[length(loc)], data$V1[1])
    )
    granges_list[[data$V1[1]]] <- gr
  }
  export(unlist(GRangesList(granges_list)), BigWigFile(bw_path))
}

# Fix an observed issue: A small genomic feature might have too high of an input
# FPKM (e.g. 1 kb wide). These fragments hinder identifying an average
# nucleosome statistic across a huge window. This is an alternative to NucTools'
# approach to filtering, which capped each fragment from only appearing in 7500
# fragment-fragment pairs (those with the smallest distances). While NucTools'
# cap could distort the repeat length autocorrelation plot, our approach of
# dropping the problematic genomic window will leave a great deal of the huge
# window being looked at untouched, with the modification being made a clean
# excision.
sliding_nucleosome_qc_for_histogram_usability <- function(gr, step=10000, mad_cap=10) {
  bins <- split(
    gr,
    cut(mid(gr), seq(-1+min(mid(gr)), max(mid(gr)), by=step))
  )
  bin_size <- sapply(bins, length)
  bin_size_cap <- median(bin_size) + mad_cap * mad(bin_size)
  unlist(
    bins[bin_size < bin_size_cap],
    use.names = FALSE
  )
}

# Sliding version of NucTools nucleosome_repeat_length.pl.
sliding_nucleosome_repeat_length_analysis <- function(gr, width, step, ...) {
  # We will create sliding windows according to the "step". Then we can resize
  # to the desired width, covering even our chr 4 which is shorter than the
  # width.
  wnd <- slidingWindows(
    GRanges(
      seqnames(gr)[1],
      IRanges(
        1,
        width = seqlengths(gr)[as.character(seqnames(gr)[1])]
      )
    ),
    step,
    step
  )[[1]]
  wnd <- wnd %>%
    GenomicRanges::resize(
      width,
      fix = "center"
    )
  scores <- as_tibble(
    sapply(
      as.character(1:10),
      \(name) integer(0)
    )
  )
  for (i in seq_along(wnd)) {
    scores <- scores %>%
      rbind(
        nucleosome_repeat_length_analysis(
          gr[to(findOverlaps(wnd[i], gr))],
          ...
        ) %>%
          c(rep(NA, 10 - length(.)))
      )
  }
  colnames(scores) <- paste0("nucpair", as.character(1:10))
  elementMetadata(wnd) <- scores
  wnd
}

# Vectorized R version of NucTools nucleosome_repeat_length.pl.
# The computation (1 bp bars of auto-correlated fragment midpoints) is followed
# by our own downstream step, mimicing the NRLcalc project. It is a Gaussian
# smoothing followed by simply finding local maxima. The NRLcalc project
# actually suggests picking out some inflection points (suggestively picking
# points on the plot with prior knowledge that they are evenly spaced). That
# could actually produce additional sequences of autocorrelation peaks which
# look good to the eye. OTOH, the peaks that are selected here will need to pass
# a QC (goodness of fit of the model).
nucleosome_repeat_length_analysis <- function(fragments, b=40) {
  if (length(fragments) < 1000) return(integer(0))
  chr_length <- seqlengths(fragments)[
    as.character(seqnames(fragments)[1])
  ]
  stopifnot(!is.na(chr_length))
  # Length of auto-correlation window to compute.
  delta <- 1500
  # From NucTools. This is a needless computational check. For deep sequencing,
  # we would need a very large delta (histogram size) size in NucTools and then
  # to chop off the range that is not of interest. One of our GSC H3 samples
  # (GC76045515) is deeply sequenced (different technology), and this delta
  # could result in chopping off only 1-3 nucleosome repeats.
  # max_lookup <- 1 + 5 * delta
  max_lookup <- Inf
  # Using NucTools defaults, 1 bed object == 1 paired-end fragment, and compute
  # the midpoint base pair.
  fragments <- fragments[order(mid(fragments))]

  # For each fragment, find the number of neighbors following the fragment that
  # we are going to tabulate.
  num_neighbors <- findInterval(
    mid(fragments) + delta,
    mid(fragments)
  ) -
    seq_along(fragments)
  num_neighbors <- num_neighbors %>% pmin(max_lookup)
  # Concatenate all of the nearby neighbors (to the right) of all of the
  # fragments. These can be compared against the monotonically increasing
  # sequence representing the left-hand fragment that we are looking at.
  neighbor_lookup <- sapply(
    seq_along(num_neighbors),
    \(i) i + seq(0, num_neighbors[i])
  ) %>%
    do.call(c, .)
  tabs <- tabulate(
    abs(
      rep(mid(fragments), 1 + num_neighbors) -
        mid(fragments)[neighbor_lookup]
    ),
    delta
  )
  # New Gaussian smooth step, which looks like the smoothing being applied in
  # the NRLcalc demo.
  signal <- ksmooth(
    seq_along(tabs),
    tabs,
    "normal",
    b = b,
    x.points = seq_along(tabs)
  )$y
  # New finding of local maxima in the smoothed plot (naive).
  points <- 1 +
    which(
      diff(sign(diff(signal))) == -2
    )
  # New filtering. Points should be at least 100 bp in autocorrelation distance
  # apart. We will suppress peaks by picking the one with the greatest y-value.
  points <- points %>%
    subset(
      sapply(
        .,
        \(p) sapply(
          .,
          \(q) abs(p - q) > 100 |
            (signal[p] == signal[q] & p <= q) |
            signal[p] > signal[q]
        ) %>%
          all()
      )
    )
}

# New criteria not found in the NucTools paper. If Pearson's R for the trend
# being fitted is too poor, then we should automatically reject the replicate
# at this locus.
nucleosome_repeat_length_calling_test_points <- function(points) {
  cor(points, 1:5) >= 0.99
}

# Fit the nucleosome model to the autocorrelation local maxima. The
# autocorrelation points were picked out by a heuristic and should pass a QC for
# goodness of fit (Pearson R for the model). This is filtered on a per-track
# basis, then we can pool the data from the tracks (replicates) to find the
# slope parameter for the model (the slope estimates the actual quantity of bp
# of interest).
nucleosome_repeat_length_calling <- function(maxima_tracks) {
  gr <- GRanges(
    seqnames(maxima_tracks[[1]]),
    ranges(maxima_tracks[[1]]),
    seqinfo = seqinfo(maxima_tracks[[1]])
  )
  gr$score <- NA
  gr$n <- NA
  peaks_take <- 2:6
  for (i in seq_along(gr)) {
    maxima_take_granges <- subset(
      maxima_tracks,
      sapply(
        maxima_tracks,
        \(tr) nucleosome_repeat_length_calling_test_points(
          unlist(elementMetadata(tr)[i, peaks_take, drop=T])
        )
      )
    )
    maxima_take <- maxima_take_granges %>%
      sapply(
        \(tr) unlist(elementMetadata(tr)[i, peaks_take, drop=T])
      ) %>%
      as.numeric()
    gr$score[i] <- cor(
      maxima_take,
      rep(peaks_take, length(maxima_take_granges))
    ) *
      sd(maxima_take) /
      sd(peaks_take)
    gr$n[i] <- length(maxima_take_granges)
  }
  gr
}

# Peak calling for NRL nucleosome accessibility track.
nucleosome_repeat_length_peak_analysis <- function(maxima_tracks) {
  single_track_nuc_expected_value <- function(track) {
    track_cor <- elementMetadata(track)[, 2:6] %>%
      t() %>%
      as.matrix() %>%
      cor(2:6) %>%
      as.numeric()
    track <- track[which(track_cor >= 0.99)]
    values <- elementMetadata(track)[, 2:6] %>%
      t() %>%
      as.matrix() %>%
      as.numeric()
    reference <- rep(2:6, length(track))
    cor(values, reference) * sd(values) / sd(reference)
  }
  offset_values <- sapply(
    names(chr.lengths),
    \(chr) sapply(
      maxima_tracks,
      \(gr) single_track_nuc_expected_value(gr[seqnames(gr) == chr])
    )
  )
  gr <- GRanges(
    seqnames(maxima_tracks[[1]]),
    ranges(maxima_tracks[[1]]),
    seqinfo = seqinfo(maxima_tracks[[1]])
  )
  pvalue <- rep(NA, length(maxima_tracks[[1]]))
  for (i in seq(length(pvalue))) {
    values <- sapply(
      maxima_tracks,
      \(track) unlist(elementMetadata(track)[i, 2:6])
    )
    values_take <- which(as.numeric(cor(values, 2:6)) >= 0.99)
    values <- values[, values_take, drop=F]
    offset_values_regression <- rep(
      offset_values[values_take, as.factor(seqnames(maxima_tracks[[1]]))[i]],
      each = 5
    ) *
      (0:4)
    if (length(values_take) > 1)
      fct_values_take_regression <- factor(rep(values_take, each = 5))
    else
      next()
    x <- rep(0:4, length(values_take))
    model <- lm(as.numeric(values) ~ x + fct_values_take_regression + offset(offset_values_regression))
    # Two-Tailed T-Test, followed by writing signs of the peaks.
    # pvalue[i] <- coef(summary(model))["x", "Pr(>|t|)"]
    # if (pvalue[i] < 0.05) {
    #   names[i] <- ifelse(coef(model)["x"] > 0, "+", "-")
    # }
    # One-Tailed T-Test
    pvalue[i] <- pt(
      coef(summary(model))["x", "t value"],
      df = summary(model)$df[match("x", names(coef(model)))],
      lower.tail = FALSE
    )
  }
  gr$pvalue <- pvalue
  gr <- gr[which(pvalue < 0.05)]
  gr
}
