fft_length <- bitwShiftL(1, 25)

cross_correlate_mapq <- function(bams, min_mapq = 20, n_corr = 500) {
  min_mapq_index <- min_mapq + 1
  signals <- bams %>% sapply(\(m) rowSums(m[, seq(min_mapq_index, ncol(m)), drop = F]))
  ffts <- signals %>%
    apply(
      2,
      \(v) fft(
        c(v, rep(0, fft_length - length(v)))
      )
    )
  corr <- fft(
    ffts[, 1] * Conj(ffts[, 2]),
    inverse = TRUE
  )
  corr <- with_options(list(warn = -1), as.numeric(corr))
  corr <- corr / length(corr) / prod(colSds(signals, center = c(0, 0))) / nrow(signals)
  rbind(
    data.frame(
      shift_forward = seq(-n_corr, -1),
      pearson = corr[(length(corr) - n_corr + 1):length(corr)]
    ),
    data.frame(
      shift_forward = 0:n_corr,
      pearson = corr[1:(n_corr + 1)]
    )
  )
}

cross_correlate_mapq2 <- function(bams, resolution = 10, min_mapq = 20, n_corr = 500) {
  fft_length <- bitwShiftL(1, 23)
  min_mapq_index <- min_mapq + 1
  signals <- bams %>% sapply(\(m) rowSums(m[, seq(min_mapq_index, ncol(m)), drop = F]))
  signals <- signals %*% diag(
    x = 1 / colSds(signals, center = rep(0, ncol(signals)))
    / sqrt(nrow(signals))
  )
  rowData <- make_coverage_rowData(resolution)
  stopifnot(nrow(rowData) == nrow(bams))
  corr <- seq(nrow(rowData)) %>%
    split(rowData$chr) %>% sapply(
      \(locs) {
        data <- signals[locs, ]
        ffts <- data %>%
          apply(
            2,
            \(v) fft(
              c(v, rep(0, fft_length - length(v)))
            )
          )
        corr_complex <- fft(
          ffts[, 1] * Conj(ffts[, 2]),
          inverse = TRUE
        ) / fft_length
        with_options(list(warn = -1), as.numeric(corr_complex))
      }
    ) %>%
    rowSums
  rbind(
    data.frame(
      shift_forward = seq(-n_corr, -1),
      pearson = corr[(length(corr) - n_corr + 1):length(corr)]
    ),
    data.frame(
      shift_forward = 0:n_corr,
      pearson = corr[1:(n_corr + 1)]
    )
  )
}

analyze_mapq <- function(bams, resolution = 10, test_mapq = c(0, 10, 20, 30, 40, 42), n_corr = 500) {
  mapq_analysis <- bind_rows(
    setNames(
      sapply(
        test_mapq,
        \(val) cross_correlate_mapq2(bams, min_mapq = val, n_corr = n_corr),
        simplify = F
      ),
      test_mapq
    ),
    .id = "mapq"
  )
  ggplot(
    mapq_analysis,
    aes(shift_forward, pearson, color = mapq, group = mapq)
  ) +
    geom_line() +
    scale_color_viridis_d(end = 0.8) +
    theme_bw() +
    labs(
      x = "Shift forward Watson strand reads to meet Crick strand reads (bp)",
      y = "Pearson R of forward-reverse alignments",
      title = paste0(
        "Repli MAPQ selection using ",
        names(bams)[1],
        "-",
        names(bams)[2],
        " held out reads cross-validation"
      )
    )
}
