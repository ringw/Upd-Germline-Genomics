psd_queries <- function(queries, granges, wavelengths) {
  stopifnot(length(unique(queries@ranges@width)) == 1)
  stopifnot(length(unique(granges@ranges@width)) == 1)
  query_wavelengths <- (queries@ranges@width[1] / rev(seq(0, queries@ranges@width[1]))) %>%
    round
  psd_wts <- pull(as.data.frame(table(query_wavelengths)), "Freq", "query_wavelengths")
  query_wavelengths <- query_wavelengths %>% unique
  query_wavelengths <- c(
    query_wavelengths[1:findInterval(min(wavelengths), query_wavelengths)],
    wavelengths,
    query_wavelengths[-(1:findInterval(max(wavelengths), query_wavelengths))]
  )
  square <- \(v) v^2

  matches <- as(findOverlaps(queries, granges), "List")
  input_x <- do.call(
    c,
    mapply(
      \(inds, range_start) granges@ranges@start[inds] - range_start,
      matches,
      queries@ranges@start,
      SIMPLIFY=FALSE
    )
  )
  input_wts <- sin(pi * (input_x - 0.5) / queries@ranges@width[1])^2 / queries@ranges@width[1]
  input_grp <- factor(rep(seq_along(matches), sapply(matches, length)), seq_along(matches))
  specs <- sapply(
    query_wavelengths,
    \(w) (
      input_wts / sqrt(2 * pi) * exp(
        complex(
          im = input_x * 2 * pi / w
        )
      )
    ) %>%
      split(input_grp) %>%
      sapply(sum)
  ) %>%
    matrix(
      nrow = length(query_wavelengths),
      ncol = length(queries),
      # We are taking column data (short and wide) and putting it into rows.
      byrow = TRUE,
      dimnames = list(as.character(query_wavelengths), NULL)
    ) %>%
    Mod %>%
    square
  # Fourier integrals
  fourier_psd_mat <- apply(
    specs,
    2,
    \(v) approx(
      query_wavelengths,
      v,
      xout = as.numeric(names(psd_wts))
    )$y
  )
  ints <- crossprod(fourier_psd_mat, psd_wts)[, 1]
  specs * rep(1/ints, each = nrow(specs))
}

granges_spacing <- function(queries, granges) {
  stopifnot(length(unique(granges@ranges@width)) == 1)
  matches <- as(findOverlaps(queries, granges), "List")
  stat <- sapply(
    matches,
    \(v) median(diff(start(granges)[v]))
  )
  tibble(
    chr = as.factor(seqnames(queries)),
    pos = start(queries) - 1 + round(0.5 * width(queries)),
    phasing = stat
  )
}
