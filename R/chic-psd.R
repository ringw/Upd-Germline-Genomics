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

granges_plot_nuc_median_phasing <- function(bed_files, legend.position = waiver()) {
  wnd <- tibble(
    enframe(chr.lengths),
    gr = mapply(
      \(n, len) slidingWindows(GRanges(n, IRanges(1, width=len)), 100000L, 100000L)[[1]],
      name,
      value
    )
  ) %>%
    pull(gr) %>%
    GRangesList %>%
    unlist
  names(bed_files) <- names(bed_files) %>% tolower
  phasing <- sapply(
    bed_files,
    \(filename) tibble(
      granges_spacing(
        wnd,
        import(filename, "bed") %>%
          attributes %>%
          with(
            GRanges(
              seqnames,
              IRanges(
                start(ranges) + round(0.5 * width(ranges)),
                width=1
              )
            )
          )
      ),
      pos.plot = ifelse(
        chr == "2R",
        pos + chr.lengths["2L"],
        ifelse(
          chr == "3R",
          pos + chr.lengths["3L"],
          pos
        )
      ),
      seqname.plot = fct_recode(
        chr,
        `2`="2L", `2`="2R", `3`="3L", `3`="3R"
      ) %>%
        relevel("X")
    ),
    simplify=F
  ) %>%
    bind_rows(.id = "celltype")
  
  pl <- ggplot(
    phasing, aes(pos.plot, phasing, group=interaction(chr, celltype),
      color=celltype, fill=celltype)
  ) + facet_grid(
    . ~ seqname.plot, space="free_x", scales="free_x"
  ) + geom_smooth(
    method="loess", span=0.25
  ) + scale_x_continuous(
    breaks = c(1, seq(1, 30) * 2000000),
    minor_breaks = 1000000 + seq(0, 30) * 2000000,
    labels = NULL,
    expand = c(0,0)
  ) + scale_color_manual(
    values = unlist(chic_line_track_colors)
  ) + scale_fill_manual(
    values = muted(unlist(chic_line_track_colors), l=75, c=60)
  ) + labs(
    x = NULL, y = "Median Nucleosome Phasing (bp)"
  ) + theme(
    legend.position = legend.position
  )
  pl <- pl %>% ggplot_build %>% ggplot_gtable
  pl$widths[c(5, 7, 9, 11, 13)] <- unit(c(12, 22, 25, 5, 5), "null")
  pl
}