approx_track <- function(score_track, fine_track) {
  stopifnot(all.equal(as.character(unique(seqnames(fine_track))), intersect(seqlevels(fine_track), unique(seqnames(fine_track)))))
  score_track <- score_track[score_track@seqnames %in% levels(fine_track@seqnames)]
  seqlevels(score_track) <- levels(seqnames(fine_track))

  # The two tracks can now be split according to the same levels (names).
  score_track <- split(score_track, seqnames(score_track))
  # First, initialize the score track from the first entry in each coarse
  # (score_track) seq's first value. Note that some of the score_track seqs have
  # more than one range in them, while the short seqs have a single window
  # covering the entire seq and benefit from this optimization.
  fine_track$score <- rep(
    sapply(score_track, \(gr) gr$score[1]) %>%
      replace(. == 0, NA),
    unlist(table(seqnames(fine_track)))
  )
  fine_track <- split(fine_track, seqnames(fine_track))
  # Now loop over the seqs that are actually affected by the finer windows that
  # are found in the fine_track.
  fine_seqs <- names(score_track)[sapply(score_track, length) > 1]
  for (n in fine_seqs) {
    x <- score_track[[n]]@ranges@start + round(score_track[[n]]@ranges@width / 2)
    y <- score_track[[n]]$score
    # rtracklayer forces us to replace NA with 0.
    x <- x[y != 0]
    y <- y[y != 0]
    if (length(y) == 0)
      fine_track[[n]]$score <- NA
    else if (length(y) == 1)
      fine_track[[n]]$score <- y
    else
      fine_track[[n]]$score <- approx(
        x,
        y,
        xout = fine_track[[n]]@ranges@start + round(fine_track[[n]]@ranges@width / 2),
        rule = 2
      )$y
  }
  unlist(fine_track %>% setNames(NULL))
}

rank_track <- function(track) {
  rnk <- track$score
  rnk[!is.na(rnk)] <- rank(rnk[!is.na(rnk)], ties="random")
  track$rank <- as.integer(rnk)
  track
}

bin_track_by_rank <- function(track, nbins = 1000) {
  # newnames <- with(
  #   attributes(track[!is.na(track$rank)]),
  #   str_glue("{seqnames}_{ranges@start}")
  # )
  newperm <- order(track$rank[!is.na(track$rank)])
  newlookup <- order(track$rank) %>% head(sum(!is.na(track$rank)))
  newtrack <- track@ranges@width[newlookup] %>%
    IRanges(
      start = cumsum(c(1, .[-length(.)])),
      width = . # ,
      # names = newnames[newperm]
    )
  bin_loc <- seq(1, sum(newtrack@width), length.out=1+nbins)[-(1+nbins)]
  bin_fct <- rep(
    seq(nbins),
    diff(c(0, findInterval(bin_loc[-1], newtrack@start), length(newtrack)))
  )
  bin <- rep(as.integer(NA), length(newtrack))
  bin[newlookup] <- bin_fct
  bin_sizes_iranges <- newtrack@width %>% split(bin_fct) %>% sapply(sum)
  size_of_bin_bp <- rep(as.integer(NA), length(newtrack))
  size_of_bin_bp[newlookup] <- bin_sizes_iranges[bin_fct]
  GRanges(track, bin = bin, size_of_bin_bp = size_of_bin_bp)
}

granges_bin_to_projection <- function(granges) {
  proj_nrow <- max(granges$bin, na.rm=T)
  granges %>%
    split(seqnames(.)) %>%
    sapply(
      \(gr) sparseMatrix(
        i = gr$bin[!is.na(gr$bin)],
        j = which(!is.na(gr$bin)),
        x = gr@ranges@width / gr$size_of_bin_bp,
        dims = c(proj_nrow, length(gr))
      ),
      simplify=F
    )
}

apply_granges_projection <- function(projs, granges) {
  rowSums(
    mapply(
      \(mat, v) (mat %*% v)[,1],
      projs,
      split(granges$score, factor(seqnames(granges), names(projs)))
    )
  )
}

plot_repli_track_raster <- function(data) {
  sample_size_mb <- round(sum(data[, "sample_size_bp"], na.rm=T) / 1000 / 1000, 1)
  data[, grep("^H", colnames(data))] <- log(data[, grep("^H", colnames(data))]) / log(2)
  data <- melt(data) %>%
    mutate(series = series %>% fct_recode(`Timing Est.`="repli")) %>%
    group_by(timing) %>%
    mutate(
      series = series %>%
        factor(c(levels(.), "Quartile")) %>%
        fct_relevel("Timing Est.", "Quartile"),
      x = value[which.max(series == "ranking")]
    )
  x_values <- unique(data$value[data$series == "ranking"])
  (
    ggplot(data, aes(x, series, fill=value))
    + geom_raster(aes(), subset(data, series == "Timing Est."))
    + scale_fill_viridis_c(
      option = "turbo", begin = 0.55, end = 1, direction = -1,
      guide = guide_colorbar(title = "timing"),
      limits = c(0.125, 0.875)
    )
    + new_scale_fill()
    + geom_raster(
      aes(fill=q),
      tribble(
        ~x, ~series, ~q,
        quantile(x_values, 0.125), "Quartile", "Q1",
        quantile(x_values, 0.375), "Quartile", "Q2",
        quantile(x_values, 0.625), "Quartile", "Q3",
        quantile(x_values, 0.875), "Quartile", "Q4"
      ) %>%
        mutate(value = 0)
    )
    + annotate(
      "text",
      x = quantile(x_values, c(0.125, 0.375, 0.625, 0.875)),
      y = 4,
      label = c("Q1", "Q2", "Q3", "Q4"),
      size = 3,
      color = "white",
      family = "Helvetica"
    )
    + scale_fill_manual(
      values = turbo(101)[c(87, 82, 78, 72)],
      guide = guide_none()
    )
    + new_scale_fill()
    + geom_raster(aes(fill=value2), mutate(subset(data, grepl("^H", series)), value2=value))
    # + scale_fill_viridis_c(
    #   option = "magma", direction = -1,
    #   guide = guide_colorbar(title = "mark/input")
    # )
    + create_direction_invert_tss_tile_matrix_gradient(limits = c(-0.5, 1.6), oob = squish)
    + scale_x_continuous(labels=percent, name = str_glue("Sorted Timing Est ({sample_size_mb} Mb)"))
    + scale_y_discrete(limits=rev)
    + coord_cartesian(expand=F)
    + guides(fill = guide_colorbar(title = "log2(mark/input)"))
    + theme(aspect.ratio = 1)
  )
}