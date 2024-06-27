track_to_heatmap <- function(
  track,
  assay.data.sc,
  mask_track = NULL,
  mask_threshold = 0.1,
  genomic_feature = "TSS",
  before = 500,
  after = 1500
) {
  if (!is.null(mask_track))
    elementMetadata(track)[, 1] <- elementMetadata(track)[, 1] %>%
      replace(which(!sapply(elementMetadata(mask_track)[, 1] < mask_threshold, isFALSE)), NA)
  assay.data.sc.orig <- assay.data.sc
  assay.data.sc$chr <- assay.data.sc$chr %>% factor(names(masked.lengths)) %>% droplevels
  assay.data.sc <- assay.data.sc %>% subset(!is.na(chr))
  track <- track[seqnames(track) %in% levels(assay.data.sc$chr)]
  seqlevels(track) <- levels(assay.data.sc$chr)

  assay.data.sc <- assay.data.sc %>%
    rowwise %>%
    mutate(
      origin = list(
        TSS = c(`+`=start, `-`=end), TES = c(`+`=end, `-`=start)
      )[[genomic_feature]][strand],
      sgn = c(`+`=1, `-`=-1)[strand]
    )

  heatmap_matrix <- mapply(
    \(chr, granges, assay.data.sc) {
      bptrack <- approx(
        granges@ranges@start + 1/2*granges@ranges@width,
        elementMetadata(granges)[,1],
        xout = seq(seqlengths(granges)[chr]),
        na.rm=F
      )$y
      ones <- matrix(1, nrow = nrow(assay.data.sc), ncol = before + after)
      pos <- as.matrix(
        Diagonal(
          x = assay.data.sc$origin
        ) %*% ones
        + Diagonal(
          x = assay.data.sc$sgn
        ) %*% ones %*% Diagonal(x = seq(-before, after - 1))
      )

      m <- bptrack[pos] %>%
        matrix(
          nrow = nrow(assay.data.sc),
          dimnames = list(
            assay.data.sc[, 1, drop=T],
            as.character(seq(-before, after - 1)) %>%
              replace(before + 1, genomic_feature)
          )
        )
    },
    levels(assay.data.sc$chr),
    split(track, seqnames(track)),
    split(assay.data.sc, assay.data.sc$chr),
    SIMPLIFY=FALSE,
    USE.NAMES=F
  ) %>%
    do.call(rbind, .)
  # This indexing replaces out-of-bounds subscripts with NA rows in the result.
  heatmap_matrix <- heatmap_matrix[
    match(assay.data.sc.orig[, 1, drop=T], rownames(heatmap_matrix)),
  ] %>%
    matrix(
      nrow = nrow(assay.data.sc.orig),
      dimnames = list(
        assay.data.sc.orig[, 1, drop=T],
        colnames(.)
      )
    )
}

inter_track_to_heatmap <- function(
  track,
  assay.data.sc,
  mask_track = NULL,
  mask_threshold = 0.1,
  bp_before,
  bp_after,
  n_sample = 99
) {
  if (!is.null(mask_track))
    elementMetadata(track)[, 1] <- elementMetadata(track)[, 1] %>%
      replace(which(!sapply(elementMetadata(mask_track)[, 1] < mask_threshold, isFALSE)), NA)
  assay.data.sc.orig <- assay.data.sc
  assay.data.sc$chr <- assay.data.sc$chr %>% factor(names(masked.lengths)) %>% droplevels
  assay.data.sc <- assay.data.sc %>% subset(!is.na(chr))
  track <- track[seqnames(track) %in% levels(assay.data.sc$chr)]
  seqlevels(track) <- levels(assay.data.sc$chr)

  assay.data.sc <- assay.data.sc %>%
    rowwise %>%
    mutate(
      origin = c(`+`=start + bp_before, `-`=end - bp_after)[strand],
      termin = c(`+`=end + bp_after, `-`=start - bp_before)[strand],
      sgn = c(`+`=1, `-`=-1)[strand]
    )

  heatmap_matrix <- mapply(
    \(chr, granges, assay.data.sc) {
      bptrack <- approx(
        granges@ranges@start + 1/2*granges@ranges@width,
        elementMetadata(granges)[,1],
        xout = seq(seqlengths(granges)[chr]),
        na.rm=F
      )$y
      pos <- t(
        apply(
          assay.data.sc[c("origin", "termin")],
          1,
          \(v) with(
            as.list(v),
            seq(origin, termin, length.out=n_sample+2)[
              -c(1, n_sample+2)
            ]
          )
        )
      )

      m <- bptrack[as.numeric(pos)] %>%
        matrix(
          nrow = nrow(assay.data.sc),
          dimnames = list(
            assay.data.sc[, 1, drop=T],
            str_glue("{seq(1, 99, length.out=n_sample)}%")
          )
        )
    },
    levels(assay.data.sc$chr),
    split(track, seqnames(track)),
    split(assay.data.sc, assay.data.sc$chr),
    SIMPLIFY=FALSE,
    USE.NAMES=F
  ) %>%
    do.call(rbind, .)
  # This indexing replaces out-of-bounds subscripts with NA rows in the result.
  heatmap_matrix <- heatmap_matrix[
    match(assay.data.sc.orig[, 1, drop=T], rownames(heatmap_matrix)),
  ] %>%
    matrix(
      nrow = nrow(assay.data.sc.orig),
      dimnames = list(
        assay.data.sc.orig[, 1, drop=T],
        colnames(.)
      )
    )
}

track_to_heatmap_with_panels <- function(
  track,
  assay.data.sc,
  mask_track = NULL,
  mask_threshold = 0.1,
  bp_flanking = 500,
  n_sample = 99,
  inter_panel_width_kb = 2
) {
  tss_data <- track_to_heatmap(
    track,
    assay.data.sc,
    mask_track=mask_track,
    mask_threshold=mask_threshold,
    genomic_feature="TSS",
    before=bp_flanking,
    after=bp_flanking+1
  )
  inter_data <- inter_track_to_heatmap(
    track,
    assay.data.sc,
    mask_track=mask_track,
    mask_threshold=mask_threshold,
    bp_before=bp_flanking,
    bp_after=bp_flanking,
    n_sample=n_sample
  )
  tes_data <- track_to_heatmap(
    track,
    assay.data.sc,
    mask_track=mask_track,
    mask_threshold=mask_threshold,
    genomic_feature="TES",
    before=bp_flanking,
    after=bp_flanking+1
  )
  colnames(tss_data) <- as.character(seq(-bp_flanking, bp_flanking)) %>%
    replace(bp_flanking + 1, "TSS")
  colnames(tes_data) <- as.character(seq(-bp_flanking, bp_flanking)) %>%
    replace(bp_flanking + 1, "TES")
  data <- cbind(tss_data, inter_data, tes_data)
  attr(data, "x") <- c(
    seq(-bp_flanking, bp_flanking),
    seq(bp_flanking, bp_flanking + 1000 * inter_panel_width_kb, length.out = n_sample + 2)[
      -c(1, n_sample + 2)
    ],
    bp_flanking + 1000 * inter_panel_width_kb + seq(
      0, 2 * bp_flanking
    )
  )
  data <- data
}

smooth_image_columns <- function(
  image, image_filter
) {
  orig_height <- nrow(image)
  image <- rbind(image, matrix(0, nrow=nrow(image), ncol=ncol(image)))
  origin <- ceiling(length(image_filter) / 2)
  image_filter <- c(
    image_filter[origin:length(image_filter)],
    rep(0, nrow(image) - length(image_filter)),
    if (origin != 1)
      image_filter[1:(origin-1)]
    else NULL
  )
  image %>%
    fft %>%
    `*`(fft(image_filter)) %>%
    fft(inverse=T) %>%
    Re %>%
    `/`(nrow(image) * ncol(image)) %>%
    head(orig_height)
}

create_direction_invert_tss_tile_matrix_gradient <- function(...) scale_fill_gradientn(
  colors = c("white", viridis(7, option = "magma", end = 0.9, direction = -1)),
  values = c(0, 1, 1.5, 2, 2.5, 3, 3.5, 4) %>% rescale,
  ...
)

display_tss_tile_matrix <- function(
  data, output_path,
  scale_image_filter = 1,
  fc_max = 5,
  fc_filter = 10,
  direction = 1
) {
  dir.create(dirname(dirname(output_path)), showW=F)
  dir.create(dirname(output_path), showW=F)
  which_tss <- which(colnames(data) == 'TSS')
  newData <- data %>%
    subset(rowMaxs(.) < fc_filter) %>%
    # Clip some extrema
    replace(. > fc_max * 4, fc_max * 4) %>%
    smooth_image_columns(scale_image_filter) %>%
    replace(. > fc_max, fc_max) %>%
    resizeImage(1000, ncol(.), "nearest") %>%
    melt(varnames = c('row', 'bp'))
  p <- (
    ggplot(newData, aes(bp, row, fill=value))
    + geom_raster()
    + scale_fill_viridis_c(
      option="magma", limits = c(0, fc_max), oob = squish,
      direction = direction,
      begin = ifelse(direction == -1, 0.25, 0),
      guide = guide_colorbar(title = "mark/input", barheight = 10)
    )
    + scale_x_continuous(
      expand = c(0, 0),
      breaks = c(1, which_tss, ncol(data)),
      labels = c(
        -(which_tss - 1),
        "TSS",
        ncol(data) - which_tss + 1
      )
    )
    + scale_y_reverse(expand = c(0, 0), breaks = NULL)
    + labs(x = 'bp from TSS', y = paste0('n = ', sum(rowMaxs(data) < fc_filter)))
    + theme(
      plot.margin = unit(c(1, 12, 1, 1), 'pt')
    )
  )
  if (direction == -1)
    p <- p + create_direction_invert_tss_tile_matrix_gradient()
  ggsave(
    output_path,
    p + theme(legend.position = "none"),
    width = 1.5,
    height = 5
  )
  ggsave(
    paste0(dirname(output_path), "/", "Heatmap-Legend.pdf"),
    get_legend(p),
    width = 1,
    height = 3
  )
  ggsave(
    paste0(dirname(output_path), "/", "Heatmap-Legend.png"),
    get_legend(p),
    width = 1,
    height = 3
  )
}