# Stack feature TSS (as in bed format with flybase id) from BigWig. The required
# features columns are:
#     flybase, chr (as in load_flybase_bed), strand, start, end.
flybase_big_matrix <- function(
  features, bw_path,
  genomic_feature = "TSS",
  before = 500, after = 1500
) {
  coverage <- import(bw_path, "bigwig") %>% coverage(weight="score")
  mat_colnames <- seq(before+after) %>%
                      replace(. == before+1, "TSS")
  chr_data <- split(features, features$chr)
  chr_data <- mapply(
    \(chr_bed, chr_coverage, chr_length) {
      if (type_sum(rownames(chr_bed)) != "chr")
        rownames(chr_bed) <- as.character(seq(nrow(chr_bed)))
      data_rows <- expand.grid(
        flybase = chr_bed$flybase,
        relative_pos = seq(-before, after-1)
      ) %>%
        left_join(
          with(
            chr_bed,
            data.frame(
              flybase = flybase,
              start = chr_bed[
                cbind(
                  rownames(chr_bed),
                  list(
                    TSS = c(`+`="start", `-`="end"),
                    TES = c(`+`="end", `-`="start")
                  )[[genomic_feature]][strand]
                )
              ] %>%
                as.numeric,
              sign = ifelse(strand == "+", 1, -1)
            )
          ),
          "flybase"
        )
      data_rows <- data_rows %>%
        within(
          {
            pos <- start + sign * relative_pos
            is_finite <- between(pos, 1, chr_length)
            pos <- pos %>% replace(!is_finite, 1)
          }
        )
      chr_data <- matrix(
        as.numeric(chr_coverage[data_rows$pos]) %>%
          replace(!data_rows$is_finite, NA),
        nrow = nrow(chr_bed),
        ncol = before + after,
        dimnames = list(chr_bed$flybase, mat_colnames)
      )
      chr_data
    },
    chr_data,
    coverage[names(chr_data)],
    chr.lengths[names(chr_data)],
    SIMPLIFY=F
  )
  chr_data <- chr_data %>% do.call(rbind, .)
  if (is.null(chr_data))
    return(matrix(nrow = 0, ncol = length(mat_colnames), dimnames = list(NULL, mat_colnames)))
  colnames(chr_data) <- mat_colnames
  chr_data[features$flybase,, drop=F]
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