# Stack feature TSS (as in bed format with flybase id) from BigWig. The required
# features columns are:
#     flybase, chr (as in load_flybase_bed), strand, start, end.
flybase_big_matrix <- function(
  features, bw_path,
  before = 500, after = 1500
) {
  coverage <- import(bw_path, "bigwig") %>% coverage(weight="score")
  mat_colnames <- seq(before+after) %>%
                      replace(. == before+1, "TSS")
  chr_data <- split(features, features$chr)
  chr_data <- mapply(
    \(chr_bed, chr_coverage, chr_length) {
      data_rows <- expand.grid(
        flybase = chr_bed$flybase,
        relative_pos = seq(-before, after-1)
      ) %>%
        left_join(
          with(
            chr_bed,
            data.frame(
              flybase = flybase,
              start = ifelse(strand == "+", start, end),
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
  colnames(chr_data) <- mat_colnames
  chr_data[features$flybase, ]
}

display_tss_tile_matrix <- function(
  data, output_path,
  fc_max = 5,
  fc_filter = 10
) {
  dir.create(dirname(dirname(output_path)), showW=F)
  dir.create(dirname(output_path), showW=F)
  which_tss <- which(colnames(data) == 'TSS')
  newData <- data %>%
    subset(rowMaxs(.) < fc_filter) %>%
    resizeImage(1000, ncol(.), "bilinear") %>%
    melt(varnames = c('row', 'bp'))
  p <- (
    ggplot(newData, aes(bp, row, fill=value))
    + geom_raster()
    + scale_fill_viridis_c(
      option="magma", limits = c(0, fc_max), oob = squish,
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