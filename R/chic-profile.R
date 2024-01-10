# Sorted BED to quartiles. For deepTools, we had sorted the BED from high-to-low
# and now will construct a feature which is in reverse order.
bed_flybase_quartile_factor <- function(bed_path, metafeatures_path, nquartiles=4) {
  bed <- load_flybase_bed(bed_path)
  bed$rank <- rev(seq(nrow(bed)))
  metafeatures <- read.csv(metafeatures_path, row.names = 1) %>%
    subset(!is.na(flybase) & !duplicated(flybase)) %>%
    rownames_to_column
  features <- inner_join(metafeatures, bed, by = "flybase")
  groups <- cut(
    features$rank,
    round(seq(0, length(features$flybase), length.out=1+nquartiles))
  ) %>%
    setNames(features$rowname)
  levels(groups) <- paste0('Q', seq(nquartiles))
  groups
}

# Analysis with TSS profile as x-axis, with ChIC tracks.
chic_average_profiles <- function(
  chic_factor,
  chic_path,
  metafeatures_path,
  chic_driver,
  legend_title,
  quartile_colors,
  output_path
) {
  dir.create(dirname(dirname(output_path)), showW=F)
  dir.create(dirname(output_path), showW=F)

  chic_table <- table(chic_factor)
  chic_factor <- chic_factor %>%
    fct_relabel(\(n) n %>% sapply(\(n) paste0(n, " (n = ", chic_table[n], ")")))

  metafeatures <- read.csv(metafeatures_path, row.names = 1)
  annotations <- metafeatures %>%
    rownames_to_column %>%
    left_join(
      data.frame(group=chic_factor) %>% rownames_to_column, "rowname"
    )
  facets <- expand.grid(
    genes = levels(chic_factor),
    marks = chic.mark.data$mark
  )
  before <- 500
  after <- 1500
  facets <- facets %>%
    cbind(
      as_tibble(
        list(
          profile = mapply(
            \(genes, mark) flybase_big_matrix(
              annotations %>% subset(group == genes & !is.na(chr) & !is.na(start) & !is.na(end)) %>% subset(!duplicated(flybase)),
              paste0(chic_path, "/", chic_driver, "_", mark, ".q5.bw"),
              before = before,
              after = after
            ) %>%
              subset(rowAlls(. > 0.001)) %>%
              log %>%
              `/`(log(2)) %>%
              colMeans,
            facets$genes,
            facets$marks,
            SIMPLIFY = F
          )
        )
      )
    )
  facets_easy = dcast(facets, genes ~ marks, value.var='profile') %>% column_to_rownames('genes')
  facet_data <- pull(facets, profile) %>%
    sapply(
      \(v) data.frame(pos=seq(-before, after-1), l2FC=v),
      simplify=F) %>%
    setNames(seq_along(.)) %>%
    bind_rows(.id = "profile_index") %>%
    left_join(
      facets %>% subset(select=-profile) %>% rownames_to_column("profile_index"),
      .,
      by = "profile_index"
    )
  facet_plot <- facet_data %>% ggplot(
    aes(x=pos, y=l2FC, color=genes, linewidth=genes, group=genes)
  ) + geom_line() + facet_wrap(
    vars(marks)
  ) + scale_color_manual(
    values = quartile_colors,
    guide = guide_legend(title = legend_title)
  ) + scale_linewidth_manual(
    values = c(0.5, 0.75, 1.5, 1.75),
    guide = guide_legend(title = legend_title)
  ) + labs(
    x = "bp (from TSS)", y = "mean(log2(mark/input))"
  )

  ggsave(
    output_path,
    facet_plot,
    width=9,
    height=3,
    dpi=300
  )
  facet_plot

  output_path
}