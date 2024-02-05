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
chic_average_profile_limits <- c(-1.25, 0.75)
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
  ) + scale_y_continuous(
    limits = chic_average_profile_limits,
    expand = c(0, 0)
  ) + labs(
    x = "bp (from TSS)", y = "mean(log2(mark/input))"
  )

  ggsave(
    output_path,
    facet_plot,
    width=9,
    height=4,
    dpi=300
  )
  facet_plot

  output_path
}

enrichr_set_names <- function(enrichr, treatment_name, control_name) {
  enrichr@names <- c(treatment_name, control_name)
  enrichr
}

enrichr_grid_ratio <- function(enrichrs, input_count, mod_count) {
  data = data.frame(
    input = sapply(
      enrichrs, \(obj) obj@names[2]
    ),
    mod = sapply(
      enrichrs, \(obj) obj@names[1]
    ),
    multiplier = sapply(
      enrichrs,
      \(obj) obj@theta[1] / (1 - obj@theta[1])
    )
  ) %>%
    mutate(
      input_filename = input,
      mod_filename = mod,
      input_size = input_count[input],
      mod_size = mod_count[mod],
      log_adjust = (log(multiplier) + log(input_size) - log(mod_size)) / log(2),
      # Can test "adjust" using acast and svd: Is the multiplier applied to each
      # sample linearly dependent? Then from the mean of RPKM-adjusted tracks,
      # we can determine what multiplier to apply to the "mark / input" track.
      adjust = exp(log_adjust * log(2))
    )
  data$input <- data$input %>% factor(input_count %>% sort(dec=T) %>% names) %>%
    fct_relabel(\(n) round(input_count[n]/1000/1000, 2) %>% paste0(" Mb") %>% make.unique)
  data$mod <- data$mod %>% factor(mod_count %>% sort(dec=T) %>% names) %>%
    fct_relabel(\(n) round(mod_count[n]/1000/1000, 2) %>% paste0(" Mb") %>% make.unique)
  ggplot(data, aes(input, mod, fill=log_adjust)) + geom_tile(
  ) + scale_x_discrete(
    position = "top"
  ) + scale_y_discrete(limits=rev) + coord_cartesian(
    expand = FALSE
  ) + labs(
    x = paste0("Input Sample (n = ", length(input_count), ")"),
    y = paste0("Mark Sample (n = ", length(mod_count), ")"),
    fill = bquote(log[2]*"(mult)")
  )
}

enrichr_ratio_per_reads_mapped_adjustment <- function(plot.chic.multiplier) {
  data <- plot.chic.multiplier$data %>% acast(
    input_filename ~ mod_filename, value.var = "multiplier"
  )
  # We are normalizing each track by number of reads mapped ("size"). "Mod" is
  # the numerator, so it gets the reciprocal applied. "Input" is in the
  # denominator, so the reciprocal cancels out. Then, we can hit the matrix of
  # normr normalization values (multiplier for treatment / control) with the
  # vectors of FPKM normalization values. This is because the actual
  # contribution of ChIC samples to the track is going to be uniform after the
  # FPKM normalization. Finally, samples contribute to the numerator or
  # denominator proportional to the number of such samples (arithmetic mean).
  my_input_weights <- t(
    plot.chic.multiplier$data$input_size
  )[, rownames(data)]
  my_mod_weights <- 1 / (
    plot.chic.multiplier$data$mod_size
  )[colnames(data)]
  my_input_weights %*% data %*% my_mod_weights / nrow(data) / ncol(data) %>% as.numeric()
}