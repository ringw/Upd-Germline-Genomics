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
chic_average_profile_limits <- c(0.37, 3.3)
chic_average_profiles <- function(
  chic_factor,
  chic_path,
  metafeatures_path,
  chic_driver,
  legend_title,
  quartile_colors
) {
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
              paste0(chic_path, "/", chic_driver, "_", mark, ".FE.bw"),
              before = before,
              after = after
            ) %>%
              subset(rowAlls(. > 0.001)) %>%
              # log %>%
              # `/`(log(2)) %>%
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
    ) %>%
    ungroup %>%
    mutate(marks = marks %>% factor(chic.mark.data$mark, ordered=TRUE))
  facet_data %>% ggplot(
    aes(x=pos, y=l2FC, color=genes, linewidth=genes, group=genes)
  ) + geom_line(
    data = cross_join(
      tribble(~pos, ~l2FC, -Inf, 1, Inf, 1),
      tibble(
        genes = NA,
        marks = factor(chic.mark.data$mark, chic.mark.data$mark, ordered=TRUE)
      )
    ),
    color = "darkred",
    linewidth = 0.25
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
  ) + coord_cartesian(expand=FALSE) + labs(
    x = "bp (from TSS)", y = "mean(mark/input)"
  )
}

chic_average_gene_list_profiles <- function(
  gene_list,
  chic_path,
  metafeatures_path,
  chic_driver
) {
  metafeatures <- read.csv(metafeatures_path, row.names = 1)
  metafeatures <- metafeatures[gene_list, ] %>%
    subset(chr %in% names(chr.lengths))

  before <- 500
  after <- 1500
  mark_tracks <- sapply(
    chic.mark.data$mark,
    \(mark) flybase_big_matrix(
      metafeatures %>% subset(!is.na(chr) & !is.na(start) & !is.na(end)) %>% subset(!duplicated(flybase)),
      paste0(chic_path, "/", chic_driver, "_", mark, ".FE.bw"),
      before = before,
      after = after
    ) %>%
      subset(rowAlls(. > 0.001)) %>%
      # log %>%
      # `/`(log(2)) %>%
      colMeans,
    simplify=FALSE
  )
  facet_data <- mark_tracks %>%
    sapply(
      \(v) data.frame(pos=seq(-before, after-1), l2FC=v),
      simplify=F) %>%
    bind_rows(.id = "marks") %>%
    mutate(marks = marks %>% factor(chic.mark.data$mark))
  facet_data %>% ggplot(
    aes(x=pos, y=l2FC)
  ) + geom_line(
    data = tribble(~pos, ~l2FC, -Inf, 1, Inf, 1),
    color = "darkred",
    linewidth = 0.25
  ) + geom_line(color = "goldenrod", linewidth = 1) + facet_wrap(
    vars(marks)
  ) + scale_y_continuous(
    limits = chic_average_profile_limits,
    expand = c(0, 0)
  ) + coord_cartesian(expand = FALSE) + labs(
    x = "bp (from TSS)", y = "mean(mark/input)"
  )
}

chic_average_gene_list_gene_body <- function(
  gene_list,
  chic_bw,
  metafeatures_path,
  after_tss = 1000,
  before_tes = 1000,
  num_points_to_sample = 100
) {
  metafeatures <- read.csv(metafeatures_path, row.names = 1)
  metafeatures <- metafeatures[gene_list, ] %>%
    subset(chr %in% names(chr.lengths))
  coverage <- import(chic_bw, "bigwig") %>% coverage(weight="score")

  metafeatures <- metafeatures %>%
    rownames_to_column %>%
    mutate(
      display_start = ifelse(
        strand == "+",
        start + after_tss,
        end - after_tss
      ),
      display_end = ifelse(
        strand == "+",
        end - before_tes,
        start + before_tes
      )
    )
  split_features <- metafeatures %>% split(.$transcript.length > (after_tss + before_tes))

  # Genes to downsample a variable-sized window from the chic track.
  display_genes <- split_features$`TRUE` %>%
    rowwise %>%
    mutate(
      values = seq(
        display_start, display_end, length.out = num_points_to_sample
      ) %>%
        round %>%
        list
    )
  display_gene_tracks <- display_genes %>%
    split(.$chr) %>%
    mapply(
      \(df, track) df$values %>%
        do.call(c, .) %>%
        `[`(track, .) %>%
        matrix(nrow = nrow(df), byrow = TRUE, dimnames = list(df$rowname, NULL)),
      .,
      coverage[names(.)],
      SIMPLIFY=FALSE
    ) %>%
    do.call(rbind, .)

  average_gene_tracks <- split_features$`FALSE` %>%
    split(.$chr) %>%
    mapply(
      \(df, track) c(df$display_start, df$display_end) %>%
        `[`(track, .) %>%
        matrix(nrow = nrow(df), ncol = 2) %>%
        rowMeans %>%
        rep(times = num_points_to_sample) %>%
        matrix(nrow = nrow(df), ncol = num_points_to_sample, dimnames = list(df$rowname, NULL)),
      .,
      coverage[names(.)],
      SIMPLIFY=FALSE
    ) %>%
    do.call(rbind, .)

  rbind(
    display_gene_tracks,
    average_gene_tracks
  )[metafeatures$rowname, ]
}

pull_chic_average_gene_list_paneled_data <- function(
  gene_list,
  chic_bw,
  metafeatures_path,
  num_tss,
  num_tes
) {
  metafeatures <- read.csv(metafeatures_path, row.names = 1)
  metafeatures <- metafeatures[gene_list, ] %>%
    subset(chr %in% names(chr.lengths))

  tss_data <- flybase_big_matrix(
    metafeatures,
    chic_bw,
    "TSS",
    before = num_tss,
    after = num_tss
  )
  tes_data <- flybase_big_matrix(
    metafeatures,
    chic_bw,
    "TES",
    before = num_tes,
    after = num_tes
  )
  inter_data = chic_average_gene_list_gene_body(
    gene_list,
    chic_bw,
    metafeatures_path,
    num_points_to_sample = 99
  )
  result <- tibble(
    labels = c(
      (-num_tss):(-1),
      "TSS",
      paste0("+", 1:(num_tss - 1)),
      paste0(1:99, "%"),
      (-num_tes):(-1),
      "TES",
      paste0("+", 1:(num_tes - 1))
    ),
    track_value = c(
        colMeans(tss_data),
        colMeans(inter_data),
      colMeans(tes_data)
    )
  )
  colnames(result)[2] <- basename(chic_bw) %>% str_replace("[.].*", "")
  result
}

pull_chic_average_gene_list_paneled_profiles_data <- function(
  gene_lists,
  chic_path,
  metafeatures_path,
  chic_driver,
  tss_size,
  inter_size,
  tes_size,
  num_tss = 500,
  num_tes = 500
) {
  num_inter <- 99
  panel_data <- tibble(
    chic_driver = chic_driver,
    chic_mark = chic.mark.data$mark,
    filename = paste0(chic_path, "/", chic_driver, "_", chic_mark, ".FE.bw")
  ) %>%
    cross_join(
      tibble(
        group = names(gene_lists) %>%
          coalesce(as.character(seq_along(gene_lists))) %>%
          factor(., ., ordered=TRUE),
        gene_list = gene_lists
      )
    ) %>%
    rowwise %>%
    mutate(
      track_data = list(
        pull_chic_average_gene_list_paneled_data(
          gene_list,
          filename,
          metafeatures_path,
          num_tss=num_tss,
          num_tes=num_tes
        )
      )
    )
  panel_data %>%
    rowwise %>%
    reframe(
      chic_mark,
      group,
      gene_list = list(gene_list),
      x = c(
        seq(0, tss_size, length.out = num_tss * 2),
        seq(tss_size, tss_size + inter_size, length.out = num_inter),
        seq(tss_size + inter_size, tss_size + inter_size + tes_size, length.out = num_tes * 2)
      ),
      x_label = track_data$labels,
      y = track_data %>% pull(2)
    ) %>%
    mutate(chic_mark = chic_mark %>% factor(chic.mark.data$mark, ordered=TRUE))
}

chic_quartile_gene_list_paneled_profiles <- function(
  quartile.factor,
  chic_path,
  metafeatures_path,
  chic_driver,
  legend_title,
  quartile_colors,
  heatmap_before_after = 500,
  # 1 kb - before and after
  tss_size = 1,
  # 2 kb long - by convention
  inter_size = 2,
  tes_size = 1
) {
  gene_lists <- names(quartile.factor) %>% split(quartile.factor)
  data <- pull_chic_average_gene_list_paneled_profiles_data(
    gene_lists,
    chic_path,
    metafeatures_path,
    chic_driver,
    tss_size,
    inter_size,
    tes_size,
    num_tss=heatmap_before_after,
    num_tes=heatmap_before_after
  ) %>%
    subset(select = -gene_list)

  label_data <- data %>%
    subset(
      as.numeric(chic_mark) == 1 & as.numeric(group) == 1,
      select = c(x, x_label)
    )
  label_data_show <- label_data %>%
    subset(
      x_label %in% c(
        paste0("-", heatmap_before_after),
        paste0("+", heatmap_before_after - 1),
        "50%",
        "TSS",
        "TES"
      )
    )
  minor_breaks_show <- label_data %>% subset(x_label %in% c("25%", "75%"))

  custom_plot <- data %>% ggplot(
    aes(x, y, color=group, linewidth=group, group=group)
  ) + geom_tile(
    # Colored tile (TSS and TES)
    data = tribble(
      ~x, ~y, ~width,
      # Full-size background color
      1,
      1,
      Inf,
      tss_size / 2,
      1,
      tss_size,
      tss_size + inter_size + (tes_size / 2),
      1,
      tes_size
    ) %>%
      cross_join(
        tibble(
          height = Inf
        )
      ),
    aes(width=width, height=height, color=NULL, linewidth=NULL, group=NULL),
    color = "transparent",
    fill = rep(
      # Use cream color from magma scale
      c("#eeeeee", viridis(2, option = "magma", end = 0.99)[2])[
        c(1, 2, 2)
      ],
      length(chic.mark.data$mark)
    )
  ) + facet_wrap(
    vars(chic_mark)
  ) + scale_color_manual(
    values = quartile_colors,
    guide = guide_legend(title = legend_title, override.aes = list(fill = "transparent"))
  ) + scale_linewidth_manual(
    values = c(0.5, 0.75, 1.5, 1.75),
    guide = guide_legend(title = legend_title)
  ) + coord_cartesian(
    c(0, tss_size + inter_size + tes_size),
    chic_average_profile_limits,
    expand=FALSE
  ) + theme(
    panel.background = element_rect(fill = NA),
    panel.ontop = TRUE
  )
  custom_plot + geom_line(
    # Dark red line at the origin (enrichment of 1)
    data = cross_join(
      tribble(~x, ~y, -Inf, 1, Inf, 1),
      tibble(
        group = NA,
        chic_mark = factor(chic.mark.data$mark, chic.mark.data$mark, ordered=TRUE)
      )
    ),
    color = "darkred",
    linewidth = 0.25
  ) + geom_line() + labs(
    x = "base pairs", y = "mean(mark/input)"
  ) + scale_x_continuous(
    breaks = label_data_show$x,
    labels = label_data_show$x_label,
    minor_breaks = minor_breaks_show$x
  ) + theme(
    panel.margin = unit(25, "pt")
  )
}

chic_custom_gene_list_paneled_profile <- function(
  gene_list,
  chic_path,
  metafeatures_path,
  chic_driver,
  heatmap_before_after = 500,
  # 1 kb - before and after
  tss_size = 1,
  # 2 kb long - by convention
  inter_size = 2,
  tes_size = 1
) {
  data <- pull_chic_average_gene_list_paneled_profiles_data(
    list(gene_list),
    chic_path,
    metafeatures_path,
    chic_driver,
    tss_size,
    inter_size,
    tes_size,
    num_tss=heatmap_before_after,
    num_tes=heatmap_before_after
  ) %>%
    subset(select = -gene_list)

  label_data <- data %>%
    subset(
      as.numeric(chic_mark) == 1 & as.numeric(group) == 1,
      select = c(x, x_label)
    )
  label_data_show <- label_data %>%
    subset(
      x_label %in% c(
        paste0("-", heatmap_before_after),
        paste0("+", heatmap_before_after - 1),
        "50%",
        "TSS",
        "TES"
      )
    )
  minor_breaks_show <- label_data %>% subset(x_label %in% c("25%", "75%"))

  custom_plot <- data %>% ggplot(
    aes(x, y)
  ) + geom_tile(
    # Colored tile (TSS and TES)
    data = tribble(
      ~x, ~y, ~width,
      # Full-size background color
      1,
      1,
      Inf,
      tss_size / 2,
      1,
      tss_size,
      tss_size + inter_size + (tes_size / 2),
      1,
      tes_size
    ) %>%
      cross_join(
        tibble(
          height = Inf
        )
      ),
    aes(width=width, height=height),
    color = "transparent",
    fill = rep(
      # Use cream color from magma scale
      c("#eeeeee", viridis(2, option = "magma", end = 0.99)[2])[
        c(1, 2, 2)
      ],
      length(chic.mark.data$mark)
    )
  ) + facet_wrap(
    vars(chic_mark)
  ) + coord_cartesian(
    c(0, tss_size + inter_size + tes_size),
    chic_average_profile_limits,
    expand=FALSE
  ) + theme(
    panel.background = element_rect(fill = NA),
    panel.ontop = TRUE
  )
  custom_plot + geom_line(
    # Dark red line at the origin (enrichment of 1)
    data = cross_join(
      tribble(~x, ~y, -Inf, 1, Inf, 1),
      tibble(
        group = NA,
        chic_mark = factor(chic.mark.data$mark, chic.mark.data$mark, ordered=TRUE)
      )
    ),
    color = "darkred",
    linewidth = 0.25
  ) + geom_line(color = "goldenrod", linewidth = 1) + labs(
    x = "base pairs", y = "mean(mark/input)"
  ) + scale_x_continuous(
    breaks = label_data_show$x,
    labels = label_data_show$x_label,
    minor_breaks = minor_breaks_show$x
  ) + theme(
    panel.margin = unit(25, "pt")
  )
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