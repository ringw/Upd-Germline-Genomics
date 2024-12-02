cluster_colors = list(
  # green
  germline=hcl(123, 113, 80),
  # magenta
  somatic="#FF5AEE",
  # cyan 
  spermatocyte=hcl(176, 62, 80),
  # old spermatocyte label
  differentiated=hcl(86, 20, 75),
  muscle=hcl(17, 112, 58),
  others=hsv(0,0,0.9),

  # purple
  somaticprecursor=hcl(277, 86, 66)
)

cluster_bar_color1 = list(
  # forest green
  germline=hcl(123, 67, 50),
  # deep magenta
  somatic=hcl(320, 75, 50)
)

cluster_bar_color2 = list(
  # muted green
  germline = hcl(123, 27, 61),
  # rose
  somatic = hcl(0, 47, 62)
)

sc_quartile_colors = c(
  low = hcl(252, 66, 37),
  hcl(200, 37, 60),
  hcl(137, 37, 80),
  high = hcl(77, 87, 87)
)
sc_quartile_colors = viridis(5)[-1]
names(sc_quartile_colors)[c(1,4)] = c('low', 'high')

sc_quartile_annotations = c(
  off = "#A8A7A7",
  hcl(200, 37, 60),
  hcl(137, 37, 80),
  high = hcl(77, 87, 87)
)

Upd_sc_plot_idents <- function(Upd_sc) {
  # Rename spermatocyte to differentiated below.
  cbind(
    Upd_sc@meta.data,
    Upd_sc[['umap']]@cell.embeddings,
    ident=FetchData(Upd_sc, "ident") # %>%
      # recode(spermatocyte='differentiated') %>%
      # fct_relevel(c('germline', 'somatic', 'muscle', 'differentiated'))
  ) %>% ggplot(
    aes(umap_1, umap_2, color=ident)
  ) + rasterize(geom_point(
    shape = 20, size = 0.001
  ), dpi=120, scale=0.5) + scale_color_manual(
    values=unlist(cluster_colors),
    guide=guide_legend(title=NULL, override.aes = list(size=4))
  ) + scale_x_continuous(
    breaks=c(-10,0,10)
  ) + scale_y_continuous(breaks=c(-10,0,5)) + theme_cowplot() + theme(
    aspect.ratio = 0.75
  )
}

Upd_sc_plot_idents_pcasubset <- function(Upd_sc) {
  # Rename spermatocyte to differentiated below.
  FetchData(
    Upd_sc,
    c("ident", "pcasubset_1", "pcasubset_2")
  ) %>%
    ggplot(
      aes(pcasubset_1, pcasubset_2, color=ident)
    ) + rasterize(geom_point(
      stroke=NA, size = 0.25
    ), dpi=80, scale=0.5) + scale_color_manual(
      values=unlist(cluster_colors),
      guide=guide_legend(title=NULL, override.aes = list(size=4))
    ) + coord_cartesian(
      c(-25, 25), c(-20, 20), expand=F
    ) + labs(
      x = "PC 1", y = "PC 2"
    ) + theme_cowplot() + theme(
      aspect.ratio = 0.75
    )
}

Upd_sc_plot_idents_pcasubset_batch <- function(Upd_sc) {
  data <- FetchData(Upd_sc, c("ident", "batch", "pcasubset_1", "pcasubset_2"))
  plot_data <- \(n) data %>%
    subset(batch == n) %>%
    ggplot(
      aes(pcasubset_1, pcasubset_2, color=ident)
    ) + rasterize(geom_point(
      stroke=NA, size = 0.25
    ), dpi=80, scale=0.5) + scale_color_manual(
      values=unlist(cluster_colors)[1:2],
      guide=guide_legend(title=NULL, override.aes = list(size=4))
    ) + coord_cartesian(
      c(-22, 22), c(-20, 20), expand=F
    ) + labs(
      x = "PC 1", y = "PC 2"
    ) + theme_cowplot() + theme(
      aspect.ratio = 2/3
    )
  fig <- gtable(
    # 1 inch margin width as well as flexible width columns.
    widths = unit(1, rep(c("null", "in"), c(2, 1))),
    heights = unit(rep(c(0.5, 0.85), c(1, 2)), rep(c("in", "null"), c(1, 2)))
  ) %>%
    gtable_add_grob(
      list(
        linesGrob(x = 0),
        textGrob("nos genotype"),
        textGrob("tj genotype"),
        as_grob(plot_data("nos.1") + theme(legend.position = "none")),
        as_grob(plot_data("nos.2") + theme(legend.position = "none")),
        as_grob(plot_data("tj.1") + theme(legend.position = "none")),
        as_grob(plot_data("tj.2") + theme(legend.position = "none")),
        as_grob(get_legend(plot_data("nos.1")))
      ),
      t = c(1, 1, 1, 2, 3, 2, 3, 2),
      l = c(2, 1, 2, 1, 1, 2, 2, 3),
      b = c(3, 1, 1, 2, 3, 2, 3, 3)
    )
}

Upd_sc_plot_subset <- function(sc.plot.idents) {
  sc.plot.idents$data <- sc.plot.idents$data %>%
    subset(ident %in% c('germline', 'somatic'))
  sc.plot.idents + scale_y_continuous(
    limits=c(-6,6), breaks=c(-5,0,5)
  ) + theme(
    aspect.ratio = 0.5
  )
}

Upd_sc_feature_plot <- function(Upd_sc, gene, cells, assay = "RNA") {
  DefaultAssay(Upd_sc) <- assay
  dplyr_rename_lookup_gene = setNames("LogNormalize", gene)
  gene.data = FetchData(Upd_sc, c("umap_1", "umap_2", gene, "ident")) %>%
    plyr::rename(all_of(dplyr_rename_lookup_gene))
  gene.data <- gene.data[cells, ]
  gene.max.intensity <- NA
  gene.data %>% ggplot(
    aes(umap_1, umap_2, color=LogNormalize)
  ) + rasterize(geom_point(
    # Scale point because we cut the print size in half for the graphic
    shape = 20, size = 0.5 * 0.5, stroke = NA
  ), dpi=600) + scale_color_viridis_c(
    begin = 0.2,
    limits = c(0, gene.max.intensity), oob = squish
  ) + scale_x_continuous(
    breaks=c(-10,0,10)
  ) + scale_y_continuous(breaks=c(-10,0,5)) + theme_cowplot(
    font_size = 14 * 0.5,
    line_size = 0.5 * 0.5
  ) + labs(
    tag=gene
  ) + theme(
    plot.tag.position = c(0.14, 0.94),
    aspect.ratio = 0.75
  )
}

Upd_sc_group_plot <- function(Upd_sc, group.by, cells) {
  dplyr_rename_lookup_gene = setNames("group", group.by)
  gene.data = FetchData(Upd_sc, c("umap_1", "umap_2", group.by, "ident")) %>%
    plyr::rename(all_of(dplyr_rename_lookup_gene))
  gene.data <- gene.data[cells, ]
  gene.data %>% ggplot(
    aes(umap_1, umap_2, color=group)
  ) + geom_point(
    # Scale point because we cut the print size in half for the graphic
    shape = 20, size = 0.5 * 0.5, stroke = NA
  ) + scale_x_continuous(
    breaks=c(-10,0,10)
  ) + scale_y_continuous(breaks=c(-10,0,5)) + theme_cowplot(
    font_size = 14 * 0.5,
    line_size = 0.5 * 0.5
  ) + theme(
    plot.tag.position = c(0.14, 0.94),
    aspect.ratio = 0.75
  ) + guides(
    color = guide_legend(override.aes = list(size = 2))
  )
}

# Single-cell dim reduc figure generator. Returns list of filenames.
Upd_sc_figures = function(figures_dir, Upd_sc) {
  Upd_sc = Upd_sc %>% NormalizeData
  dir.create(figures_dir, showW = F)
  # Rename spermatocyte to differentiated
  cbind(
    Upd_sc@meta.data,
    Upd_sc[['umap']]@cell.embeddings,
    ident=Idents(Upd_sc) %>%
      recode(spermatocyte='differentiated') %>%
      fct_relevel(c('germline', 'somatic', 'muscle', 'differentiated'))
  ) %>% ggplot(
    aes(umap_1, umap_2, color=ident)
  ) + rasterize(geom_point(
    size = 0.01
  ), dpi=120, scale=0.5) + scale_color_manual(
    values=unlist(cluster_colors),
    guide=guide_legend(title=NULL, override.aes = list(size=4))
  ) + scale_x_continuous(
    breaks=c(-10,0,10)
  ) + scale_y_continuous(breaks=c(-10,0,5)) + theme_cowplot()
  filenames <- NULL
  filenames <- c(paste(figures_dir, 'UMAP.svg', sep='/'), filenames)
  ggsave(filenames[1], width=8, height=4.5)
  filenames <- c(paste(figures_dir, 'UMAP.pdf', sep='/'), filenames)
  ggsave(filenames[1], width=8, height=4.5)
  umap_plot <- last_plot()
  umap_plot$data <- umap_plot$data %>%
    subset(ident %in% c('germline', 'somatic'))
  print(umap_plot + scale_y_continuous(limits=c(-6,NA), breaks=c(-5,0,5)))
  filenames <- c(paste(figures_dir, 'UMAP-Subset.svg', sep='/'), filenames)
  ggsave(filenames[1], width=8, height=3)
  filenames <- c(paste(figures_dir, 'UMAP-Subset.pdf', sep='/'), filenames)
  ggsave(filenames[1], width=8, height=3)

  for (gene in c('AGO3','RpL22-like','vas','nos','tj','lncRNA:roX1','lncRNA:roX2','Mst87F','soti','sunz','w-cup','Act57B','Dl')) {
    gene.save = gene %>% str_replace('lncRNA:', '')
    gene.data = Upd_sc@meta.data %>%
      cbind(
        Upd_sc[['umap']]@cell.embeddings,
        ident=Idents(Upd_sc),
        LogNormalize=Upd_sc[['RNA']][gene,] %>% as.numeric
      )
    gene.max.intensity = quantile(gene.data$LogNormalize, 0.99)
    gene.max.intensity = c(Mst87F=5, Act57B=4, soti=3, sunz=2)[gene] %>% replace(is.na(.), gene.max.intensity)
    gene.data %>% ggplot(
      aes(umap_1, umap_2, color=LogNormalize)
    ) + rasterize(geom_point(
      size = 0.01
    ), dpi=240) + # + scale_color_viridis_c(
      # option='magma', limits=c(0,gene.max.intensity), oob=squish
    scale_color_viridis_c(
      begin = 0.2,
      limits = c(0, gene.max.intensity), oob = squish
    ) + scale_x_continuous(
      breaks=c(-10,0,10)
    ) + scale_y_continuous(breaks=c(-10,0,5)) + theme_cowplot(
    ) + labs(
      tag=gene.save
    ) + theme(
      plot.tag.position = c(0.1, 0.94)
    )
    filenames <- c(paste0(figures_dir, '/UMAP-', gene.save, '.svg'), filenames)
    ggsave(filenames[1], width=8, height=4.5)
    filenames <- c(paste0(figures_dir, '/UMAP-', gene.save, '.pdf'), filenames)
    ggsave(filenames[1], width=8, height=4.5)
  }
  filenames
}

load_cell_cycle_score_drosophila <- function(cell_cycle_drosophila_path, metafeatures_path) {
  cell_cycle_drosophila <- cell_cycle_drosophila_path %>% read.csv
  metafeatures <- metafeatures_path %>% read.csv(row.names=1)
  colnames(cell_cycle_drosophila) <- colnames(cell_cycle_drosophila) %>%
    replace(. == "geneID", "flybase")
  cell_cycle_drosophila <- cell_cycle_drosophila %>%
    inner_join(metafeatures %>% rownames_to_column, by="flybase")
  cell_cycle_drosophila %>% pull(rowname) %>% split(cell_cycle_drosophila %>% pull(phase))
}

fpkm_quartile_factor_make_cutoffs <- function(log_data) {
  apply(log_data, 2, \(v) v %>% subset(. >= cutoff) %>% quantile(c(1/3, 2/3)), simplify=F)
}

fpkm_third_density <- function(
  log_fpkm, cutoff = log(5) / log(10), ylim = c(-5, NA), y_label = bquote(log[10]*"(CPM)"), clusters = c('germline', 'somatic'),
  inter_cutoffs = apply(log_fpkm, 2, \(v) v %>% subset(. >= cutoff) %>% quantile(c(1/3, 2/3)), simplify=F)
) {
  log_fpkm <- log_fpkm %>% subset(rowAlls(is.finite(.)))
  names(dimnames(log_fpkm)) <- c('gene', 'cluster')
  density_data <- apply(
    log_fpkm,
    2,
    \(v) density(v, n=1000),
    simplify = F
  )
  density_cut <- mapply(
    \(name, v) c(-Inf, cutoff, inter_cutoffs[[name]], Inf),
    colnames(log_fpkm),
    apply(log_fpkm, 2, identity, simplify=FALSE)
  )
  density_melt <- mapply(
      \(name, d, density_cut) d %>%
        with(data.frame(x, y, quartile=cut(x, density_cut))) %>%
        as_tibble %>%
        mutate(
          quartile = quartile %>% `levels<-`(
            value = if (length(levels(.)) <= 2) c("off", "on") else paste0('Q', 1:4)
          )
        ) %>%
        group_by(quartile) %>%
        filter(
          sum(between(log_fpkm[, name], min(x), max(x))) > 0
        ),
      colnames(log_fpkm),
      density_data,
      density_cut %>% asplit(2),
      SIMPLIFY = F) %>%
    bind_rows(.id = "cluster") %>%
    subset(between(x, -10, 10))
  density_melt$cluster <- density_melt$cluster %>% factor(sce.clusters$cluster)
  density_melt <- density_melt %>%
    left_join(density_melt %>% group_by(cluster) %>% summarise(max_density = max(y), res = diff(x)[1]), "cluster") %>%
    within(violinwidth <- y / max_density) %>%
    subset(cluster %in% clusters)
  polygon_melt <- density_melt %>%
    group_by(cluster, quartile) %>%
    arrange(cluster, quartile, x) %>%
    reframe(
      x = c(x[1] - res[1], x, rev(x), x[1] - res[1]),
      violinwidth = c(violinwidth[1], violinwidth, rev(-violinwidth), -violinwidth[1])
    )
  polygon_melt %>%
    within(quartile <- quartile %>% factor(rev(sort(unique(.))))) %>%
    ggplot(
      aes(x = as.numeric(cluster) + violinwidth * 0.95 / 2, y = x, fill = quartile)
    ) + geom_polygon(aes(group = interaction(cluster, quartile))) + scale_fill_manual(
      values = sc_quartile_colors %>% replace(1, "#aaaaaa") %>% setNames(NULL) %>% rev
    ) + coord_cartesian(
      NULL, ylim
    ) + scale_x_continuous(
      breaks = c(1, 2),
      labels = sapply(
        clusters,
        \(n) paste0(toupper(substr(n, 1, 1)), substr(n, 2, str_length(n)))
      )
    ) + theme_bw() + theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) + labs(
      x = 'Cluster', y = y_label
    )
}

fpkm_simple_violin <- function(
  data, ylim = c(-2.75, 4.5), ygray = log(5) / log(10)
) {
  gg <- melt(data, variable.name = "cluster") %>%
    subset(value > ylim[1] - 1) %>%
    # We would like to set bounds on ydensity... Just add some pseudo-counts
    # that are out of bounds to be displayed, so that the violins will have the
    # same bounds.
    rbind(tibble(cluster = colnames(data), value = ylim[1] - 1)) %>%
    rbind(tibble(cluster = colnames(data), value = 10)) %>%
    ggplot(aes(cluster, value, fill=cluster)) +
    geom_violin(linewidth = NA) +
    scale_fill_manual(
      values = unlist(
        cluster_colors[c("spermatocyte", "somaticprecursor", "muscle")],
        use.names = FALSE
      )
    )
  data <- ggplot_build(gg)$data[[1]]
  data %>%
    subset(y >= ygray) %>%
    ggplot(
      aes(x, y, violinwidth=violinwidth, fill=fill, group=group)
    ) +
    geom_violin(
      data=data %>% subset(y <= ygray + 0.02),
      stat="identity",
      fill="#cccccc",
      color="transparent"
    ) +
    geom_violin(stat="identity", color="transparent") +
    coord_cartesian(c(0.4, 3.6), ylim, expand=F) +
    scale_x_continuous(
      breaks = 1:3,
      minor_breaks = NULL,
      labels = levels(gg$data$cluster)
    ) +
    scale_y_continuous(breaks = seq(-2, 4)) +
    scale_fill_identity() +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "none"
    ) +
    labs(x = 'Cluster', y = bquote(log[10]*"(CPM)"))
}

gene_group_bar_plot <- function(
  quartile.factor_Germline, quartile.factor_Somatic, Upd_cpm
) {
  data <- rbind(
    tibble(
      label = "Both - GSC Level",
      group = quartile.factor_Germline %>% subset(
        Upd_cpm[, "germline"] >= 5 & Upd_cpm[, "somatic"] >= 5
      )
    ),
    tibble(
      label = "Both - CySC Level",
      group = quartile.factor_Somatic %>% subset(
        Upd_cpm[, "germline"] >= 5 & Upd_cpm[, "somatic"] >= 5
      )
    ),
    tibble(
      label = "GSC Exclusive",
      group = quartile.factor_Germline %>% subset(
        Upd_cpm[, "germline"] >= 5 & Upd_cpm[, "somatic"] < 5
      )
    ),
    tibble(
      label = "CySC Exclusive",
      group = quartile.factor_Somatic %>% subset(
        Upd_cpm[, "germline"] < 5 & Upd_cpm[, "somatic"] >= 5
      )
    )
  ) %>%
    mutate(
      label = label %>% factor(unique(.)),
      group = group %>% factor() %>% recode(Q2="low", Q3="medium", Q4="high")
    )
  data <- as_tibble(
    data.frame(
      label = ifelse(
        Upd_cpm[, "germline"] >= 5,
        ifelse(
          Upd_cpm[, "somatic"] >= 5,
          "Both",
          "GSC Exclusive"
        ),
        ifelse(
          Upd_cpm[, "somatic"] >= 5,
          "CySC Exclusive",
          NA
        )
      ) %>%
        setNames(NULL) %>%
        factor(c("Both", "GSC Exclusive", "CySC Exclusive")),
      celltype = factor(
        rep(c("GSC", "CySC"), each = length(quartile.factor_Germline)),
        c("GSC", "CySC")
      ),
      group = c(
        quartile.factor_Germline,
        quartile.factor_Somatic
      )
    )
  ) %>%
    subset(!is.na(label)) %>%
    mutate(
      label = label %>% factor(unique(.)),
      group = group %>% factor() %>% recode(Q1="off", Q2="low", Q3="medium", Q4="high")
    )
  data %>%
    ggplot(aes(celltype, fill = group)) +
    facet_grid(cols=vars(label), switch="x") +
    geom_bar(position = position_fill(rev = TRUE)) +
    scale_fill_manual(values = sc_quartile_colors %>% setNames(NULL)) +
    scale_y_continuous(labels = percent) +
    coord_cartesian(
      c(0.4, 2.6),
      c(0, 1),
      expand = FALSE
    ) +
    labs(
      x = "Gene Classification",
      y = "% of Selected Genes"
    ) +
    theme(
      aspect.ratio = 3,
      axis.line = element_line(color = "black", linewidth = 0.5),
      panel.border = element_rect(color = NA),
      panel.spacing = unit(-0.01, "cm"),
      strip.background = element_rect(fill = NA, color = NA),
      strip.placement = "outside",
    )
}

write_Upd_sc_cell_cycle_phases = function(Upd_sc, cell_cycle_drosophila, metafeatures) {
  Upd_sc = Upd_sc %>% NormalizeData
  cell_cycle = load_cell_cycle_score_drosophila(cell_cycle_drosophila, metafeatures)
  Upd_sc = Upd_sc %>% CellCycleScoring(s. = cell_cycle$S, g2m. = cell_cycle$`G2/M`)
  Upd_sc$Phase <- Upd_sc$Phase %>% factor(c("G1", "S", "G2M"))
  append(
    list(
      total = table(Idents(Upd_sc), Upd_sc$Phase) %>%
        `/`(rowSums(.)) %>%
        round(3)
    ),
    sapply(
      c("nos.1", "nos.2", "tj.1", "tj.2"),
      \(n) table(
        Idents(Upd_sc)[Upd_sc$batch == n],
        Upd_sc$Phase[Upd_sc$batch == n]
      ) %>% `/`(rowSums(.)) %>% round(3),
      simplify = F
    ) %>%
      # Prettify the sce.data batch names.
      setNames(str_extract(sce.data$tenx_path, "/([^/]+)/", group=1))
  ) %>%
    sapply(
      # table to data frame of the table's entries
      \(t) as.data.frame(matrix(t, nrow=nrow(t), ncol=ncol(t), dimnames=dimnames(t))) %>%
        rownames_to_column("cluster"),
      simplify = FALSE
    )
}

plot_genes_on_chrs <- function(
  regression_table,
  coef,
  assay.data.sc,
  do_smooth=FALSE
) {
  assay.data.sc <- assay.data.sc %>%
    read.csv(row.names = 1) %>%
    subset(
      rownames(.) %in% rownames(regression_table)
      & chr %in% names(chr.lengths)
    ) %>%
    mutate(
      pos = start + c(
        `2L`=0, `2R`=chr.lengths["2L"], `3L`=0, `3R`=chr.lengths["3L"],
        `4`=0, `X`=0, `Y`=0
      )[chr],
      chromosome = chr %>%
        factor %>%
        recode(
          `2L`="2", `2R`="2", `3L`="3", `3R`="3"
        )
    ) %>%
    arrange(chromosome, pos)
  print_data <- cbind(
    assay.data.sc,
    log2FC = regression_table[rownames(assay.data.sc), coef] / log(2)
  )
  print_data$log2FC <- print_data$log2FC %>%
    split(print_data$chromosome) %>%
    sapply(\(v) rollmean(v, pmin(10, length(v) * 0.1), c("extend", "extend", "extend"))) %>%
    do.call(c, .)
  ggplot(
    print_data,
    aes(pos, log2FC, group=chromosome)
  ) + facet_grid(
    . ~ chromosome,
    scales = "free_x", space = "free_x"
  ) + geom_line(
    data = cross_join(
      data.frame(pos = c(-Inf, Inf), log2FC = 0),
      data.frame(chromosome = levels(print_data$chromosome))
    ),
    color = "red"
  ) + (if (do_smooth) geom_smooth(method="loess") else geom_line()) + coord_cartesian(
    NULL, c(-1.25, 1.25), expand=F
  ) + scale_x_continuous(
    labels = NULL,
    breaks = c(1, seq(2000000, 50000000, by=2000000)),
    minor_breaks = seq(1000000, 51000000, by=2000000)
  )
}

plot_chr_ratio_on_clusters <- function(Upd_sc) {
  X_genes <- (Upd_sc[["RNA"]]@meta.data$chr == "X") %>%
    replace_na(0) %>%
    `/`(sum(.))
  AA_genes <- (Upd_sc[["RNA"]]@meta.data$chr %in% c("2L", "2R", "3L", "3R", "4")) %>%
    replace_na(0) %>%
    `/`(sum(.))
  feature_matrix <- t(GetAssayData(Upd_sc, assay = "RNA", layer = "counts"))
  Upd_sc$X_AA <- (
    feature_matrix %*% X_genes
    /
    feature_matrix %*% AA_genes
  ) %>%
    as.data.frame %>%
    rownames_to_column %>%
    deframe
  ggplot(
    FetchData(Upd_sc, c("ident", "X_AA", "batch")),
    aes(ident, X_AA)
  ) + geom_segment(
    aes(xend = xend, yend = yend),
    data.frame(ident = -Inf, xend = Inf, X_AA = 1, yend = 1),
    linetype = "dashed"
  ) + geom_boxplot(
    aes(fill = ident), outlier.shape=NA
  ) + scale_fill_manual(
    values = cluster_colors, guide = NULL
  ) + theme_bw() + theme(
    aspect.ratio = 4/3,
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 12)
  ) + coord_cartesian(
    NULL,
    c(0.45, 1.85)
  ) + scale_y_continuous(
    breaks = c(0.5, 1, 1.5)
  ) + labs(
    title = "X ratio per cell: avg_X(Gene UMI) / avg_AA(Gene UMI)",
    x = NULL, y = NULL
  )
}

# Cell cycle scoring ----
apply_cell_cycle_score <- function(Upd_sc, cell_cycle_drosophila, assay.data.sc) {
  cell_cycle <- load_cell_cycle_score_drosophila(cell_cycle_drosophila, assay.data.sc)
  Upd_sc %>%
    NormalizeData %>%
    CellCycleScoring(s. = cell_cycle$S, g2m. = cell_cycle$`G2/M`)
}

# Standardize the single-cell regression so that for each cluster, the
# classified phases of the cells are 33% G1, 33% S, 33% G2M (uniformly).
subsample_ident_normalize_phase <- function(df) {
  df %>%
    group_by(ident, batch) %>%
    mutate(
      n = min(table(phase))
    ) %>%
    group_by(ident, batch, phase) %>%
    dplyr::slice(sample(length(rowname), n[1])) %>%
    arrange(match(rowname, df$rowname))
}

plot_volcano_apeglm <- function(
  Upd_regression_somatic, log2Threshold = 1.5,
  color_column = NULL
) {
  quant = tibble(
    rowname = rownames(Upd_regression_somatic$map),
    log2FC = -Upd_regression_somatic$map[,2] / log(2),
    log10.svalue = log(Upd_regression_somatic$svalue[match(rownames(Upd_regression_somatic$map), rownames(Upd_regression_somatic$svalue)),1]) / log(10)
  ) %>%
    arrange(runif(nrow(.))) %>%
    mutate(
      color = if (is.character(color_column))
        color_column[rowname]
      else
        ifelse(
          log2FC >= log2Threshold & log10.svalue <= -4,
          chic_line_track_colors$germline,
          ifelse(
            log2FC <= -log2Threshold & log10.svalue <= -4,
            chic_line_track_colors$somatic,
            "#dddddd"
          )
        )
    )

  ggplot(
    quant, aes(log2FC, log10.svalue, color=color)
  ) + geom_point(size=0.1) + geom_line(
    aes(group=group),
    tribble(
      ~log2FC, ~log10.svalue, ~group,
      -Inf, -4, "1",
      -log2Threshold, -4, "1",
      -log2Threshold, -Inf, "1",
      log2Threshold, -Inf, "2",
      Inf, -4, "2",
      log2Threshold, -4, "2"
    ),
    color = "#8a2311",
    linetype = "longdash"
  ) + scale_color_identity() + scale_x_continuous(
    limits=c(-1,1) * 7.5, oob=scales::squish, expand=c(0,0)
  ) + scale_y_reverse(
    limits=c(0, -50), oob=scales::squish, expand=rep(0.01,2)
  ) + theme_bw() + theme(aspect.ratio = 3/2)
}

l2fc_bar_plot <- function(l2fc_data) {
  l2fc_data <- l2fc_data %>% sort(dec=T)
  ggplot(
    enframe(l2fc_data) %>% cbind(y = seq(nrow(.))),
    aes(y, value, fill=factor(y, y))
  ) + geom_bar(
    stat="identity"
  ) + geom_text(
    aes(label=name),
    hjust = l2fc_data > 0,
    size = 3
  ) + coord_flip() + scale_x_reverse(
    expand = rep(0.01, 2),
    labels = NULL
  ) + scale_y_continuous(
    expand = rep(0.02, 2),
    breaks = \(v) union(-2.5, scales::pretty_breaks()(v))
  ) + scale_fill_manual(
    values = hcl(
      -30 + seq(0, 360, length.out=length(l2fc_data)+1)[-(length(l2fc_data)+1)],
      60,
      90
    ),
    guide = NULL
  ) + labs(
    y = bquote(log[2]*"(CySC/GSC)"), x = NULL
  ) + theme_bw() + theme(
    panel.grid = element_blank()
  )
}

# Scatter plot of the main regression L2FC values (asynchronous, cell phase not
# considered) vs G1-phase L2FC values (remove effect of cell phase).
plot_apeglm_async_phase_model <- function(orig_params, new_params_list) {
  data <- tibble(
    all = orig_params$map[,2] %>% replace(is.na(.), 0) / -log(2),
    G1 = sapply(new_params_list, \(lst) lst$map[2]) %>% replace(is.na(.), 0) / -log(2)
  )
  coef <- lm(G1 ~ all, data) %>% coef()
  ggplot(
    data,
    aes(all, G1)
  ) +
    rasterise(geom_point(stroke=NA, size=0.15), dpi=300) +
    geom_abline(
      aes(intercept=intercept, slope=slope, color=color),
      tribble(
        ~intercept, ~slope,
        0, 1,
        coef["(Intercept)"], coef["all"]
      ),
      color=c("blue", "#cc0000"),
      linetype=c("dotted", "solid")
    ) +
    coord_cartesian(
      c(-8, 8),
      c(-8, 8),
      expand=FALSE
    ) +
    annotate(
      "text",
      -6, 6,
      label = str_glue("R = {with(data, round(cor(all, G1), 2))}"),
      hjust = 0
    ) +
    labs(
      x = bquote(log[2]*"(GSC/CySC) Async"),
      y = bquote(log[2]*"(GSC/CySC) 33/33/33 Phase Classification    ")
    ) +
    theme(
      aspect.ratio = 1,
      legend.position = "none"
    )
}
