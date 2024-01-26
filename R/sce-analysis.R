cluster_colors = list(
  germline=hsv(0.35,0.9,0.6),
  somatic=hsv(0.5,0.75,0.95),
  spermatocyte=hcl(86, 93, 85),
  differentiated=hcl(86, 20, 75),
  muscle=hcl(17, 112, 58),
  others=hsv(0,0,0.9)
)

cluster_contrast = list(
  germline = hcl(124, 55, 38),
  somatic = hcl(214, 37, 45)
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
  low = hcl(252, 66, 50),
  hcl(200, 37, 60),
  hcl(137, 37, 80),
  high = hcl(77, 87, 87)
)

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
    aes(UMAP_1, UMAP_2, color=ident)
  ) + rasterize(geom_point(
    size = 0.01
  ), dpi=120, scale=0.5) + scale_color_manual(
    values=unlist(cluster_colors),
    guide=guide_legend(title=NULL, override.aes = list(size=4))
  ) + scale_x_continuous(
    breaks=c(-10,0,10)
  ) + scale_y_continuous(breaks=c(-7,0,10)) + theme_cowplot()
  filenames <- NULL
  filenames <- c(paste(figures_dir, 'UMAP.svg', sep='/'), filenames)
  ggsave(filenames[1], width=8, height=4.5)

  for (gene in c('RpL22-like','vas','tj','lncRNA:roX1','lncRNA:roX2','Mst87F','soti','sunz','w-cup','Act57B')) {
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
      aes(UMAP_1, UMAP_2, color=LogNormalize)
    ) + rasterize(geom_point(
      size = 0.01
    ), dpi=240) + # + scale_color_viridis_c(
      # option='magma', limits=c(0,gene.max.intensity), oob=squish
    scale_color_viridis_c(
      begin = 0.2,
      limits = c(0, gene.max.intensity), oob = squish
    ) + scale_x_continuous(
      breaks=c(-10,0,10)
    ) + scale_y_continuous(breaks=c(-7,0,10)) + theme_cowplot(
    ) + labs(
      tag=gene.save
    ) + theme(
      plot.tag.position = c(0.1, 0.94)
    )
    filenames <- c(paste0(figures_dir, '/UMAP-', gene.save, '.svg'), filenames)
    ggsave(filenames[1], width=8, height=4.5)
  }
  filenames
}

# Plots for the replicates, using our filter by pct.mito and pct.ribo but no
# clustering or downstream doublet removal.
gene_pseudobulk_fpkm <- function(tx_files, seurats, seurat_names, gtf_path, metadata_path, metafeatures_path) {
  metadata <- read.csv(metadata_path, row.names = 1)
  size_factors <- seurat_names %>%
    sapply(\(n) metadata %>% subset(batch == n & nCount_RNA_filter == "nCount_RNA_pass") %>% pull(nCount_RNA) %>% sum)
  cpm_table <- pseudobulk_cpm(data.frame(batch = 'nos.1', tx_file = tx_files), data.frame(batch = 'nos.1', size_factor = size_factors), gtf_path)
  fpkm_table <- cpm_table %>% tx_cpm_to_fpkm(metafeatures_path)
  data.frame(
    fpkm = fpkm_table[rownames(seurats[[1]][['RNA']]), 1],
    percent_expressed = rowMeans(
      seurats %>%
        mapply(
          \(seur, n) seur[['RNA']]@counts[
            ,
            metadata %>%
              subset(batch == n & nCount_RNA_filter == "nCount_RNA_pass") %>%
              rownames %>%
              str_replace(paste0(n, "_"), "")
          ] %>%
            `!=`(0) %>%
            rowMeans,
          .,
          seurat_names)
    )
  )
}

gene_cluster_fpkm <- function(Upd_fpkm, seurats, ident.quantify, metadata_path) {
  metadata <- read.csv(metadata_path, row.names = 1)
  data.frame(
    fpkm = Upd_fpkm[rownames(seurats[[1]]), ident.quantify],
    percent_expressed = rowMeans(
      seurats %>%
        mapply(
          \(seur, n) seur[['RNA']]@counts[
            ,
            metadata %>%
              subset(batch == n & ident == ident.quantify & nCount_RNA_filter == "nCount_RNA_pass") %>%
              rownames %>%
              str_replace(paste0(n, "_"), "")
          ] %>%
            `!=`(0) %>%
            rowMeans,
          .,
          names(.))
    )
  )
}

load_cell_cycle_score_orthologs = function(Upd_sc, supplemental_data_xlsx) {
  cell_cycle_features = read_excel(
    supplemental_data_xlsx, sheet='Gene Sets Used in Analysis'
  )
  cell_cycle_features = as.list(
    cell_cycle_features
  ) %>% subset(
    grepl('^(S|G2/M)$', names(.))
  ) %>% sapply(
    \(v) data.frame(input_gene = na.omit(v), row.names = na.omit(v)),
    simplify=F
  ) %>% bind_rows(.id = 'feature_set')
  cell_cycle_features = cell_cycle_features %>% convert_orthologs(
    input_species='human',
    output_species='dmelanogaster'
  )
  cell_cycle_features = cell_cycle_features %>% subset(
    rownames(.) %in% rownames(Upd_sc)
  )
  rownames(cell_cycle_features) %>% split(cell_cycle_features$feature_set)
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

plot_Upd_pca_components = function(Upd_sc, cell_cycle) {
  Upd_sc = Upd_sc %>% NormalizeData
  Upd_sc = Upd_sc %>% CellCycleScoring(s. = cell_cycle$S, g2m. = cell_cycle$`G2/M`)
  Upd_subset = Upd_sc[, !is.na(Upd_sc[['pca.subset']]@cell.embeddings[,1])]

  meta.data = cbind(
    Upd_subset@meta.data,
    Mst84Da=Upd_subset[['RNA']]@data['Mst84Da',],
    MtnA=Upd_subset[['RNA']]@data['MtnA',],
    Hsromega=Upd_subset[['RNA']]@data['lncRNA:Hsromega',]
  )
  mn.expl = function(mn) (
    colVars(mn$fitted.values, useNames = T) / (
      colVars(mn$fitted.values, useNames = T)
      + colVars(mn$residuals, useNames = T)
    )
  )
  Upd_sc@meta.data = Upd_sc@meta.data %>% cbind(logUMI = log(Upd_sc$nCount_RNA) / log(2))
  expl.stats = sapply(
    c('Mst84Da','MtnA','lncRNA:Hsromega','lncRNA:roX2','AGO3','pct.mito','pct.ribo','logUMI','phase','driver','batch'),
    \(name) mn.expl(
      manova(
        pca ~ gene,
        list(
          pca = Upd_subset[['pca.subset']]@cell.embeddings[,1:5],
          gene = append(
            list(
              pct.mito=Upd_subset$pct.mito,
              pct.ribo=Upd_subset$pct.ribo,
              logUMI=log(Upd_subset$nCount_RNA),
              phase=model.matrix(~ 0 + Phase, Upd_subset@meta.data),
              driver=factor(ifelse(grepl("nos", Upd_subset$batch), "nos", "tj")),
              batch=factor(Upd_subset$batch)
            ),
            (Upd_subset[['integrated']]@scale.data
            %>% subset(rownames(.) == name) %>% t %>% as.data.frame)
          )[[name]]
        )
      )
    )
  ) %>% cbind(
    between.sum.squares = mn.expl(
      manova(
        pca ~ ident,
        list(
          pca = Upd_subset[['pca.subset']]@cell.embeddings[,1:5],
          ident=Idents(Upd_subset)
        )
      )
    ),
    germline.var = colVars(
      Upd_subset[['pca.subset']]@cell.embeddings[,1:5]
      %>% subset(Idents(Upd_subset) == 'germline')
    ) * mean(
      Idents(Upd_subset) == 'germline'
    ) / colVars(Upd_subset[['pca.subset']]@cell.embeddings[,1:5]),
    somatic.var = colVars(
      Upd_subset[['pca.subset']]@cell.embeddings[,1:5]
      %>% subset(Idents(Upd_subset) == 'somatic')
    ) * mean(
      Idents(Upd_subset) == 'somatic'
    ) / colVars(Upd_subset[['pca.subset']]@cell.embeddings[,1:5]),
    germline.phase.var = mn.expl(manova(pca ~ Phase, list(pca = Upd_subset[["pca.subset"]]@cell.embeddings[, 1:5], Phase = model.matrix(~ Phase, Upd_subset@meta.data)), subset = Idents(Upd_subset) == "germline"))
    * mean(Idents(Upd_subset) == 'germline'),
    somatic.phase.var = mn.expl(manova(pca ~ Phase, list(pca = Upd_subset[["pca.subset"]]@cell.embeddings[, 1:5], Phase = model.matrix(~ Phase, Upd_subset@meta.data)), subset = Idents(Upd_subset) == "germline"))
    * mean(Idents(Upd_subset) == 'somatic')
  ) %>% as.data.frame
  DimPlot(Upd_subset, red='pca.subset')

  bar.data = expl.stats %>% with(
    data.frame(
      component=paste0('PC_', 1:5),
      `SSB(clusters)`=between.sum.squares,
      `within(germline)`=germline.var,
      `within(somatic)`=somatic.var,
      `SSB(G1/S/G2M) germline`=germline.phase.var,
      `SSB(G1/S/G2M) somatic`=somatic.phase.var,
      roX2=`lncRNA:roX2`,
      Hsromega=`lncRNA:Hsromega`,
      Mst84Da=Mst84Da,
      check.names=F
    )
  ) %>% reshape2::melt(
    id.vars = 'component',
    variable.name = 'variance'
  )
  bar.data$variance = bar.data$variance %>% factor(., unique(.)) %>% fct_relevel(
    'within(germline)', 'within(somatic)', 'SSB(clusters)'
  )
  bar.display.data = data.frame(
    variance=character(),
    group=character()
  ) %>% add_row(
    variance='SSB(clusters)', group='clusters'
  ) %>% add_row(
    variance='within(germline)', group='clusters'
  ) %>% add_row(
    variance='within(somatic)', group='clusters'
  ) %>% add_row(
    variance='SSB(G1/S/G2M) germline', group='Phase'
  ) %>% add_row(
    variance='SSB(G1/S/G2M) somatic', group='Phase'
  )
  bar.display.data = bar.display.data %>% rbind(
    sapply(
      setdiff(unique(bar.data$variance), bar.display.data$variance),
      \(n) c(variance=n, group=n)
    ) %>% t
  )
  # Maintain order of appearance in bar.display.data
  bar.display.data$group = bar.display.data$group %>% factor(
    ., unique(.)
  )

  pc_expl_bar_values1 <- c(
    `within(germline)`=cluster_bar_color1$germline,
    `within(somatic)`=cluster_bar_color1$somatic,
    `SSB(clusters)`=cluster_bar_color1$germline,
    `SSB(G1/S/G2M) somatic`=cluster_bar_color2$somatic,
    `SSB(G1/S/G2M) germline`=cluster_bar_color2$germline,
    roX2=hcl(42, 96, 68),
    Hsromega=hcl(13, 138, 49),
    Mst84Da=hcl(275, 114, 49)
  )
  pc_expl_bar_values2 <- c(
    `within(germline)`=cluster_bar_color1$germline,
    `within(somatic)`=cluster_bar_color1$somatic,
    `SSB(clusters)`=cluster_bar_color1$somatic,
    `SSB(G1/S/G2M) somatic`=cluster_bar_color2$somatic,
    `SSB(G1/S/G2M) germline`=cluster_bar_color2$germline,
    roX2=hcl(42, 96, 68),
    Hsromega=hcl(13, 138, 49),
    Mst84Da=hcl(275, 114, 49)
  )
  pc_expl_plot <- function(pc) ggplot(
    # bar.data %>% subset(component == pc) %>% left_join(bar.display.data, 'variance') %>%
    # within(variance <- variance %>% factor(c("within(germline)", "within(somatic)", "SSB(clusters)", "SSB(G1/S/G2M germline)", "SSB(G1/S/G2M) somatic", "roX2", "Hsromega", "Mst84Da")))
    bar.display.data %>% left_join(bar.data %>% subset(component == pc), by = "variance") %>%
    within(variance <- variance %>% factor(c("within(germline)", "within(somatic)", "SSB(clusters)", "SSB(G1/S/G2M) germline", "SSB(G1/S/G2M) somatic", "roX2", "Hsromega", "Mst84Da")))
    ,
    aes(fill=variance, pattern_fill=variance, x=value, y=group)
  ) + geom_bar_pattern(
    # Odd quirk in horizontal bars. We want the reverse order to stack from left to right.
    position=position_stack(reverse = T),
    stat='identity',
    pattern_color=NA,
    pattern_angle=45,
    pattern_spacing=0.05,
    pattern_density=0.5,
    pattern_key_scale_factor=0.25
  ) + scale_fill_manual(
    values = pc_expl_bar_values1,
    labels = names(pc_expl_bar_values1) %>%
      replace(. == "Hsromega", "Hsr\u03C9")
  ) + scale_pattern_fill_manual(
    values = pc_expl_bar_values2,
    labels = names(pc_expl_bar_values1) %>%
      replace(. == "Hsromega", "Hsr\u03C9")
  ) + scale_x_continuous(
    expand=c(0,0), labels=percent
  ) + scale_y_discrete(
    limits=rev, labels=\(v) v %>% replace(. == "Hsromega", "Hsr\u03C9")
  ) + labs(
    x = '% variance explained in PC',
    y = ''
  ) + theme_bw() + theme(
    # Might want to edit figure to delete the legend and blow up the bar chart.
    # Make the bar chart not overly tall.
    plot.margin = margin(t = 35, b = 35, l = 5.5, r = 5.5),
    legend.margin = margin(l = 15, t = 5.5, b = 5.5, r = 5.5),
    legend.background = element_rect(fill = "transparent")
  )

  # "somatic" level comes after "germline" level, want this sign to be 1
  Upd_subset[['pca.subset']]@cell.embeddings[,1] = (
    Upd_subset[['pca.subset']]@cell.embeddings[,1]
    * sign(
      coef(
        lm(
          idents ~ pcasubset_1,
          cbind(data.frame(idents=Idents(Upd_subset) == "somatic"), Upd_subset[['pca.subset']]@cell.embeddings[,'pcasubset_1',drop=F])
        )
      )['pcasubset_1']
    )
  )
  DefaultAssay(Upd_subset) = 'RNA'
  pairs_borders = theme(
    plot.background = element_rect(color='black', linewidth=1)
  )
  ggarrange(
    # Top left.
    pc_expl_plot('PC_1')
    + pairs_borders,
    # Top center.
    DimPlot(Upd_subset, c(2,1), red='pca.subset', com=F)[[1]] %>%
      rasterise(dpi = 160)
    + scale_color_manual(values=unlist(cluster_colors[1:2]))
    + scale_x_continuous(name = NULL, labels = NULL)
    + scale_y_continuous(name = NULL, labels = NULL)
    + theme(plot.margin = margin(l = 70, r = 70, t = 5.5, b = 5.5))
    + pairs_borders,
    # Top right.
    DimPlot(Upd_subset, c(3,1), red='pca.subset', com=F)[[1]] %>%
      rasterise(dpi = 160)
    + scale_color_manual(values=unlist(cluster_colors[1:2]))
    + scale_x_continuous(name = NULL, labels = NULL)
    + scale_y_continuous(name = NULL, labels = NULL)
    + theme(plot.margin = margin(l = 70, r = 70, t = 5.5, b = 5.5))
    + pairs_borders,
    # Center left: a lot of features that we need to explain using PC 1 and 2.
    plot_grid(
      FeaturePlot(Upd_subset, 'AGO3', red='pca.subset', pt.size = 1e-3, com=F)[[1]] %>%
        rasterise(dpi = 160)
      + labs(x = NULL, y = NULL)
      + guides(color = guide_colorbar(barwidth = 0.5, barheight = 3))
      + scale_x_continuous(labels=NULL) + scale_y_continuous(labels=NULL)
      + scale_color_viridis_c(begin=0.2, breaks=pretty_breaks(3)),
      FeaturePlot(Upd_subset, 'lncRNA:roX2', red='pca.subset', pt.size = 1e-3, com=F)[[1]] %>%
        rasterise(dpi = 160)
      + labs(x = NULL, y = NULL, title = "roX2")
      + guides(color = guide_colorbar(barwidth = 0.5, barheight = 3))
      + scale_x_continuous(labels=NULL) + scale_y_continuous(labels=NULL)
      + scale_color_viridis_c(begin=0.2, breaks=pretty_breaks(3)),
      FeaturePlot(Upd_subset, 'lncRNA:Hsromega', red='pca.subset', pt.size = 1e-3, com=F)[[1]] %>%
        rasterise(dpi = 160)
      + labs(x = NULL, y = NULL, title = "Hsr\u03C9")
      + guides(color = guide_colorbar(barwidth = 0.5, barheight = 3))
      + scale_x_continuous(labels=NULL) + scale_y_continuous(labels=NULL)
      + scale_color_viridis_c(begin=0.2, breaks=pretty_breaks(2)),
      DimPlot(Upd_subset, gr='Phase', red='pca.subset', pt.size = 1e-3, com=F)[[1]] %>%
        rasterise(dpi = 160)
      + scale_x_continuous(labels=NULL) + scale_y_continuous(labels=NULL)
      + labs(x = NULL, y = NULL)
    )
    + pairs_borders,
    # Center: PC 2.
    pc_expl_plot('PC_2')
    + pairs_borders,
    # Center right: scatter plot.
    DimPlot(Upd_subset, c(3,2), red='pca.subset', com=F)[[1]] %>%
      rasterise(dpi = 160)
    + scale_color_manual(values=unlist(cluster_colors[1:2]))
    + scale_x_continuous(name = NULL, labels = NULL)
    + scale_y_continuous(name = NULL, labels = NULL)
    + theme(plot.margin = margin(l = 70, r = 70, t = 5.5, b = 5.5))
    + pairs_borders,
    # Bottom left: another PC 1-related feature (germline).
    FeaturePlot(Upd_subset, 'vas', red='pca.subset', dim=c(1,3), com=F)[[1]] %>%
      rasterise(dpi = 160)
    + scale_color_viridis_c(begin=0.2)
    + scale_x_continuous(name = NULL, labels = NULL)
    + scale_y_continuous(name = NULL, labels = NULL)
    + pairs_borders,
    # Bottom center: A PC 3-related feature (unwanted cells in this subset).
    FeaturePlot(Upd_subset, 'Mst84Da', red='pca.subset', dim=c(2,3), com=F)[[1]] %>%
      rasterise(dpi = 160)
    + scale_color_viridis_c(begin=0.2)
    + scale_x_continuous(name = NULL, labels = NULL)
    + scale_y_continuous(name = NULL, labels = NULL)
    + pairs_borders,
    # Bottom right: PC 3 scatter plot.
    pc_expl_plot('PC_3')
    + pairs_borders,
    nrow=3,ncol=3
  )
}

fpkm_quarter_density <- function(Upd_fpkm, output_file, clusters = c('germline', 'somatic')) {
  Upd_fpkm <- Upd_fpkm %>% subset(rowAlls((. > 0) %>% replace(is.na(.), FALSE)))
  names(dimnames(Upd_fpkm)) <- c('gene', 'cluster')
  log_fpkm <- log(Upd_fpkm) / log(10)
  density_data <- apply(
    log_fpkm,
    2,
    \(v) density(v, n=1000),
    simplify = F
  )
  density_cut <- apply(
    log_fpkm,
    2,
    \(v) c(-Inf, quantile(v, c(0.25, 0.50, 0.75)), Inf)
  )
  density_melt <- mapply(
      \(d, density_cut) d %>%
        with(data.frame(x, y, quartile=cut(x, density_cut))) %>%
        within(levels(quartile) <- factor(paste0('Q', 1:4))),
      density_data,
      density_cut %>% asplit(2),
      SIMPLIFY = F) %>%
    bind_rows(.id = "cluster")
  density_melt$cluster <- density_melt$cluster %>% factor(sce.clusters$cluster)
  density_melt <- density_melt %>%
    left_join(density_melt %>% group_by(cluster) %>% summarise(max_density = max(y), res = diff(x)[1]), "cluster") %>%
    within(violinwidth <- y / max_density) %>%
    subset(cluster %in% clusters)
  density_melt %>%
    ggplot(
      aes(cluster, x, height=res, width=y, fill=cluster, group=cluster)
    ) + geom_tile()
  density_melt %>%
    ggplot(
      aes(cluster, x, violinwidth=violinwidth, fill=quartile, group=interaction(quartile, cluster))
    ) + geom_violin(stat='identity')
  density_melt %>%
    ggplot(
      aes(cluster, x, height=res, width=violinwidth * 0.95, fill=quartile)
    ) + rasterise(geom_tile(), dpi=240) + scale_fill_manual(values = sc_quartile_colors %>% setNames(NULL)) + coord_cartesian(
      NULL, c(-5,NA)
    ) + theme_bw() + labs(
      x = 'Cluster', y = bquote(log[10]*"(FPKM)")
    )
  ggsave(output_file, width = 4, height = 4)
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
    )
  ) %>%
    sapply(
      # table to data frame of the table's entries
      \(t) as.data.frame(matrix(t, nrow=nrow(t), ncol=ncol(t), dimnames=dimnames(t))),
      simplify = FALSE
    )
}