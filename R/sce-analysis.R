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
  low = hcl(252, 66, 50),
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
    rename(all_of(dplyr_rename_lookup_gene))
  gene.data <- gene.data[cells, ]
  gene.max.intensity = gene.data %>% group_by(ident) %>% summarise(quantile(LogNormalize, 0.95)) %>% pull(2) %>% max
  gene.max.intensity = c(Mst87F=5, Act57B=4, soti=3, sunz=2, `Amy-d`=2, `scpr-B`=3, ey=1.5)[gene] %>% replace(is.na(.), gene.max.intensity)
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
    rename(all_of(dplyr_rename_lookup_gene))
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

# Plots for bulk expression of genes in each batch.
gene_pseudobulk_cpm <- function(tx_files, seurats, seurat_names, gtf_path, metadata_path, metafeatures_path) {
  metadata <- read.csv(metadata_path, row.names = 1)
  size_factors <- seurat_names %>%
    sapply(\(n) metadata %>% subset(batch == n & nCount_RNA_filter == "nCount_RNA_pass") %>% pull(nCount_RNA) %>% sum)
  cpm_table <- pseudobulk_cpm(
    data.frame(batch = 'nos.1', tx_file = tx_files),
    data.frame(batch = 'nos.1', size_factor = size_factors),
    gtf_path
  )
  data.frame(
    cpm = cpm_table[rownames(seurats[[1]][['RNA']]), 1],
    percent_expressed = rowMeans(
      seurats %>%
        mapply(
          \(seur, n) seur[['RNA']]@layers$counts[
            ,
            metadata %>%
              subset(batch == n & nCount_RNA_filter == "nCount_RNA_pass") %>%
              rownames %>%
              match(Cells(seur))
          ] %>%
            `!=`(0) %>%
            rowMeans,
          .,
          seurat_names)
    ),
    row.names = rownames(seurats[[1]][['RNA']])
  )
}

gene_cluster_cpm <- function(Upd_cpm, seurats, ident.quantify, metadata_path) {
  metadata <- read.csv(metadata_path, row.names = 1)
  percent_expressed <- enframe(seurats) %>%
    rowwise %>%
    summarise(
      gene_percents = subset(
        value,
        cells = metadata %>%
          subset(batch == name & ident == ident.quantify & nCount_RNA_filter == "nCount_RNA_pass") %>%
          rownames
      ) %>%
        GetAssayData %>%
        `!=`(0) %>%
        rowMeans %>%
        list
    ) %>%
    pull(1) %>%
    simplify2array %>%
    rowMeans
  data.frame(
    cpm = Upd_cpm[rownames(seurats[[1]]), ident.quantify],
    percent_expressed
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
  Upd_subset = Upd_sc[, !is.na(Upd_sc[['pcasubset']]@cell.embeddings[,1])]

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
          pca = Upd_subset[['pcasubset']]@cell.embeddings[,1:5],
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
          pca = Upd_subset[['pcasubset']]@cell.embeddings[,1:5],
          ident=Idents(Upd_subset)
        )
      )
    ),
    germline.var = colVars(
      Upd_subset[['pcasubset']]@cell.embeddings[,1:5]
      %>% subset(Idents(Upd_subset) == 'germline')
    ) * mean(
      Idents(Upd_subset) == 'germline'
    ) / colVars(Upd_subset[['pcasubset']]@cell.embeddings[,1:5]),
    somatic.var = colVars(
      Upd_subset[['pcasubset']]@cell.embeddings[,1:5]
      %>% subset(Idents(Upd_subset) == 'somatic')
    ) * mean(
      Idents(Upd_subset) == 'somatic'
    ) / colVars(Upd_subset[['pcasubset']]@cell.embeddings[,1:5]),
    germline.phase.var = mn.expl(manova(pca ~ Phase, list(pca = Upd_subset[["pcasubset"]]@cell.embeddings[, 1:5], Phase = model.matrix(~ Phase, Upd_subset@meta.data)), subset = Idents(Upd_subset) == "germline"))
    * mean(Idents(Upd_subset) == 'germline'),
    somatic.phase.var = mn.expl(manova(pca ~ Phase, list(pca = Upd_subset[["pcasubset"]]@cell.embeddings[, 1:5], Phase = model.matrix(~ Phase, Upd_subset@meta.data)), subset = Idents(Upd_subset) == "germline"))
    * mean(Idents(Upd_subset) == 'somatic')
  ) %>% as.data.frame
  DimPlot(Upd_subset, red='pcasubset')

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
  Upd_subset[['pcasubset']]@cell.embeddings[,1] = (
    Upd_subset[['pcasubset']]@cell.embeddings[,1]
    * sign(
      coef(
        lm(
          idents ~ pcasubset_1,
          cbind(data.frame(idents=Idents(Upd_subset) == "somatic"), Upd_subset[['pcasubset']]@cell.embeddings[,'pcasubset_1',drop=F])
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
    DimPlot(Upd_subset, c(2,1), red='pcasubset', com=F)[[1]] %>%
      rasterise(dpi = 160)
    + scale_color_manual(values=unlist(cluster_colors[1:2]))
    + scale_x_continuous(name = NULL, labels = NULL)
    + scale_y_continuous(name = NULL, labels = NULL)
    + theme(plot.margin = margin(l = 70, r = 70, t = 5.5, b = 5.5))
    + pairs_borders,
    # Top right.
    DimPlot(Upd_subset, c(3,1), red='pcasubset', com=F)[[1]] %>%
      rasterise(dpi = 160)
    + scale_color_manual(values=unlist(cluster_colors[1:2]))
    + scale_x_continuous(name = NULL, labels = NULL)
    + scale_y_continuous(name = NULL, labels = NULL)
    + theme(plot.margin = margin(l = 70, r = 70, t = 5.5, b = 5.5))
    + pairs_borders,
    # Center left: a lot of features that we need to explain using PC 1 and 2.
    plot_grid(
      FeaturePlot(Upd_subset, 'AGO3', red='pcasubset', pt.size = 1e-3, com=F)[[1]] %>%
        rasterise(dpi = 160)
      + labs(x = NULL, y = NULL)
      + guides(color = guide_colorbar(barwidth = 0.5, barheight = 3))
      + scale_x_continuous(labels=NULL) + scale_y_continuous(labels=NULL)
      + scale_color_viridis_c(begin=0.2, breaks=pretty_breaks(3)),
      FeaturePlot(Upd_subset, 'lncRNA:roX2', red='pcasubset', pt.size = 1e-3, com=F)[[1]] %>%
        rasterise(dpi = 160)
      + labs(x = NULL, y = NULL, title = "roX2")
      + guides(color = guide_colorbar(barwidth = 0.5, barheight = 3))
      + scale_x_continuous(labels=NULL) + scale_y_continuous(labels=NULL)
      + scale_color_viridis_c(begin=0.2, breaks=pretty_breaks(3)),
      FeaturePlot(Upd_subset, 'lncRNA:Hsromega', red='pcasubset', pt.size = 1e-3, com=F)[[1]] %>%
        rasterise(dpi = 160)
      + labs(x = NULL, y = NULL, title = "Hsr\u03C9")
      + guides(color = guide_colorbar(barwidth = 0.5, barheight = 3))
      + scale_x_continuous(labels=NULL) + scale_y_continuous(labels=NULL)
      + scale_color_viridis_c(begin=0.2, breaks=pretty_breaks(2)),
      DimPlot(Upd_subset, gr='Phase', red='pcasubset', pt.size = 1e-3, com=F)[[1]] %>%
        rasterise(dpi = 160)
      + scale_x_continuous(labels=NULL) + scale_y_continuous(labels=NULL)
      + labs(x = NULL, y = NULL)
    )
    + pairs_borders,
    # Center: PC 2.
    pc_expl_plot('PC_2')
    + pairs_borders,
    # Center right: scatter plot.
    DimPlot(Upd_subset, c(3,2), red='pcasubset', com=F)[[1]] %>%
      rasterise(dpi = 160)
    + scale_color_manual(values=unlist(cluster_colors[1:2]))
    + scale_x_continuous(name = NULL, labels = NULL)
    + scale_y_continuous(name = NULL, labels = NULL)
    + theme(plot.margin = margin(l = 70, r = 70, t = 5.5, b = 5.5))
    + pairs_borders,
    # Bottom left: another PC 1-related feature (germline).
    FeaturePlot(Upd_subset, 'vas', red='pcasubset', dim=c(1,3), com=F)[[1]] %>%
      rasterise(dpi = 160)
    + scale_color_viridis_c(begin=0.2)
    + scale_x_continuous(name = NULL, labels = NULL)
    + scale_y_continuous(name = NULL, labels = NULL)
    + pairs_borders,
    # Bottom center: A PC 3-related feature (unwanted cells in this subset).
    FeaturePlot(Upd_subset, 'Mst84Da', red='pcasubset', dim=c(2,3), com=F)[[1]] %>%
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

fpkm_quarter_density <- function(
  log_fpkm, ylim = c(-5, NA), y_label = bquote(log[10]*"(CPM)"), clusters = c('germline', 'somatic')
) {
  log_fpkm <- log_fpkm %>% subset(rowAlls(is.finite(.)))
  names(dimnames(log_fpkm)) <- c('gene', 'cluster')
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
      values = sc_quartile_colors %>% setNames(NULL) %>% rev
    ) + coord_cartesian(
      NULL, ylim
    ) + theme_bw() + labs(
      x = 'Cluster', y = y_label
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

analyze_pcasubset_batch_effect <- function(Upd_sc) {
  pcasubset <- Upd_sc[["pcasubset"]]@cell.embeddings
  genotype <- recode(Upd_sc$batch, nos.1='nos', nos.2='nos', tj.1='tj', tj.2='tj')
  batch_effect <- cbind(
    nos_batch = as.numeric(Upd_sc$batch == "nos.1"),
    tj_batch = as.numeric(Upd_sc$batch == "tj.1")
  )
  manova(pcasubset ~ 0 + genotype + batch_effect)
}

analyze_pcasubset_ident <- function(Upd_sc, cell_cycle_drosophila, assay.data.sc) {
  pcasubset <- Upd_sc[["pcasubset"]]@cell.embeddings
  ident <- Idents(Upd_sc)

  Upd_sc = Upd_sc %>% NormalizeData
  cell_cycle = load_cell_cycle_score_drosophila(cell_cycle_drosophila, assay.data.sc)
  Upd_sc = Upd_sc %>% CellCycleScoring(s. = cell_cycle$S, g2m. = cell_cycle$`G2/M`)
  Phase <- Upd_sc$Phase %>% factor(c("G1", "S", "G2M"))

  manova(pcasubset ~ 0 + ident + Phase)
}

analyze_manova <- function(mn, term_names, ncol_mn) {
  # Review the manova display of the terms, source code here:
  # print(stats:::print.aov)
  # The "effects" matrix is a change of basis of the dependent variable (in
  # MANOVA, a matrix of response variables in the columns). The sum of
  # squares of entries (rows) which correspond to model coefficients relating to
  # the term of interest is the sum of squares explained.
  cross_join(
    data.frame(term = term_names),
    data.frame(response = colnames(mn$effects)[seq(ncol_mn)])
  ) %>%
    rowwise %>%
    mutate(
      sse = sum(mn$effects[grepl(term, rownames(mn$effects)), response]^2),
      ssr = sum(mn$fitted.values[, response]^2 + mn$residuals[, response]^2) - sse,
      R2 = sse / (sse + ssr)
    )
}

plot_multiple_umap_data <- function(data) {
  ggplot(
    data,
    aes(umap_1, umap_2, color=ident)
  ) + facet_wrap(
    vars(batch),
    nrow = 2,
    scales = "free"
  ) + rasterize(
    geom_point(
      shape = 20, size = 0.001
    ),
    dpi = 240,
    scale = 0.5
  ) + scale_color_manual(
    values = c(unlist(cluster_colors), doublet="#cccccc"),
    guide = guide_none()
  ) + theme_cowplot() + theme(
    # For facets with no labels indicating what the facet is.
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
}

report_cpm_in_gene_sets <- function(Upd_cpm, gene_sets, output_html, output_pdf) {
  tibble(
    name = names(gene_sets),
    num_genes = sapply(gene_sets, length),
    marker_genes = sapply(
      gene_sets,
      \(v) v %>%
        intersect(c("nos", "vas", "tj", "zfh1", "lncRNA:roX1", "lncRNA:roX2", "AGO1", "AGO3")) %>%
        append(list(sep = ", ")) %>%
        do.call(paste, .)
    )
  ) %>%
    gt %>%
    gtsave(output_html)
  pdf(output_pdf)
  for (n in names(gene_sets)) {
    data <- Upd_cpm[gene_sets[[n]], c("germline", "somatic")]
    names(dimnames(data)) <- c("gene", "cluster")
    data <- log(data) %>% subset(rowAlls(is.finite(.)))
    # melt(data, value.var = "log10CPM") %>%
    print(
      data %>%
        as.data.frame %>%
        ggplot(aes(germline, somatic)) + geom_point(
        ) + geom_density_2d_filled(alpha=0.5) + scale_x_continuous(
          labels = exp,
          breaks = log(c(1, 10, 30, 50, 100, 1000))
        ) + scale_y_continuous(
          labels = exp,
          breaks = log(c(1, 10, 30, 50, 100, 1000))
        ) + coord_cartesian(
          log(c(1, 1000)), log(c(1, 1000)), expand=F
        ) + theme_bw() + guides(
          fill = guide_legend(title = "density")
        ) + labs(
          title = n
        )
    )
  }
  dev.off()
  return(c(output_html, output_pdf))
}