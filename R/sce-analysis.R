cluster_colors = list(
  germline=hsv(0.35,0.9,0.6),
  somatic=hsv(0.5,0.75,0.95),
  spermatocyte=hcl(86, 93, 85),
  muscle=hcl(17, 112, 58),
  others=hsv(0,0,0.9)
)

sc_quartile_colors = c(
  low = hcl(252, 66, 37),
  hcl(200, 37, 60),
  hcl(137, 37, 80),
  high = hcl(77, 87, 87)
)

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
  cbind(Upd_sc@meta.data, Upd_sc[['umap']]@cell.embeddings, ident=Idents(Upd_sc)) %>% ggplot(
    aes(UMAP_1, UMAP_2, color=ident)
  ) + rasterize(geom_point(
    size = 0.01
  ), dpi=120, scale=0.5) + scale_color_manual(
    values=c(
      germline=hsv(0.35,0.9,0.6),
      somatic=hsv(0.5,0.75,0.95),
      spermatocyte=hcl(86, 93, 85),
      muscle=hcl(17, 112, 58),
      others=hsv(0,0,0.9)),
    guide=guide_legend(title=NULL, override.aes = list(size=4))
  ) + scale_x_continuous(
    breaks=c(-10,0,10)
  ) + scale_y_continuous(breaks=c(-7,0,10)) + theme_cowplot()
  filenames <- NULL
  filenames <- c(paste(figures_dir, 'UMAP.svg', sep='/'), filenames)
  ggsave(filenames[1], width=8, height=4.5)

  for (gene in c('vas','tj','lncRNA:roX2','Mst87F','soti','sunz','Act57B')) {
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
    scale_color_gradientn(
      colors = c(
        low = hcl(252, 66, 37),
        hcl(200, 37, 60),
        hcl(137, 37, 80),
        hcl(77, 87, 87)
      ),
      limits = c(0, gene.max.intensity), oob = squish
    ) + scale_x_continuous(
      breaks=c(-10,0,10)
    ) + scale_y_continuous(breaks=c(-7,0,10)) + theme_cowplot(
    ) + labs(
      tag=gene.save
    ) + theme(
      plot.tag.position = c(0.1, 0.94)
    )
    filenames <- c(paste0(figures_dir, '/UMAP', gene.save, '.svg'), filenames)
    ggsave(filenames[1], width=8, height=4.5)
  }
  filenames
}

load_cell_cycle_score = function(Upd_sc, supplemental_data_xlsx) {
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

Upd_pca_figures = function(figures_dir, Upd_sc) {
  Upd_sc = Upd_sc %>% NormalizeData
  cell_cycle = load_cell_cycle_score(Upd_sc, supplemental_data_xlsx)
  Upd_sc = Upd_sc %>% CellCycleScoring(s. = cell_cycle$S, g2m. = cell_cycle$`G2/M`)
  Upd_subset = Upd_sc[, !is.na(Upd_sc[['pca.subset']]@cell.embeddings[,1])]

  meta.data = cbind(
    Upd_subset@meta.data,
    Mst84Da=Upd_subset[['integrated']]@scale.data['Mst84Da',],
    MtnA=Upd_subset[['integrated']]@scale.data['MtnA',],
    Hsromega=Upd_subset[['integrated']]@scale.data['lncRNA:Hsromega',]
  )
  mn.expl = function(mn) (
    colVars(mn$fitted.values, useNames = T) / (
      colVars(mn$fitted.values, useNames = T)
      + colVars(mn$residuals, useNames = T)
    )
  )
  Upd_sc@meta.data = Upd_sc@meta.data %>% cbind(logUMI = log(Upd_sc$nCount_RNA) / log(2))
  expl.stats = sapply(
    c('Mst84Da','MtnA','lncRNA:Hsromega','lncRNA:roX2','RpL22-like','pct.mito','pct.ribo','logUMI'),
    \(name) mn.expl(
      manova(
        pca ~ gene,
        list(
          pca = Upd_subset[['pca.subset']]@cell.embeddings[,1:5],
          gene = append(
            list(
              pct.mito=Upd_subset$pct.mito,
              pct.ribo=Upd_subset$pct.ribo,
              logUMI=log(Upd_subset$nCount_RNA)
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
    ) / colVars(Upd_subset[['pca.subset']]@cell.embeddings[,1:5])
  ) %>% as.data.frame
  DimPlot(Upd_subset, red='pca.subset')

  bar.data = expl.stats %>% with(
    data.frame(
      component=paste0('PC_', 1:5),
      SSB=between.sum.squares,
      `within(germline)`=germline.var,
      `within(somatic)`=somatic.var,
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
    'within(germline)', 'within(somatic)', 'SSB'
  )
  bar.display.data = data.frame(
    variance=character(),
    group=character()
  ) %>% add_row(
    variance='SSB', group='clusters'
  ) %>% add_row(
    variance='within(germline)', group='clusters'
  ) %>% add_row(
    variance='within(somatic)', group='clusters'
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

  pc_expl_plot <- function(pc) ggplot(
    bar.data %>% subset(component == pc) %>% left_join(bar.display.data, 'variance'),
    aes(fill=variance, pattern_fill=variance, x=value, y=group)
  ) + geom_bar_pattern(
    position='stack', stat='identity',
    pattern_color=NA,
    pattern_angle=45,
    pattern_spacing=0.05,
    pattern_density=0.5,
    pattern_key_scale_factor=0.25
  ) + scale_fill_manual(
    values=c(
      `within(germline)`=cluster_colors$germline,
      `within(somatic)`=cluster_colors$somatic,
      `SSB`=cluster_colors$germline,
      roX2=hcl(298, 100, 65),
      Hsromega=hcl(346, 100, 65),
      Mst84Da=hcl(50, 100, 65)
    )
  ) + scale_pattern_fill_manual(
    values=c(
      `within(germline)`=cluster_colors$germline,
      `within(somatic)`=cluster_colors$somatic,
      `SSB`=cluster_colors$somatic,
      roX2=hcl(298, 100, 65),
      Hsromega=hcl(346, 100, 65),
      Mst84Da=hcl(50, 100, 65)
    )
  ) + scale_x_continuous(
    expand=c(0,0), labels=percent
  ) + scale_y_discrete(
    limits=rev
  ) + labs(
    x = '% variance explained in PC (Pearson R)',
    y = ''
  )

  # "somatic" level comes after "germline" level, want this sign to be 1
  Upd_subset[['pca.subset']]@cell.embeddings[,1] = (
    Upd_subset[['pca.subset']]@cell.embeddings[,1]
    * sign(
      coef(
        lm(
          idents ~ pcasubset_1,
          cbind(idents=Idents(Upd_subset), Upd_subset[['pca.subset']]@cell.embeddings[,'pcasubset_1',drop=F])
        )
      )['pcasubset_1']
    )
  )
  DefaultAssay(Upd_subset) = 'RNA'
  pairs_borders = theme(
    plot.background = element_rect(color='black', linewidth=1)
  )
  ggarrange(
    pc_expl_plot('PC_1')
    + theme(plot.background = element_rect(color='black', linewidth=1))
    ,
    pc_expl_plot('PC_2')
    + theme(plot.background = element_rect(color='black', linewidth=1))
  )
  ggarrange(
    pc_expl_plot('PC_1')
    + theme_bw() + pairs_borders,
    DimPlot(Upd_subset, red='pca.subset', com=F)[[1]]
    + scale_color_manual(values=unlist(cluster_colors[1:2]))
    + pairs_borders,
    DimPlot(Upd_subset, c(1,3), red='pca.subset', com=F)[[1]]
    + scale_color_manual(values=unlist(cluster_colors[1:2]))
    + pairs_borders,

    FeaturePlot(Upd_subset, 'lncRNA:roX2', red='pca.subset', com=F)[[1]]
    + scale_color_viridis_c(option='magma')
    + pairs_borders,
    pc_expl_plot('PC_2')
    + pairs_borders,
    DimPlot(Upd_subset, c(2,3), red='pca.subset', com=F)[[1]]
    + scale_color_manual(values=unlist(cluster_colors[1:2]))
    + pairs_borders,
    FeaturePlot(Upd_subset, 'lncRNA:Hsromega', red='pca.subset', dim=c(1,3), com=F)[[1]]
    + scale_color_viridis_c(option='magma')
    + pairs_borders,
    FeaturePlot(Upd_subset, 'Mst84Da', red='pca.subset', dim=c(2,3), com=F)[[1]]
    + scale_color_viridis_c(option='magma')
    + pairs_borders,
    pc_expl_plot('PC_3')
    + pairs_borders,
    nrow=3,ncol=3
  )
  ggsave(
    paste(figures_dir, 'Germline-Somatic-Pairs.png', sep='/'),
    width=15, height=10, dpi=100,
    bg='white')
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