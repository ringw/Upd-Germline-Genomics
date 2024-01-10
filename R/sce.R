read_mt_genome <- function(features_path) {
  features = read.table(features_path, row.names=1)
  features = rownames(features) %>% str_replace('^Dmel_', '')
  suppressWarnings(
    mapIds(
      drosophila2.db,
      features %>% subset(
        mapIds(drosophila2.db, features, 'CHR', 'FLYBASECG') == 'MT'
      ),
      'SYMBOL',
      'FLYBASECG'
    )
  )
}

map_feature_names <- function(feature_names) {
  # SingleCellExperiment does not like a names vector which is a named vector.
  feature_names = setNames(
    mapIds(drosophila2.db, str_replace(feature_names, '^Dmel_', ''), 'SYMBOL', 'FLYBASECG')
    %>% replace(is.na(.) | duplicated(.), names(.)[is.na(.) | duplicated(.)]),
    NULL
  )
  feature_names[feature_names == 'UrbanChen_H3-GFP'] = 'H3-GFP'
  feature_names
}

select_features <- function(paths) {
  feature_names = (
    paths[1]
    %>% paste0('/features.tsv.gz')
    %>% read.table(row.names=1)
    %>% rownames
    %>% map_feature_names
  )
  features_keep = rep(TRUE, length(feature_names)) %>% setNames(feature_names)
  for (tenx_path in paths) {
    matrix_data = read.table(paste0(tenx_path, '/matrix.mtx.gz'), comment='%', header=T)
    colnames(matrix_data)[1] = 'gene'
    matrix_data$gene = factor(matrix_data$gene, seq(1, length(feature_names)))
    gene_num_distinct_cells = table(matrix_data$gene)
    features_keep = features_keep & (gene_num_distinct_cells >= 50)
  }
  features_keep = feature_names[features_keep]
}

create_meta_features <- function(
    path, fbgn_path, gtf_path, present.symbols, mt.symbols, output_path) {
  cg_names = (
    path
    %>% paste0('/features.tsv.gz')
    %>% read.table(row.names=1)
    %>% rownames
  )
  fbgn_ann = read.table(fbgn_path, quote='', sep='\t', header=F)
  colnames(fbgn_ann)[c(1,3)] = c('symbol','flybase')
  fbgn_ann = as.data.frame(fbgn_ann)
  fbgn_ann_mapping = data.frame(
    annotation = fbgn_ann$V5,
    flybase = fbgn_ann$flybase,
    symbol = fbgn_ann$symbol
  )
  for (row in seq(nrow(fbgn_ann)))
    if (str_length(fbgn_ann[row,'V6'])) fbgn_ann_mapping = fbgn_ann_mapping %>% rbind(
    data.frame(
      annotation = strsplit(fbgn_ann[row, 'V6'], ',')[[1]],
      flybase = fbgn_ann[row, 'flybase'],
      symbol = fbgn_ann[row, 'symbol']
    )
  )
  feature_names = cg_names %>% map_feature_names
  meta.features = data.frame(
    row.names = feature_names,
    flybasecg = c('', str_replace(cg_names[-1], 'Dmel_', '')),
    flybase = c('', fbgn_ann_mapping$flybase[match(str_replace(cg_names[-1], 'Dmel_', ''), fbgn_ann_mapping$annotation)]),
    pass.min.cells = feature_names %in% present.symbols,
    is.mito = feature_names %in% mt.symbols
  )

  # Now build a length column, for FPKM calculation.
  gtf = read.table(
    gtf_path,
    sep = '\t',
    col.names = c('chr', 'source', 'type', 'start', 'end', 'sc', 'strand', 'fr', 'annotation'),
    header = F,
    quote = ''
  ) %>% subset(grepl('RNA', type)) # include "transcript" types in the reference
  gtf$flybase = gtf$annotation %>% str_extract(
    'gene_id "([^"]+)"',
    group=1
  )
  gtf = gtf %>% group_by(flybase) %>% summarise(
    chr = max(chr),
    start = min(start),
    end = max(end),
    strand = max(strand),
    transcript.length = max(abs(end - start)) + 1
  )
  meta.features = meta.features %>% left_join(gtf, 'flybase')
  rownames(meta.features) = feature_names

  write.csv(meta.features, output_path)
  output_path
}

load_flybase <- function(tenx_path, batch, features_keep, mt_features, meta_path, mt_pct=10, ribo_pct=40) {
  sce = Read10X(tenx_path)

  rownames(sce) = map_feature_names(rownames(sce))
  seur = CreateSeuratObject(sce)
  seur[['RNA']]@meta.features = read.csv(meta_path, row.names=1)
  seur$pct.mito = seur %>% PercentageFeatureSet(features = mt_features)
  seur$pct.ribo = seur %>% PercentageFeatureSet('^Rp[SL]')
  seur = seur[features_keep, seur$pct.mito < mt_pct & seur$pct.ribo < ribo_pct]
  seur$batch = batch
  seur = seur %>% SCTransform(
    vst.flavor = 'v2', do.correct.umi = F, return.only.var.genes = F,
    verbose = F
  )
  # Include H3-GFP as a scale.data gene, past the end of the original 3000 HVGs
  # so that it is easy to remove when integrating the data.
  VariableFeatures(seur) = VariableFeatures(seur) %>% union('H3-GFP')
  seur[['SCT']]@scale.data = seur[['SCT']]@scale.data %>% subset(
    rownames(.) %in% VariableFeatures(seur)
  )
  seur = seur %>% RunPCA(verbose = F)
}

call_nos.1 <- function(nos.1) {
  nos.1 = nos.1 %>% FindNeighbors(dims=1:8, verb=F) %>% FindClusters(res = 0.1, verb=F)
  DotPlot(nos.1, 'SCT', c('vas','bam', 'tj','lncRNA:roX2', 'Mst77F','Act57B', 'H3-GFP'))
  nos.1 = nos.1 %>% RenameIdents(`1`='germline', `2`='somatic', `4`='doublet',
                                 `3`='spermatocyte', `7`='muscle')
}

call_nos.2 <- function(nos.2) {
  # In nos.2, there is a cluster of spermatocytes (Mst84Da+) which is being
  # confused for GSC-like and which has very low nCount_RNA (75%ile = 795).
  # These cells were not integrated with the cells from other samples - they
  # were integrated into the "germline" cluster. These cells could even be a
  # spermatocyte-germline doublet cluster.
  # We are also going to call a mamo+ doublet cluster (likely another
  # germline-tumor - somatic doublet) which is very close to the somatic
  # cluster.
  nos.2 = nos.2 %>% FindNeighbors(dims=1:8, verb=F) %>% FindClusters(res = 0.15, verb=F)
  DotPlot(nos.2, 'SCT', c('vas','bam', 'tj','lncRNA:roX2','mamo', 'H3-GFP'))
  nos.2 = nos.2 %>% RenameIdents(
    # Cluster 3 is vas+, bam+, GFP+ but has lower nCount_RNA.
    `1`='germline', `3`='germline',
    `0`='somatic',
    `5`='doublet', `2`='doublet.spermatocyte', `4`='doublet.mamo+',
    `6`='spermatocyte', `7`='muscle')
}

call_tj.1 <- function(tj.1) {
  tj.1 = tj.1 %>% FindNeighbors(dims=1:8, verb=F) %>% FindClusters(res = 0.15, verb=F)
  tj.1 = tj.1 %>% RenameIdents(`1`='germline', `0`='somatic', `2`='doublet',
                               `3`='spermatocyte', `6`='muscle')
}

call_tj.2 <- function(tj.2) {
  tj.2 = tj.2 %>% FindNeighbors(dims=1:8, verb=F) %>% FindClusters(res = 0.1, verb=F)
  tj.2 = tj.2 %>% RenameIdents(`2`='germline', `0`='somatic', `1`='doublet',
                               `4`='spermatocyte', `5`='muscle')
}

gfp_pc = function(seurats) {
  features = Reduce(intersect, lapply(seurats, VariableFeatures)) %>% setdiff('H3-GFP')
  cells = do.call(
    cbind,
    lapply(
      seurats,
      \(seur) manova(
        genes ~ 0 + GFP,
        list(
          genes = seur[['SCT']]@scale.data[features, ] %>% t,
          GFP = seur[['SCT']]@scale.data['H3-GFP', ]
        )
      )$fitted.values %>% t
    )
  )
  irlba.cells = irlba(cells, length(seurats))
  rownames(irlba.cells$u) = features
  rownames(irlba.cells$v) = do.call(
    c,
    mapply(
      \(name, seurat) paste(name, Cells(seurat), sep='_'),
      names(seurats),
      seurats
    )
  )
  irlba.cells
}

filter_integrate_data = function(seurats) {
  # Remove doublets. Remove clusters that have not been renamed (numeric name).
  seurats = seurats %>% sapply(
    \(sce) sce %>% subset(
      cells = Cells(sce) %>% subset(
        !grepl('doublet', as.character(Idents(sce)))
        & !grepl('[0-9]', as.character(Idents(sce)))
      )
    ),
    simplify=F
  )
  gfp_dim_reduc = gfp_pc(seurats)
  seurats = seurats %>% sapply(
    \(sce) sce %>% subset(
      cells = Cells(sce) %>% subset(
        between(sce$nCount_RNA, 1300, 7500)
        & between(sce$nFeature_RNA, 500, 2500)
      )
    ),
    simplify=F
  )
  for (i in seq_along(seurats)) {
    Idents(seurats[[i]]) = Idents(seurats[[i]]) %>% fct_relabel(
      \(n) ifelse(grepl('[0-9]', n),
                  paste(seurats[[i]]$batch[1], n, sep='.'),
                  n))
    seurats[[i]] = seurats[[i]] %>% RenameCells(add.cell.id = seurats[[i]]$batch[1])
  }
  rna.key = Key(seurats[[1]]@assays[['RNA']])
  merge_RNA = CreateAssayObject(
    counts = do.call(
      cbind,
      seurats %>% sapply(\(sce) sce[['RNA']]@counts)
    )
  )
  meta.features = seurats[[1]][['RNA']]@meta.features
  Key(merge_RNA) = rna.key
  for (i in seq_along(seurats)) {
    # Integration features will be selected using RNA LogNormalize std variance.
    seurats[[i]] = seurats[[i]] %>% NormalizeData(assay='RNA', verb=F) %>% FindVariableFeatures(assay='RNA', verb=F)
    rna.key = Key(seurats[[i]]@assays[['RNA']])
    seurats[[i]]@assays[['RNA']] = (
      seurats[[i]]@assays[['RNA']]
      %>% SetAssayData(
        'counts',
        sparseMatrix(
          i=1,j=1,x=0,
          dims=dim(.),
          dimnames=dimnames(.)
        )
      ) %>% SetAssayData(
        'data',
        sparseMatrix(
          i=1,j=1,x=0,
          dims=dim(.),
          dimnames=dimnames(.)
        )
      )
    )
    Key(seurats[[i]]@assays[['RNA']]) = rna.key
  }

  integ.features = SelectIntegrationFeatures(object.list = seurats, nfeatures = 2001)
  # H3-GFP is not a high-quality feature for integration, as the batches differ
  # in their transgene and H3-GFP can be negatively correlated between batches.
  integ.features = integ.features %>% setdiff('H3-GFP')
  Upd_sc = (seurats
            %>% PrepSCTIntegration(anchor.features = integ.features)
            %>% FindIntegrationAnchors(normalization.method = 'SCT', anchor.features = integ.features)
            %>% IntegrateData(normalization.method = 'SCT'))

  # Pseudo-PC from power method of H3-GFP feature
  Upd_sc@misc$gfp_dim_reduc = gfp_dim_reduc
  gfp_pc1 = gfp_dim_reduc$u[,1]
  Upd_sc$somatic.score = (
    t(Upd_sc[['SCT']]@scale.data[names(gfp_pc1), ]) %*% gfp_pc1
  ) * -1

  # Dim reduction
  Upd_sc = Upd_sc %>% RunPCA(verb=F) %>% RunUMAP(dims=1:15)
  # PCA of germline and somatic clusters only
  Upd_subset = Upd_sc %>% subset(idents = c('germline','somatic'))
  # Now that we removed other clusters, regressing out pct.ribo is more
  # effective and is useful for producing interpretable principal components.
  Upd_subset[['integrated']]@scale.data = (
    manova(
      scale.data ~ pct.ribo,
      list(scale.data=t(Upd_subset[['integrated']]@scale.data),
           pct.ribo=Upd_subset$pct.ribo)
    )$residuals
    %>% scale %>% t
  )
  Upd_subset = Upd_subset %>% RunPCA(verb=F)
  # Create a pca object and re-key it
  pca.subset = Upd_sc[['pca']]
  Key(pca.subset) = 'pcasubset_'
  pca.subset@cell.embeddings[,] = NA
  pca.subset@cell.embeddings[colnames(Upd_subset),] = Upd_subset[['pca']]@cell.embeddings
  pca.subset@feature.loadings = Upd_subset[['pca']]@feature.loadings
  pca.subset@stdev = Upd_subset[['pca']]@stdev
  Upd_sc[['pca.subset']] = pca.subset

  # merge_RNA = merge_RNA %>% NormalizeData(verb=F)
  # Upd_sc@assays = list(RNA=merge_RNA)
  # DefaultAssay(Upd_sc) = 'RNA'

  # Include RNA, but kill SCT which is quite large
  Upd_sc = Upd_sc %>% SetAssayData(
    assay = 'SCT', slot = 'counts', Upd_sc[['RNA']]@counts
  ) %>% SetAssayData(
    assay = 'SCT', slot = 'data', Upd_sc[['RNA']]@counts
  ) %>% SetAssayData(
    assay = 'SCT',
    slot = 'scale.data',
    matrix(0, nrow=0, ncol=ncol(Upd_sc), dimnames=list(NULL, colnames(Upd_sc)))
  ) %>% SetAssayData(
    assay = 'RNA', slot = 'counts', merge_RNA@counts
  )
  Upd_sc[['RNA']]@meta.features = meta.features

  DefaultAssay(Upd_sc) = 'RNA'

  Idents(Upd_sc) = Idents(Upd_sc) %>% factor(c('germline','somatic','spermatocyte','muscle'))

  Upd_sc
}

celltype_Upd_sc = function(Upd_sc) {
  # integrated_assay = Upd_sc[['RNA']]
  # Key(integrated_assay) = 'integrated_'
  # Upd_sc[['integrated']] = integrated_assay
  Upd_sc = Upd_sc %>% FindNeighbors(verb = F)
  Upd_sc$celltype = (Upd_sc %>% FindClusters(graph.name = 'integrated_snn', res=0.5, verb=F))$seurat_clusters
  levels(Upd_sc$celltype) = c('germline','somatic',levels(Upd_sc$celltype))
  Upd_sc$celltype[Idents(Upd_sc) == 'germline'] <- 'germline'
  Upd_sc$celltype[Idents(Upd_sc) == 'somatic'] <- 'somatic'
  Upd_sc$celltype = Upd_sc$celltype %>% recode(
    # Mst87F+
    `7`='spermatocyte', `8`='spermatocyte',
    # Act57B+
    `9`='muscle',
    `0`='others', `1`='others', `2`='others', `3`='others', `4`='others',
    `5`='others', `6`='others', `10`='others', `11`='others'
  ) %>% fct_relevel('germline', 'somatic', 'spermatocyte', 'muscle', 'others')
  Upd_sc$celltype
}
