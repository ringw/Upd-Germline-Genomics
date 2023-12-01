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

load_flybase <- function(tenx_path, batch, features_keep, mt_features, mt_pct=10, ribo_pct=40) {
  sce = Read10X(tenx_path)

  rownames(sce) = map_feature_names(rownames(sce))
  seur = CreateSeuratObject(sce)
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
  DotPlot(nos.1, 'SCT', c('vas','bam', 'tj','lncRNA:roX2', 'H3-GFP'))
  nos.1 = nos.1 %>% RenameIdents(`1`='germline', `2`='somatic', `4`='doublet')
}

call_nos.2 <- function(nos.2) {
  nos.2 = nos.2 %>% FindNeighbors(dims=1:8, verb=F) %>% FindClusters(res = 0.15, verb=F)
  nos.2 = nos.2 %>% RenameIdents(`1`='germline', `0`='somatic', `5`='doublet')
}

call_tj.1 <- function(tj.1) {
  tj.1 = tj.1 %>% FindNeighbors(dims=1:8, verb=F) %>% FindClusters(res = 0.15, verb=F)
  tj.1 = tj.1 %>% RenameIdents(`1`='germline', `0`='somatic', `2`='doublet')
}

call_tj.2 <- function(tj.2) {
  tj.2 = tj.2 %>% FindNeighbors(dims=1:8, verb=F) %>% FindClusters(res = 0.1, verb=F)
  tj.2 = tj.2 %>% RenameIdents(`2`='germline', `0`='somatic', `1`='doublet')
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
  irlba.cells = irlba.cells %>% subset(names(.) %in% c('u','d'))
  rownames(irlba.cells$u) = features
  irlba.cells
}

filter_integrate_data = function(seurats) {
  seurats = seurats %>% sapply(
    \(sce) sce %>% subset(
      cells = Cells(sce) %>% subset(
        Idents(sce) != 'doublet'
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

  integ.features = SelectIntegrationFeatures(object.list = seurats, nfeatures = 2000)
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
  Upd_subset = Upd_subset %>% RunPCA(verb=F)
  Upd_sc[['pca.subset']] = Upd_sc[['pca']]
  Upd_sc[['pca.subset']]@cell.embeddings[,] = NA
  Upd_sc[['pca.subset']]@cell.embeddings[colnames(Upd_subset),] = Upd_subset[['pca']]@cell.embeddings
  Upd_sc[['pca.subset']]@feature.loadings = Upd_subset[['pca']]@feature.loadings
  Upd_sc[['pca.subset']]@stdev = Upd_subset[['pca']]@stdev

  merge_RNA = merge_RNA %>% NormalizeData(verb=F)
  Upd_sc@assays = list(RNA=merge_RNA)
  DefaultAssay(Upd_sc) = 'RNA'
  Idents(Upd_sc) = Idents(Upd_sc) %>% fct_relevel('germline', 'somatic')
  Upd_sc
}
