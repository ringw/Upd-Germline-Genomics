# Differentially expressed genes.
analyze_sce_to_csv <- function(seurats, output_path) {
  for (name in names(seurats)) seurats[[name]] = seurats[[name]] %>% RenameCells(
    add.cell.id = name
  )
  for (i in seq_along(seurats)) {
    Idents(seurats[[i]]) = Idents(seurats[[i]]) %>% fct_relabel(
      \(n) ifelse(grepl('[0-9]', n),
                  paste(seurats[[i]]$batch[1], n, sep='.'),
                  n))
    seurats[[i]]@meta.data = seurats[[i]]@meta.data %>% subset(
      select = -grep('SCT_snn_res', colnames(.))
    )
  }
  sce = SingleCellExperiment(
    list(counts = do.call(cbind, seurats %>% sapply(\(seurat) seurat[['RNA']]@counts, simplify=F))),
    colData = do.call(rbind, seurats %>% setNames(NULL) %>% sapply(\(seurat) seurat@meta.data, simplify=F))
  )
  colData(sce) = cbind(
    ident = do.call(
      c,
      sapply(seurats, Idents, simplify=F)
    ) %>% setNames(colnames(sce)),
    colData(sce)
  )

  quickClusters = quickCluster(sce, min.size = 1000)
  # Choose the large germline-like cluster as reference cluster with an average
  # size factor of 1, as this is the cluster that we are going to use to
  # establish the baseMean quantification of the gene.
  ref.clust = assay(sce)['vas',] %>% split(quickClusters) %>% sapply(sum) %>% which.max
  colData(sce)$size_factor = pooledSizeFactors(
    sce,
    clusters = quickClusters,
    ref.clust = ref.clust
  )
  write.csv(colData(sce), file=output_path)
  output_path
}

fit_glm <- function(seurats, metadata_file) {
  metadata = read.csv(metadata_file, row.names=1)
  for (name in names(seurats)) seurats[[name]] = seurats[[name]] %>% RenameCells(
    add.cell.id = name
  )
  for (i in seq_along(seurats)) {
    Idents(seurats[[i]]) = Idents(seurats[[i]]) %>% fct_relabel(
      \(n) ifelse(grepl('[0-9]', n),
                  paste(seurats[[i]]$batch[1], n, sep='.'),
                  n))
    seurats[[i]] = seurats[[i]] %>% subset(
      idents = setdiff(levels(Idents(.)), 'doublet')
    )
  }
  counts = do.call(cbind, seurats %>% sapply(\(seurat) seurat[['RNA']]@counts, simplify=F))
  Upd_glm <- glm_gp(
    counts,
    ~ ident + batch,
    metadata[, colnames(counts)],
    size_factors = metadata$size_factor,
    on_disk = F,
    verbose = T
  )
  assay(Upd_glm$data) = assay(Upd_glm$data) %>% as('sparseMatrix')
  Upd_glm
}