# Load all of the flybase annotation id tsv files. The first one in the
# character vector should contain most of the features (the release that we used
# -- FlyBase r6.47). The character vector can contain additional FlyBase
# references, to cover the FlyBase features that are now out of date, but
# nevertheless appeared in the GTF file that was built for genomic alignment
# (CR31927, CR32123, CR31485). Therefore, when we search for the FlyBase id,
# then the tsv file line coming from the first tsv file in the list will be
# expected to take priority in resolving the gene name.
load_flybase_annotation_mapping <- function(fbgn_paths) {
 apply_fbgn_ref <- function(fbgn_path) {
    fbgn_ann = read.table(fbgn_path, quote='', sep='\t', header=F)
    colnames(fbgn_ann)[c(1,3)] = c('symbol','flybase')
    fbgn_ann = as.data.frame(fbgn_ann)
    fbgn_ann_mapping = data.frame(
      annotation = fbgn_ann$V5,
      flybase = fbgn_ann$flybase,
      symbol = fbgn_ann$symbol
    )
    for (row in seq(nrow(fbgn_ann)))
      if (str_length(fbgn_ann[row,'V6']))
        fbgn_ann_mapping = fbgn_ann_mapping %>% rbind(
            data.frame(
              annotation = strsplit(fbgn_ann[row, 'V6'], ',')[[1]],
              flybase = fbgn_ann[row, 'flybase'],
              symbol = fbgn_ann[row, 'symbol']
            )
          )
    fbgn_ann_mapping
  }
  sapply(fbgn_paths, apply_fbgn_ref, simplify=FALSE) %>% do.call(rbind, .)
}

load_feature_names <- function(paths, fbgn_paths) {
  fbgn_ann_mapping <- load_flybase_annotation_mapping(fbgn_paths)
  annotation_names <- paths[1] %>%
    paste0('/features.tsv.gz') %>%
    read.table(row.names=1) %>%
    rownames
  sce.features <- c(
    "H3-GFP",
    fbgn_ann_mapping$symbol[
      match(
        annotation_names[-1] %>% str_replace("Dmel_", ""),
        fbgn_ann_mapping$annotation
      )
    ]
  )
  sce.features[duplicated(sce.features)] <- c(
    "H3-GFP",
    fbgn_ann_mapping$flybase[
      match(
        annotation_names[-1] %>% str_replace("Dmel_", ""),
        fbgn_ann_mapping$annotation
      )
    ]
  )[duplicated(sce.features)]
  stopifnot(all(!duplicated(sce.features)))
  sce.features
}

select_features <- function(paths, fbgn_paths) {
  feature_names <- paths %>% load_feature_names(fbgn_paths)
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

create_assay_data_sc <- function(
    path, sce.features, fbgn_paths, gtf_path, h3_gfp_gtf_path, present.symbols, output_path) {
  cg_names = (
    path
    %>% paste0('/features.tsv.gz')
    %>% read.table(row.names=1)
    %>% rownames
  )
 
  fbgn_ann_mapping <- load_flybase_annotation_mapping(fbgn_paths)
  meta.data = tibble(
    rowname = sce.features,
    flybasecg = c('', str_replace(cg_names[-1], 'Dmel_', '')),
    flybase = c(
      '',
      fbgn_ann_mapping$flybase[match(sce.features[-1], fbgn_ann_mapping$symbol)]
    ),
    pass.min.cells = sce.features %in% present.symbols,
    is.mito = grepl("^mt:", rowname)
  )

  # Now build a length column, for FPKM calculation.
  col.names <- c('chr', 'source', 'type', 'start', 'end', 'sc', 'strand', 'fr', 'annotation')
  # Include "transcript" types in the reference (RNA). What would be left out
  # are some rRNAs, which have a "pseudogene" annotation for the putative
  # transcript. Genes that have an "mRNA" may have multiple isoforms, but the
  # Cell Ranger workflow is including introns, so we can report the longest of
  # these isoforms from 5' of the first exon to 3' of the last exon. This is
  # generally the "CDS" which is another annotation that we get for these genes.
  gtf = read.table(
    gtf_path,
    sep = '\t',
    col.names = col.names,
    header = F,
    quote = ''
  ) %>% subset(grepl('RNA|pseudogene', type))
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
  meta.data = meta.data %>% left_join(gtf, 'flybase') %>% column_to_rownames
  meta.data["H3-GFP", "transcript.length"] <- with(
    read.table(h3_gfp_gtf_path, sep = "\t", col.names = col.names, header = F, quote = "") %>%
      subset(type == "mRNA") %>%
      head(1) %>%
      as.list,
    abs(end - start) + 1
  )

  write.csv(meta.data, output_path)
  output_path
}

read_seurat_sctransform <- function(
  tenx_path, batch, features_keep, assay_path, mt_pct=10, ribo_pct=40,
  return.only.var.genes = TRUE, run_pca = TRUE, sctransform = TRUE
) {
  sce = Read10X(tenx_path)

  assay.data <- read.csv(assay_path, row.names=1)
  stopifnot(
    all.equal(
      rownames(sce) %>% str_replace("UrbanChen.*|Dmel_", ""),
      assay.data$flybasecg
    )
  )
  rownames(sce) <- rownames(assay.data)
  seur = CreateSeuratObject(sce)
  # Every cell will have a nos.1_ prefix, etc.
  seur = seur %>% RenameCells(add.cell.id = batch)
  seur[['RNA']]@meta.data = read.csv(assay_path, row.names=1)
  mt_features <- rownames(seur[['RNA']]@meta.data) %>%
    subset(seur[['RNA']]@meta.data$is.mito)
  seur$pct.mito = seur %>% PercentageFeatureSet(features = mt_features)
  seur$pct.ribo = seur %>% PercentageFeatureSet('^Rp[SL]')
  seur = seur[features_keep, seur$pct.mito < mt_pct & seur$pct.ribo < ribo_pct]
  seur$batch = batch
  if (!sctransform) return(seur)
  seur = seur %>% SCTransform(
    vst.flavor = 'v2', do.correct.umi = F, return.only.var.genes = F,
    verbose = F, min_cells = 0
  )
  # Include H3-GFP as a scale.data gene, past the end of the original 3000 HVGs
  # so that it is easy to remove when integrating the data.
  VariableFeatures(seur) = VariableFeatures(seur) %>% union('H3-GFP')
  # Subset the residuals matrix now that we have selected our custom
  # VariableFeatures.
  if (return.only.var.genes)
    seur[['SCT']]@scale.data = seur[['SCT']]@scale.data %>% subset(
      rownames(.) %in% VariableFeatures(seur)
    )
  if (run_pca) seur <- seur %>% RunPCA(verbose = F)
  seur
}

call_nos.1 <- function(nos.1) {
  nos.1 = nos.1 %>% FindNeighbors(dims=1:8, verb=F) %>% FindClusters(res = 0.1, verb=F)
  DotPlot(nos.1, c('vas','bam', 'tj','lncRNA:roX2', 'Mst77F','Act57B', 'H3-GFP'), assay='SCT')
  nos.1 = nos.1 %>% RenameIdents(
    `1`='germline', `2`='somatic', `3`='doublet',
    `4`='spermatocyte', `8`='muscle',

    `5`='somaticprecursor',
    `6`='somaticprecursor'
  )
}

call_nos.2 <- function(nos.2) {
  # Doublets in Nos-Rep2 the most difficult to tease apart (including doublets
  # of the differentiated cell types). We applied a greater clustering res,
  # and then recombined somatic and germline cells that actually each form a
  # large (e.g. blob, elliptical) population on a nos.2 UMAP. Then, the
  # "doublet.spermatocyte" cells, when left in the Upd_sc object, are very
  # spread-out germ cells, some of which are adjacent to either the germline
  # ident or the somatic ident, and should be removed as they do have some %
  # expression of vas and tj in this cluster. (along with medium % w-cup expr).
  nos.2 = nos.2 %>% FindNeighbors(dims=1:9, verb=F) %>% FindClusters(res = 0.5, verb=F)
  # mamo, wb may be somaticprecursor markers
  DotPlot(nos.2, c('vas','bam', 'tj','lncRNA:roX2','mamo','wb', 'H3-GFP', 'w-cup', 'Act57B'))
  nos.2 = nos.2 %>% RenameIdents(
    `0`='somatic',
    `4`='somatic',
    `5`='somatic',

    `1`='germline',
    `2`='germline',
    `6`='germline',

    `11`='spermatocyte',
    `13`='muscle',

    `8`='doublet.stemlike',
    # Cluster "9" has quantile(vas[assay=SCT], 0.90) > 1 so it looks like it is
    # a differentiated cluster significantly contaminated by vasa (stem cells),
    # may be doublet cells.
    `9`='doublet.spermatocyte',

    `12`='somaticprecursor'
  )
}

call_tj.1 <- function(tj.1) {
  tj.1 = tj.1 %>% FindNeighbors(dims=1:8, verb=F) %>% FindClusters(res = 0.15, verb=F)
  tj.1 = tj.1 %>% RenameIdents(
    `1`='germline', `0`='somatic', `2`='doublet',
    `3`='spermatocyte', `6`='muscle',
    `4`='somaticprecursor'
  )
}

call_tj.2 <- function(tj.2) {
  tj.2 = tj.2 %>% FindNeighbors(dims=1:8, verb=F) %>% FindClusters(res = 0.1, verb=F)
  tj.2 = tj.2 %>% RenameIdents(
    `2`='germline', `0`='somatic', `1`='doublet',
    `4`='spermatocyte', `5`='muscle',
    `3`='somaticprecursor'
  )
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
  assay_data = seurats[[1]][["RNA"]]@meta.data
  # Remove doublets. Use relabeling of the numeric (unidentified) clusters so
  # that we learn a cluster mean independently for each gene and for each
  # unidentified cluster.
  seurats = seurats %>% sapply(
    \(sce) sce %>%
      subset(
        cells = Cells(sce) %>% subset(
          !grepl('doublet', as.character(Idents(sce)))
        )
      ) %>%
      `Idents<-`(
        value = Idents(.) %>%
          fct_relabel(
            \(n) ifelse(grepl('[0-9]', n),
                        paste(sce$batch[1], n, sep='.'),
                        n))
      ),
    simplify=F
  )
  gfp_dim_reduc = gfp_pc(seurats)
  seurats = seurats %>% sapply(
    \(sce) sce %>% subset(
      cells = Cells(sce) %>% subset(
        between(sce$nCount_RNA, 2200, 7500)
        & between(sce$nFeature_RNA, 500, 2500)
      )
    ),
    simplify=F
  )

  integ.features = SelectIntegrationFeatures(object.list = seurats, nfeatures = 2001)
  # H3-GFP is not a high-quality feature for integration, as the batches differ
  # in their transgene and H3-GFP can be negatively correlated between batches.
  integ.features = integ.features %>% setdiff('H3-GFP')
  # Now need scale.data to contain the integration features.
  seurats <- seurats %>%
    sapply(
      \(sce) sce %>%
        SCTransform(
          vst.flavor = 'v2', do.correct.umi = F,
          residual.features = integ.features,
          verbose = F, min_cells = 0
        )
    )
  Upd_sc = (seurats
            %>% PrepSCTIntegration(anchor.features = integ.features)
            %>% FindIntegrationAnchors(normalization.method = 'SCT', anchor.features = integ.features)
            %>% IntegrateData(normalization.method = 'SCT'))

  # Pseudo-PC from power method of H3-GFP feature
  Upd_sc@misc$gfp_dim_reduc = gfp_dim_reduc
  gfp_pc1 = gfp_dim_reduc$u[,1]
  Upd_sc$somatic.score = (
    t(Upd_sc[['SCT']]@scale.data[names(gfp_pc1), ]) %*% gfp_pc1
  )

  # Dim reduction
  Upd_sc = Upd_sc %>% RunPCA(verb=F) %>% RunUMAP(dims=1:15, seed.use=2)
  # Update cell embedding so that GSC is on left. We would also like
  # differentiated germline clusters (Mst77F+) in the +y direction while the
  # largest contaminant cluster (muscle, Act57B+) is in the -y.
  Upd_sc[['umap']]@cell.embeddings[,1] = (
    Upd_sc[['umap']]@cell.embeddings[,1]
    * sign(
      coef(
        lm(
          identsomatic ~ umap_1,
          FetchData(Upd_sc, c("ident", "umap_1")) %>%
            mutate(identsomatic = ident == "somatic")
        )
      )['umap_1']
    )
  )
  # We noticed regardless of random seed that y needs to be flipped.
  Upd_sc[['umap']]@cell.embeddings[,2] <- -Upd_sc[['umap']]@cell.embeddings[,2]

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
  DefaultAssay(Upd_subset) <- "integrated"
  Upd_subset <- Upd_subset %>% RunPCA(verb=F)
  # Create a pca object and re-key it
  pcasubset = Upd_sc[['pca']]
  pcasubset.key = 'pcasubset_'
  cell.embeddings <- pcasubset@cell.embeddings
  cell.embeddings[,] <- NA
  cell.embeddings[colnames(Upd_subset),] <- Upd_subset[["pca"]]@cell.embeddings
  colnames(cell.embeddings) <- paste0(pcasubset.key, 1:ncol(cell.embeddings))
  feature.loadings <- Upd_subset[["pca"]]@feature.loadings
  colnames(feature.loadings) <- colnames(cell.embeddings)
  pcasubset <- CreateDimReducObject(
    cell.embeddings,
    feature.loadings,
    stdev = Upd_subset[["pca"]]@stdev,
    key = pcasubset.key,
    assay = "integrated"
  )
  Upd_sc[['pcasubset']] = pcasubset

  # Include RNA, but kill SCT which is quite large. We will use RNA assay for
  # feature plots, and integrated assay for relating the transcriptome between
  # different cells, but no longer need the individual runs of SCT (all genes).
  Upd_sc = Upd_sc %>% SetAssayData(
    assay = 'SCT',
    layer = 'counts',
    sparseMatrix(i=1, j=1, x=0, dims=dim(Upd_sc[["SCT"]]), dimnames=dimnames(Upd_sc[["SCT"]]))
  ) %>% SetAssayData(
    assay = 'SCT',
    layer = 'data',
    sparseMatrix(i=1, j=1, x=0, dims=dim(Upd_sc[["SCT"]]), dimnames=dimnames(Upd_sc[["SCT"]]))
  ) %>% SetAssayData(
    assay = 'SCT',
    layer = 'scale.data',
    matrix(0, nrow=0, ncol=ncol(Upd_sc), dimnames=list(NULL, Cells(Upd_sc)))
  )
  # Join RNA matrices (counts) together
  Upd_sc[['RNA']] <- Upd_sc[['RNA']] %>% JoinLayers
  Upd_sc[['RNA']]@meta.data = assay_data

  Upd_sc$indep_idents = Idents(Upd_sc)
  DefaultAssay(Upd_sc) <- "integrated"
  Upd_sc <- Upd_sc %>%
    FindNeighbors(dims = 1:7, nn.method = "rann") %>%
    FindClusters(res = 0.1, random.seed = 0)
  DefaultAssay(Upd_sc) <- "RNA"
  Idents(Upd_sc) <- Idents(Upd_sc) %>%
    fct_recode(
      germline="0",
      somatic="1",
      germline="2",
      spermatocyte="3",
      somaticprecursor="4",
      muscle="5"
    ) %>%
    fct_relevel(
      c("germline", "somatic", "spermatocyte", "somaticprecursor", "muscle")
    )

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

run_umap_on_batch <- function(seurat, metadata) {
  metadata <- metadata %>% read.csv(row.names = 1)
  seurat <- seurat %>%
    subset(
      cells = Cells(seurat) %>% subset(
        between(seurat$nCount_RNA, 2200, 7500)
        & between(seurat$nFeature_RNA, 500, 2500)
      )
    )
  Idents(seurat) <- metadata[Cells(seurat), "ident"] %>%
    replace_na("doublet") %>%
    factor(c(sce.clusters$cluster, "doublet"))
  seurat <- seurat %>% RunUMAP(dims = 1:25)
  DefaultAssay(seurat) <- "RNA"
  seurat <- seurat %>% NormalizeData
  cols_to_fetch <- c("ident", "umap_1", "umap_2", "tj", "vas", "H3-GFP")
  data <- FetchData(seurat, cols_to_fetch)
  # Flip UMAP. We want germline in -x direction and muscle in -y direction.
  if (order(data %>% subset(ident == "germline", select=c(umap_1, umap_2)) %>% summarise_all(~ abs(mean(.))) %>% unlist)[1] != 2) {
    seurat[["umap"]]@cell.embeddings <- seurat[["umap"]]@cell.embeddings %*% (
      matrix(
        c(0, 1, 1, 0),
        nrow = 2,
        dimnames = rep(list(colnames(seurat[["umap"]]@cell.embeddings)), 2)
      )
    )
  }
  data <- FetchData(seurat, cols_to_fetch)
  seurat[["umap"]]@cell.embeddings <- seurat[["umap"]]@cell.embeddings %*% (
    matrix(
      c(
        data %>% subset(ident == "germline") %>% with(-umap_1) %>% mean %>% sign,
        0,
        0,
        data %>% subset(ident == "muscle") %>% with(-umap_2) %>% mean %>% sign
      ),
      nrow = 2,
      dimnames = rep(list(colnames(seurat[["umap"]]@cell.embeddings)), 2)
    )
  )
  FetchData(seurat, cols_to_fetch)
}

unintegrated_report_cluster_expression <- function(
  seurat_obj, doublet_names, gene_list, Upd_metadata
) {
  idents <- Upd_metadata[Cells(seurat_obj), "ident"] %>%
    replace(is.na(.), "excluded") %>%
    ordered(c(sce.clusters$cluster, "excluded"))
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- seurat_obj %>% NormalizeData(verb=F)
  cluster_data <- seurat_obj %>%
    FetchData(gene_list) %>%
    split(seurat_obj$seurat_clusters) %>%
    sapply(\(m) colMeans(m != 0) %>% round(5)) %>%
    t %>%
    as.data.frame %>%
    dplyr::rename(roX2="lncRNA:roX2")
  cluster_data <- FetchData(seurat_obj, "nCount_RNA")[, 1] %>%
      split(seurat_obj$seurat_clusters) %>%
      sapply(\(v) quantile(v, c(0.25, 0.5, 0.75))) %>%
      matrix(
        ncol = 3,
        byrow = TRUE,
        dimnames = list(
          levels(seurat_obj$seurat_clusters),
          c("nUMI_IQR-", "nUMIMedian", "nUMI_IQR+")
        )
      ) %>%
      as.data.frame %>%
      cbind(
        tibble(
          n = as.data.frame(table(seurat_obj$seurat_clusters))$Freq,
          removed_doublet = levels(seurat_obj$seurat_clusters) %in% doublet_names,
          is_high_nUMI = .$`nUMI_IQR-` >= (
            seurat_obj$nCount_RNA %>% subset(. >= 2200) %>% median
          ) * 0.9,
          is_doublet_scDblFinder = SingleCellExperiment(list(counts=GetAssayData(seurat_obj, "RNA", "counts"))) %>%
            `colLabels<-`(value = seurat_obj$seurat_clusters) %>%
            findDoubletClusters %>%
            as.data.frame %>%
            rownames_to_column %>%
            pull(num.de, rowname) %>%
            `[`(levels(seurat_obj$seurat_clusters)) %>%
            isOutlier(type="lower", log=TRUE) %>%
            sapply(isTRUE),
          mode = idents %>%
            split(seurat_obj$seurat_clusters) %>%
            mapply(
              \(chrs, filter_cells) c(
                with(
                  as.data.frame(table(chrs)),
                  chrs[which.max(Freq)]
                ),
                with(
                  as.data.frame(table(chrs[filter_cells])),
                  if (max(Freq) > 25) Var1[which.max(Freq)] else tail(Var1, 1)
                )
              ) %>%
                droplevels %>%
                levels %>%
                head(1),
              .,
              with(
                FetchData(seurat_obj, c("nCount_RNA", "nFeature_RNA")),
                between(nCount_RNA, 2200, 7500)
                & between(nFeature_RNA, 500, 2500)
              ) %>%
                split(seurat_obj$seurat_clusters)
            )
        ),
        cluster_data,
        .
      ) %>%
    rownames_to_column("Seurat Cluster ID")
}