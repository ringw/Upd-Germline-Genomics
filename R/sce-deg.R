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

fit_glm <- function(tenx_paths, metadata_file) {
  metadata = read.csv(metadata_file, row.names=1) %>% as.data.frame
  counts = mapply(
    \(name, path) {
      mat = Read10X(path)
      colnames(mat) = paste(name, colnames(mat), sep='_')
      mat
    },
    names(tenx_paths),
    tenx_paths,
    SIMPLIFY=F,
    USE.NAMES=F
  )
  counts = do.call(cbind, counts)
  rownames(counts) = rownames(counts) %>% map_feature_names
  # Remove counts that are missing from metadata, and remove doublets
  counts.keep = colnames(counts) %>% subset(
    metadata[., 'ident'] %in% c(
      'germline','somatic','spermatocyte','muscle'
    )
  )
  gene.names = rownames(counts)
  counts = counts[, counts.keep]
  metadata = metadata[counts.keep, ]
  metadata$ident = metadata$ident %>% fct_relevel(
    'germline','somatic','spermatocyte','muscle'
  )
  counts = writeHDF5Array(counts, with.dimnames = T)
  Upd_glm <- glm_gp(
    counts,
    ~ ident + batch,
    metadata,
    size_factors = metadata[, 'size_factor'],
    subsample = 10000,
    on_disk = T,
    verbose = T
  )
  # Move data from hdf5 to memory:
  # Store data as sparseMatrix.
  assay(Upd_glm$data) = assay(Upd_glm$data) %>% as('sparseMatrix')
  # We did not specify an offset per gene. Grab the offset per cell (column).
  # We can reconstruct a full-size Offset matrix later.
  Upd_glm$Offset = as.matrix(Upd_glm$Offset[1,, drop=F])
  # Mu is quite large and uncompressible, and can be recalculated using
  # glmGamPoi:::calculate_mu.
  Upd_glm$Mu = NULL
  Upd_glm
}

create_fpkm_reference_longest_isoform <- function(
    meta_path,
    Upd_glm,
    contrast,
    output_path) {
  col_data = colData(Upd_glm$data)[
    colAlls(t(Upd_glm$model_matrix) == contrast),
  ]
  # Predict for the sum of size factors, divide by the sum of counts. We will
  # divide by base pairs reference (bp), so the scale factor is 1B, not 1MM.
  offset_predictions = (
    log(sum(col_data$size_factor))
    - log(sum(col_data$nCount_RNA))
    + log(1000)*3
  )
  predictions = predict(
    Upd_glm,
    matrix(contrast, nrow=1),
    type='response',
    offset=offset_predictions
  )

  features = read.csv(meta_path, row.names=1) %>% rownames_to_column %>% subset(
    !duplicated(flybase) & !is.na(flybase)
  )
  features = features %>% inner_join(
    data.frame(
      rowname = rownames(predictions),
      fragments_per_kilo_mm = predictions
    ),
    by='rowname'
  )

  # Include introns mode - use entire length of mRNA (first exon to last exon)
  features = features %>% arrange(desc(transcript.length))
  features$fpkm = features$fragments_per_kilo_mm / (abs(features$end - features$start) + 1)

  bed = features %>% subset(
    select=c(chr,start,end,flybase,fpkm,strand)
  )
  bed = bed %>% arrange(desc(fpkm))
  write.table(bed, output_path, row.names=F, col.names=F, quote=F, sep='\t')
  output_path
}

apeglm_coef_table <- function(g, gene_list, coef=2, num_cells=20000) {
  message('Construct dense matrices in memory for regression')
  g$Offset = matrix(
    g$Offset[1,, drop=T],
    byrow=T,
    nrow=nrow(g$data),
    ncol=ncol(g$data),
    dimnames=dimnames(g$data)
  )
  g$Mu = glmGamPoi:::calculate_mu(
    Beta=g$Beta,
    model_matrix=g$model_matrix,
    offset_matrix=g$Offset
  )
  message('F-test of regression without shrinkage')
  mle = test_de(g, colnames(g$Beta)[coef])
  message('Calculate regression s.e. using QR decomposition')
  lfc = predict(
    g,
    matrix(
      seq(ncol(g$Beta)) == coef,
      nrow=1
    ),
    offset=0,
    se.fit=T
  ) %>% with(cbind(fit, se.fit))
  message('Construct a smaller Offset matrix for Approximate Posterior')
  cell_logical = with_seed(
    1,
    seq(ncol(g$data)) %in% sample(ncol(g$data), num_cells)
  )
  g$Mu = NULL
  Offset = matrix(
    g$Offset[1, cell_logical, drop=T],
    byrow=T,
    nrow=length(gene_list),
    ncol=num_cells,
    dimnames=list(gene_list, colnames(g$data)[cell_logical])
  )
  g$Offset = NULL
  message('Optimize Approximate Posterior')
  lfc = apeglm(
    assay(g$data)[gene_list, cell_logical],
    g$model_matrix[cell_logical,],
    log.lik=NULL,
    param=g$overdispersion_shrinkage_list$dispersion_trend[
      match(gene_list, rownames(g$data))
    ],
    coef=coef,
    mle=lfc,
    method='nbinomCR',
    offset=Offset
  )
  lfc$mle.test = mle
  lfc
}

# coef: ident - germline
# Contrast: Subtract 0.5*somatic (take geometric mean of germline and somatic
# mean as the denominator)
apeglm_contrasts = list(
  spermatocyte=c(0,-0.5,1,0,0,0,0),
  muscle=c(0,-0.5,0,1,0,0,0)
)

apeglm_coef_table_contrast <- function(g, gene_list, contrast, num_cells=20000) {
  message('Construct dense matrices in memory for regression')
  g$Offset = matrix(
    g$Offset[1,, drop=T],
    byrow=T,
    nrow=nrow(g$data),
    ncol=ncol(g$data),
    dimnames=dimnames(g$data)
  )
  g$Mu = glmGamPoi:::calculate_mu(
    Beta=g$Beta,
    model_matrix=g$model_matrix,
    offset_matrix=g$Offset
  )
  message('F-test of regression without shrinkage')
  mle = test_de(g, contrast)
  mle = NULL
  message('Calculate regression s.e. using QR decomposition')
  lfc = predict(
    g,
    matrix(
      contrast,
      nrow=1
    ),
    offset=0,
    se.fit=T
  ) %>% with(cbind(fit, se.fit))

  message('Construct a smaller Offset matrix for Approximate Posterior')
  cell_logical = with_seed(
    1,
    seq(ncol(g$data)) %in% sample(ncol(g$data), num_cells)
  )
  g$Mu = NULL
  Offset = matrix(
    g$Offset[1, cell_logical, drop=T],
    byrow=T,
    nrow=length(gene_list),
    ncol=num_cells,
    dimnames=list(gene_list, colnames(g$data)[cell_logical])
  )
  g$Offset = NULL

  message('Construct a transformed design matrix')
  # Y = beta X^T (beta coefficients in row vectors)
  # beta -> beta * betamat.transform (one column takes the cluster estimate
  #                                   coefficient and transforms it to be the
  #                                   desired contrast)
  # X^T -> solve(betamat.transform) * X^T
  # X -> X * solve(betamat.transform^T)
  betamat.transform = diag(nrow=ncol(g$Beta))
  dimnames(betamat.transform) = list(
    colnames(g$Beta)
  ) %>% rep(2)
  betamat.transform[,2] = contrast
  colnames(betamat.transform)[2] = 'contrast'
  model_matrix = g$model_matrix %*% solve(t(betamat.transform))

  message('Optimize Approximate Posterior')
  lfc = apeglm(
    assay(g$data)[gene_list, cell_logical],
    model_matrix[cell_logical,],
    log.lik=NULL,
    param=g$overdispersion_shrinkage_list$dispersion_trend[
      match(gene_list, rownames(g$data))
    ],
    coef=2,
    mle=lfc,
    method='nbinomCR',
    offset=Offset
  )
  lfc$mle.test = mle
  lfc
}

download_bam_transcripts <- function(download_script, metadata, sample, cluster, output_path) {
  dir.create(dirname(output_path), showW = FALSE)
  with_options(
    list(warn=2),
    system2(
      c(
        download_script,
        metadata,
        sample,
        cluster,
        output_path
      )
    )
  )
  output_path
}

pseudobulk_cpm <- function(
  tx_files, size_factors, gtf_file
) {
  cluster_name = factor(
    tx_files$tx_file %>% str_extract('([a-z]+)\\.txt', group=1),
    sce.clusters$cluster
  )
  if (all(!is.na(cluster_name)))
  pseudobulk_colData = cbind(size_factors, ident = cluster_name)
  else
    pseudobulk_colData = cbind(size_factors, ident = 'intercept')
  # Use nos.2 as reference level for FPKM calculation
  pseudobulk_colData$batch = pseudobulk_colData$batch %>%
    factor(
      c('nos.2', 'nos.1', 'tj.1', 'tj.2')
    )
  pseudobulk_counts <- tx_files$tx_file %>%
    sapply(\(f) read.table(f, col.names=c(f, 'tx')), simplify=F) %>%
    purrr::reduce(\(a,b) full_join(a,b, 'tx')) %>%
    column_to_rownames('tx') %>%
    as.matrix %>%
    replace(is.na(.), 0)
  colnames(pseudobulk_counts) = rownames(pseudobulk_colData)
  if (length(unique(pseudobulk_colData$batch)) > 1 && length(unique(pseudobulk_colData$ident)) > 1)
  dds <- DESeqDataSetFromMatrix(
    pseudobulk_counts, pseudobulk_colData, ~ 0 + ident + batch)
  else
    dds <- DESeqDataSetFromMatrix(pseudobulk_counts, pseudobulk_colData, ~ 1)
  sizeFactors(dds) <- pseudobulk_colData$size_factor
  dds <- dds %>% DESeq(fitType = 'local')

  gtf <- read.table(
    gtf_file,
    sep = '\t',
    col.names = c('chr', 'source', 'type', 'start', 'end', 'sc', 'strand', 'fr', 'annotation'),
    header = F,
    quote = ''
  ) %>% subset(grepl('RNA', type)) # include "transcript" types in the reference
  gtf$tx <- gtf$annotation %>% str_extract(
    'transcript_id "([^"]+)"',
    group = 1
  )
  gtf$gene_id <- gtf$annotation %>% str_extract(
    'gene_id "([^"]+)"',
    group = 1
  )

  cpm_table <- data.frame(
    tx = gtf$tx,
    gene_id = gtf$gene_id,
    length = abs(gtf$end - gtf$start + 1)
  )
  if (length(resultsNames(dds)) > 1)
  for (cluster_name in sce.clusters$cluster) {
    res <- results(dds, name = paste0("ident", cluster_name))
    cpm_table <- cpm_table %>%
      left_join(
        data.frame(
          tx = rownames(res),
          clust = exp(res$log2FoldChange * log(2) + log(1000) * 2)
        ),
        "tx"
      )
    colnames(cpm_table)[ncol(cpm_table)] <- cluster_name
    }
  else {
    cpm_table <- cpm_table %>%
      left_join(
        data.frame(
          tx = rownames(results(dds)),
          clust = exp(results(dds)$log2FoldChange * log(2) + log(1000) * 2)
        ),
        "tx"
      )
    colnames(cpm_table)[4] = "cpm"
  }
  cpm_table %>%
    column_to_rownames("tx") %>% replace(is.na(.), 0)
}

tx_cpm_to_fpkm <- function(cpm_table, metafeatures_path) {
  metafeatures <- read.csv(metafeatures_path, row.names = 1)
  fpkm_table <- cpm_table %>%
    subset(select=-c(gene_id)) %>%
    as.matrix %>%
    `*`(1000 / .[, "length"]) %>%
    subset(select=-length) %>%
    cbind(cpm_table %>% subset(select=gene_id), .)
  fpkm_table <- fpkm_table %>%
    group_by(gene_id) %>%
    summarise_all(max)
  colnames(fpkm_table)[1] <- "flybase"
  fpkm_table <- metafeatures %>%
    rownames_to_column %>%
    subset(select=c(rowname,flybase)) %>%
    left_join(fpkm_table, "flybase") %>%
    column_to_rownames %>%
    subset(select = -flybase) %>%
    as.matrix
  fpkm_table
}

longest_isoform_fpkm <- function(
    Upd_glm, meta_path) {
  fpkm_table <- matrix(
    nrow = nrow(Upd_glm$Beta),
    ncol = nrow(sce.clusters),
    dimnames = list(rownames(Upd_glm$Beta), sce.clusters$cluster)
  )
  for (i in seq(nrow(sce.clusters))) {
    contrast <- sce.clusters[i, "contrast"] %>% unlist
    col_data = colData(Upd_glm$data)[
      colAlls(t(Upd_glm$model_matrix) == contrast),
    ]
    # Predict for the sum of size factors, divide by the sum of counts. We will
    # divide by base pairs reference (bp), so the scale factor is 1B, not 1MM.
    offset_predictions = (
      log(sum(col_data$size_factor))
      - log(sum(col_data$nCount_RNA))
      + log(1000)*3
    )
    fpkm_table[, i] <- predict(
      Upd_glm,
      matrix(contrast, nrow=1),
      type='response',
      offset=offset_predictions
    )
  }

  metafeatures <- read.csv(meta_path, row.names=1)
  features = metafeatures %>% rownames_to_column %>% subset(
    !duplicated(flybase) & !is.na(flybase)
  )
  fpkm_table <- fpkm_table %>%
    as.data.frame %>%
    rownames_to_column %>%
    left_join(data.frame(rowname = features$rowname), ., "rowname") %>%
    column_to_rownames %>%
    as.matrix
  fpkm_table <- fpkm_table / features$transcript.length
  fpkm_table %>%
    as.data.frame %>%
    rownames_to_column %>%
    right_join(data.frame(rowname=rownames(metafeatures)), "rowname") %>%
    column_to_rownames %>%
    as.matrix
}

reference_sort_by_fpkm_table <- function(
    fpkm_table,
    cluster_name,
    meta_path,
    output_path) {
  fpkm_table <- fpkm_table %>% subset(rowAlls(is.finite(fpkm_table)))

  features = read.csv(meta_path, row.names=1) %>% rownames_to_column %>% subset(
    !duplicated(flybase) & !is.na(flybase)
  ) %>%
    inner_join(fpkm_table %>% as.data.frame %>% rownames_to_column, "rowname")
  colnames(features) <- colnames(features) %>%
    replace(. == cluster_name, "fpkm")

  bed = features %>% subset(
    fpkm > 0,
    select=c(chr,start,end,flybase,fpkm,strand)
  )
  bed = bed %>% arrange(desc(fpkm))
  write.table(bed, output_path, row.names=F, col.names=F, quote=F, sep='\t')
  output_path
}

# Read the bed file (no col names) as produced by reference_sort_by_fpkm_table.
# Note that we are only going to analyze the chromosomes given by
# names(chr.lengths).
load_flybase_bed <- function(bed_path) {
  read.table(
    bed_path, col.names = c("chr","start","end","flybase","fpkm","strand"),
    quote = "",
    sep = "\t"
  ) %>%
    subset(chr %in% names(chr.lengths))
}

dot_plot_fpkm <- function(genotypes, gene_list, logfpkm_min=0, logfpkm_max=3, oob_squish=FALSE) {
  data <- genotypes %>%
    sapply(\(df) df[gene_list, ] %>% rownames_to_column("gene"), simplify = FALSE) %>%
    bind_rows(.id = "genotype")
  # Show gene in order of appearance in gene_list
  data$gene <- data$gene %>% factor(levels = unique(.)) %>%
    fct_relabel(\(v) v %>% str_replace('lncRNA:', '') %>% str_replace('Hsromega', 'Hsr\u03C9'))
  data$genotype <- data$genotype %>% factor(levels = names(genotypes))
  g <- data %>% ggplot(
    aes(x = genotype, y = gene, size = percent_expressed, color = log(fpkm)/log(10))
  ) + geom_point()
  if (oob_squish)
    g <- g + scale_color_gradient(
      limits = c(logfpkm_min, logfpkm_max),
      oob = squish
    )
  else
    g <- g + scale_color_gradient(
      limits = c(logfpkm_min, logfpkm_max)
    )
  g + scale_y_discrete(limits=rev) + scale_size_continuous(
    labels=percent, limits=c(0.01,0.99), range=c(1, 10)
  ) + labs(
    x = NULL, y = NULL, color = bquote(log[10]*"(FPKM)"), size = "expression"
  ) + theme_cowplot()
}

# Another cluster-gene quantification method: Take the quantile of interest for
# the scaled data.
quantify_quantiles <- function(
  seurat, metadata, assay = "SCT", slot = "scale.data",
  per_million = FALSE
) {
  metadata <- metadata %>% read.csv(row.names = 1) %>%
    rownames_to_column %>%
    right_join(data.frame(rowname = colnames(seurat)), by="rowname")
  clusters_of_interest <- c("germline", "somatic")
  cells <- clusters_of_interest %>%
    sapply(
      \(n) GetAssayData(
        seurat[[assay]],
        slot
      )[
        ,
        as.character(metadata$ident) == n
        & metadata$nCount_RNA_filter == "nCount_RNA_pass"
      ],
      simplify = FALSE
    )
  qs <- sapply(
    cells,
    \(m) rowQuantiles(m, p=seq(1, 100) / 100, na.rm = TRUE),
    simplify = FALSE
  )
  if (per_million) {
    qs <- qs %>%
      sapply(
        \(mat) {
          scale_multiplier <- diag(
            1000 * 1000 / colSums(mat)
          )
          rownames(scale_multiplier) = colnames(mat)
          colnames(scale_multiplier) = colnames(mat)
          mat %*% scale_multiplier
        },
        simplify = FALSE
      )
  }
  qs
}

combine_gene_quantiles <- function(nested_list) {
  clusters_of_interest <- c("germline", "somatic")
  sapply(
    clusters_of_interest,
    \(n) purrr::reduce(
      sapply(nested_list, \(l) l[[n]], simplify=F),
      `+`
    ) / length(nested_list),
    simplify=F
  )
}

venn_q_parameter_search <- function(qs, q = seq(5, 95), v = seq(100, 1000, by = 25)) {
  grid <- expand.grid(
    q = paste0(q, "%"),
    v = v
  )
  q_parameter <- q
  grid <- grid %>% split(
    cut(seq(nrow(grid)), c(seq(0, nrow(grid), by=10), Inf))
  )
  grid <- grid %>% sapply(
    \(grid) grid %>%
      rowwise %>%
      mutate(
        germline = list(qs$germline[, q %>% as.character] %>% replace(is.na(.), 0) %>% `>=`(v)),
        somatic = list(qs$somatic[, q %>% as.character] %>% replace(is.na(.), 0) %>% `>=`(v))
      ) %>%
      mutate(
        smaller = min(sum(germline), sum(somatic)),
        larger = max(sum(germline), sum(somatic)),
        overlap = sum(germline & somatic) / smaller,
        smaller_unique = min(sum(germline & !somatic), sum(somatic & !germline))
      ) %>%
      ungroup %>%
      subset(select = -c(germline, somatic)) %>%
      mutate(q = q_parameter[as.numeric(q)]),
    simplify = FALSE
  ) %>%
    do.call(rbind, .)
  grid %>% ggplot(aes(q, v, fill=overlap)) + geom_tile() %>% rasterise(dpi = 300) + scale_fill_viridis_c() + coord_cartesian(expand=FALSE)
  # germline_genes = names(which(sctransform_quantile$germline[, "90%"] >= 1))
  grid
}