# Contrasts matrix for the cell type idents.
Upd_celltype_contrasts <- matrix(
  c(0.5, 0.5, 0, 0, 0,
    -1, 1, 0, 0, 0,
    -1, 0, 1, 0, 0,
    0, -1, 0, 1, 0,
    0, -1, 0, 0, 1),
  nrow = 5,
  byrow = TRUE,
  dimnames = list(
    c("germline", "somatic", "spermatocyte", "somaticprecursor", "muscle"),
    c("Mean", "CySCoverGSC", "TidOverGSC", "SPOverCySC", "MuscleOverCySC")
  )
)
Upd_celltype_model_matrix <- matrix(
  Upd_celltype_contrasts %>% solve,
  nrow = 5,
  dimnames = dimnames(Upd_celltype_contrasts)
)

# Find an nCount_RNA threshold so that % expressed is more comparable between
# the Nos and tj-driven samples. This is saved in the csv metadata and is
# applied later for a supplemental dot plot.
plot_nCount_RNA_threshold <- function(Upd_sc, ident, gene) {
  Upd_sc$batch <- Upd_sc$batch %>% factor
  data <- levels(Upd_sc$batch) %>%
    sapply(
      \(n) Upd_sc@meta.data %>%
        subset(batch == n & Idents(Upd_sc) == ident) %>%
        cbind(is.expressed = Upd_sc[['RNA']]@layers$counts[
          gene, Upd_sc$batch == n & Idents(Upd_sc) == ident
        ] != 0) %>%
        arrange(desc(nCount_RNA)) %>%
        cbind(num.cells = seq_along(.$is.expressed), num.expressed = cumsum(.$is.expressed)) %>%
        cbind(expression.rate = .$num.expressed / .$num.cells),
      simplify = FALSE)
  thresholds = seq(1000, 5000, length.out = 101) %>% round
  threshold_stat <- thresholds %>%
    sapply(
      \(t) data %>% sapply(
        \(df) df[
          nrow(df) + 1 - findInterval(t, c(0, rev(df$nCount_RNA))),
          "expression.rate"
        ]
      )
    )
  colnames(threshold_stat) <- thresholds
  names(dimnames(threshold_stat)) <- c("batch", "nCount_RNA")
  melt(threshold_stat) %>% ggplot(aes(nCount_RNA, value, group=batch, color=batch)) + geom_line()
  data.frame(thresholds, sd = colSds(threshold_stat))
}
plot_nCount_RNA_threshold_summary <- function() {
  mapply(plot_nCount_RNA_threshold, rep(list(Upd_sc), 2), c("germline","germline","germline","somatic","somatic","somatic"), c("vas","RpL22-like","blanks","tj","lncRNA:roX1","zfh1"), SIMPLIFY=F) %>% bind_rows(.id = "gene") %>% ggplot(aes(thresholds,sd,group=gene,color=gene)) + geom_line()
}

pooled_size_factors_seurat_ident_autosome <- function(seurat) {
  sce <- GetAssayData(seurat, assay = "RNA", layer = "counts") %>%
    list(counts = .) %>%
    SingleCellExperiment
  sce <- sce[seurat[["RNA"]]@meta.data$chr %in% c("2L", "2R", "3L", "3R"), ]
  pooledSizeFactors(
    sce,
    clusters = Idents(seurat)
  )
}

extract_upd_metadata_to_csv <- function(Upd_sc, size_factors, output_path) {
  Upd_sc$size_factors <- size_factors
  meta.data <- cbind(
    FetchData(Upd_sc, "ident"),
    Upd_sc@meta.data,
    nCount_RNA_filter = "nCount_RNA_pass"
  )
  write.csv(meta.data, output_path)
  output_path
}

build_model_matrix <- function(data_frame) {
  contrasts(data_frame$ident) <- Upd_celltype_model_matrix[, -1]
  model.matrix(~ ident + batch, data_frame)
}

fit_glm <- function(Upd_sc, Upd_model_matrix, Upd_metadata) {
  Upd_metadata <- Upd_metadata %>% read.csv(row.names = 1)
  Upd_glm <- glm_gp(
    GetAssayData(Upd_sc, assay="RNA", layer="counts"),
    Upd_model_matrix,
    size_factors = Upd_metadata[Cells(Upd_sc), "size_factors"],
    on_disk = F,
    verbose = T
  )
  # Quite useful to keep the metadata in the GLM object.
  colData(Upd_glm$data)[, colnames(Upd_metadata)] <- (
    Upd_metadata[colnames(Upd_glm$data), ]
  )
  # We did not specify an offset per gene. Grab the offset per cell (column).
  # We can reconstruct a full-size Offset matrix later.
  Upd_glm$Offset = as.matrix(Upd_glm$Offset[1,, drop=F])
  # Mu is quite large and uncompressible, and can be recalculated using
  # glmGamPoi:::calculate_mu.
  Upd_glm$Mu = NULL
  Upd_glm
}

apeglm_coef_table <- function(g, coef=2, test_mle=TRUE) {
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
  if (test_mle) {
    message('F-test of regression without shrinkage')
    mle <- test_de(g, colnames(g$Beta)[coef])
  } else {
    mle <- NULL
  }
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
  # Our Upd_glm generally has enough cells in each cluster that the se.fit is
  # never NA. If we are testing an individual cluster other than the 4 factor
  # levels used so far for fitting Upd_glm, then we might have NA or Inf and
  # this could cause an error but can easily be excluded.
  lfc[rowAnys(!is.finite(lfc)), ] <- NA
  # Garbage-collect Mu.
  g$Mu = NULL
  message('Optimize Approximate Posterior')
  lfc = apeglm(
    assay(g$data),
    g$model_matrix,
    log.lik=NULL,
    param=g$overdispersion_shrinkage_list$dispersion_trend,
    coef=coef,
    mle=lfc,
    method='nbinomCR',
    offset=g$Offset
  )
  if (test_mle) lfc$mle.test <- mle
  lfc
}

small_cluster_test_glm <- function(Upd_glm, metadata, seurat_qc_obj, batch, level) {
  metadata <- metadata %>% read.csv(row.names = 1)
  known_cells_take <- metadata %>%
    rownames_to_column %>%
    filter(
      ident %in% c("germline", "somatic"), batch == {{batch}}
    ) %>%
    pull(rowname) %>%
    # The cells to include are determined by the filtering done in the
    # filter_integrate_data function.
    intersect(colnames(Upd_glm$data))
  if (length(intersect(known_cells_take, Cells(seurat_qc_obj)[Idents(seurat_qc_obj) == level])) > 0)
    return(NA)

  sce <- SingleCellExperiment(
    list(counts=assay(Upd_glm$data, "counts")[, known_cells_take]),
    colData = metadata[known_cells_take, ]
  )
  rowData(sce)$overdispersions = Upd_glm$overdispersions
  # Counts matrix for the level of interest.
  counts_level <- GetAssayData(seurat_qc_obj, assay="RNA", layer="counts")[
    rownames(sce),
    Idents(seurat_qc_obj) == level
  ]
  # Only interested in genes where we can show positively as a marker for the
  # level of interest.
  genes_take <- names(which(rowSums(counts_level != 0) >= 20))
  sce <- sce[genes_take, ]
  counts_level <- counts_level[genes_take, ]

  sce <- cbind(
    sce,
    SingleCellExperiment(
      list(counts = counts_level),
      colData = metadata[Cells(seurat_qc_obj)[Idents(seurat_qc_obj) == level], ]
    )
  )
  sce$ident <- sce$ident %>% fct_relevel("germline", "somatic")
  mm <- model.matrix(~ 0 + ident, colData(sce))
  # Matrix of the contrasts
  xform <- matrix(
    c(1, 1, 1,
      -0.5, 0.5, 0,
      0, 0, 1),
    ncol = 3,
    dimnames = list(
      # Rownames - disappear after multiplying by xform on the right
      c("germline", "somatic", "other"),
      # Names of the contrasts
      c("Mean", "CySCoverGSC", "OtherOverMean")
    )
  )
  mm <- mm %*% xform
  glm_gp(
    sce,
    ~ 0 + mm,
    size_factors = sce$size_factor,
    overdispersion = rowData(sce)$overdispersions,
    on_disk = F,
    verbose = T
  )
}

small_cluster_test_apeglm <- function(Upd_glm, metadata, seurat_qc_obj, batch, level, ...) {
  result <- small_cluster_test_glm(Upd_glm, metadata, seurat_qc_obj, batch, level)
  if (!is.list(result))
    return(NA)
  # Contrasts: GSC+CySC, CySCoverGSC, OtherOverMean
  apeglm_coef_table(result, coef=3, ...)
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

# Pseudobulk the filtered exonic UMIs from the Cell Ranger BAM file.
pseudobulk_cpm <- function(tx_files, size_factors, gtf_file) {
  cluster_name = factor(
    tx_files$tx_file %>% str_extract('([a-z]+)\\.txt', group=1),
    sce.clusters$cluster
  )
  if (all(!is.na(cluster_name)))
    pseudobulk_colData = cbind(size_factors, ident = cluster_name)
  else
    pseudobulk_colData = cbind(size_factors, ident = 'intercept')
  # Use nos.2 as reference level for CPM calculation. All 4 batches have plenty
  # of GSC and CySC cells, but nos.2 has the most cyte/tid and eya+ somatic
  # cells so we want to quantify those cells clearly.
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

  gtf.all.types <- read.table(
    gtf_file,
    sep = '\t',
    col.names = c('chr', 'source', 'type', 'start', 'end', 'sc', 'strand', 'fr', 'annotation'),
    header = F,
    quote = ''
  )
  gtf.all.types$tx <- gtf.all.types$annotation %>% str_extract(
    'transcript_id "([^"]+)"',
    group = 1
  )
  gtf.all.types$gene_id <- gtf.all.types$annotation %>% str_extract(
    'gene_id "([^"]+)"',
    group = 1
  )
  # include "transcript" types in the reference
  gtf <- gtf.all.types %>% subset(grepl('RNA', type))
  cpm_table <- data.frame(
    tx = gtf$tx,
    gene_id = gtf$gene_id,
    cds_length = abs(gtf$end - gtf$start) + 1
  )

  # Now calculate the exon "hitbox" for Cell Ranger UMIs of type "E" (exonic).
  # First, try summing the exon lengths. Second, using Cell Ranger documentation
  # that the read (known length 50 bp) must have at least 50% overlap with the
  # exon, we create a "hitbox" by removing half the read length from each end.
  exons <- gtf.all.types %>% subset(type == "exon") %>% group_by(tx)
  cpm_table <- cpm_table %>%
    left_join(
      exons %>% summarise(exon_length = sum(abs(end - start) + 1)),
      "tx"
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

glm_make_cpm_table <- function(Upd_glm) {
  cpm_table <- matrix(
    nrow = nrow(Upd_glm$Beta),
    ncol = nrow(sce.clusters),
    dimnames = list(rownames(Upd_glm$Beta), sce.clusters$cluster)
  )
  for (i in seq(nrow(sce.clusters))) {
    # Contrast is a "short" vector and is missing the later model matrix columns
    # such as batch effect.
    contrast <- sce.clusters[i, "contrast"] %>% unlist
    contrast_cols <- seq(length(contrast))
    col_data <- cbind(colData(Upd_glm$data), size_factor = Upd_glm$size_factors)
    col_data = col_data[
      colAlls(
        t(Upd_glm$model_matrix)[contrast_cols, ]
        == contrast[contrast_cols]
      ),
    ]
    # Now we do not need to model batch effect when predicting cpm.
    contrast <- c(
      contrast, rep(0, ncol(Upd_glm$model_matrix) - length(contrast))
    )
    # Predict the sum of counts where the offset is the sum of size factors for
    # the cells. Normalize by sum of UMI count. Multiply by 1MM yielding counts
    # per million.
    offset_predictions = (
      log(sum(col_data$size_factor))
      - log(sum(col_data$nCount_RNA))
      + log(1000) * 2
    )
    cpm_table[, i] <- predict(
      Upd_glm,
      matrix(contrast, nrow=1),
      type='response',
      offset=offset_predictions
    )
  }
  cpm_table
}

cpm_to_fpkm_using_cds <- function(cpm_table, assay_data_sc) {
  assay_data_sc <- assay_data_sc %>% read.csv(row.names = 1)
  cpm_table %>%
    `/`(assay_data_sc[rownames(cpm_table), "transcript.length"] / 1000) %>%
    replace(!is.finite(.), 0)
}

cpm_select_transcript_to_quantify <- function(
  cpms, assay_data, correct_for_exon_length = T
) {
  assay_data <- assay_data %>% read.csv(row.names = 1)
  cpms$gene_id <- cpms$gene_id %>% factor(assay_data$flybase)
  cpms %>%
    rownames_to_column %>%
    group_by(gene_id) %>%
    summarise(
      germline = rowname[which.max(germline / if (correct_for_exon_length) exon_length else 1)],
      somatic = rowname[which.max(somatic / if (correct_for_exon_length) exon_length else 1)],
      spermatocyte = rowname[which.max(spermatocyte / if (correct_for_exon_length) exon_length else 1)],
      somaticprecursor = rowname[which.max(somaticprecursor / if (correct_for_exon_length) exon_length else 1)],
      muscle = rowname[which.max(muscle / if (correct_for_exon_length) exon_length else 1)]
    ) %>%
    mutate(gene_id = gene_id %>% replace_na("")) %>%
    column_to_rownames(var = "gene_id")
}

lookup_mat_transcripts_exon_length <- function(Upd_cpm_transcripts, Upd_cpm_transcript_to_use, assay.data.sc) {
  assay.data.sc <- assay.data.sc %>% read.csv(row.names = 1)
  exon_length <- matrix(
    Upd_cpm_transcripts[Upd_cpm_transcript_to_use[assay.data.sc$flybase, ] %>% as.matrix %>% as.character %>% cbind("exon_length")] %>%
      # character to numeric
      as.numeric,
    ncol=5,
    dimnames=list(rownames(assay.data.sc), colnames(Upd_cpm_transcript_to_use))
  )
}

join_cpm_data <- function(
  assay_data,
  Upd_cpm_transcripts,
  Upd_cpm_transcript_to_use,
  Upd_cpm_regression,
  correct_for_exon_length = T
) {
  assay_data <- read.csv(assay_data, row.names = 1)
  Upd_cpm_transcript_to_use <- left_join(
    assay_data %>% rownames_to_column %>% reframe(rowname, flybase),
    Upd_cpm_transcript_to_use %>%
      rownames_to_column(var = "flybase") %>%
      subset(flybase != ""),
    "flybase"
  )
  Upd_cpm <- mapply(
    \(tnames, transcripts, transcript_lookup) transcripts[
      match(transcript_lookup, tnames)
    ],
    list(rownames(Upd_cpm_transcripts)),
    Upd_cpm_transcripts %>%
      subset(select = c(germline, somatic, spermatocyte, somaticprecursor, muscle)) %>%
      `/`(if (correct_for_exon_length) Upd_cpm_transcripts$exon_length else 1),
    Upd_cpm_transcript_to_use %>%
      subset(select = c(germline, somatic, spermatocyte, somaticprecursor, muscle))
  )
  dimnames(Upd_cpm) = list(
    Upd_cpm_transcript_to_use$rowname,
    colnames(Upd_cpm_transcript_to_use)[3:7]
  )

  # Find genes that didn't have a transcript_id in the 10X BAM file, but which
  # were present in enough cells so that we included the gene when computing
  # Upd_cpm_regression.
  replace_genes <- is.na(Upd_cpm) %>%
    rowAnys(useNames = TRUE) %>%
    which %>%
    names
  replace_genes <- intersect(replace_genes, rownames(Upd_cpm_regression))
  Upd_cpm[replace_genes,] <- Upd_cpm_regression[
    replace_genes,
    c("germline", "somatic", "spermatocyte", "somaticprecursor", "muscle")
  ]

  Upd_cpm
}

table_to_tpm <- function(quant_table) {
  quant_table %>%
    t %>%
    `*`(1000 * 1000 / rowSums(., na.rm=T)) %>%
    t
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

dot_plot_cpm <- function(genotypes, gene_list, logcpm_min=0, logcpm_max=3, oob_squish=FALSE) {
  data <- genotypes %>%
    sapply(\(df) df[gene_list, ] %>% rownames_to_column("gene"), simplify = FALSE) %>%
    bind_rows(.id = "genotype")
  # Show gene in order of appearance in gene_list
  data$gene <- data$gene %>% factor(levels = unique(.)) %>%
    fct_relabel(\(v) v %>% str_replace('lncRNA:', '') %>% str_replace('Hsromega', 'Hsr\u03C9'))
  data$genotype <- data$genotype %>% factor(levels = names(genotypes))
  g <- data %>% ggplot(
    aes(x = genotype, y = gene, size = percent_expressed, color = log(cpm)/log(10))
  ) + geom_point()
  if (oob_squish)
    g <- g + scale_color_viridis_c(
      begin = 0.15,
      limits = c(logcpm_min, logcpm_max),
      oob = squish
    )
  else
    g <- g + scale_color_viridis_c(
      begin = 0.15,
      limits = c(logcpm_min, logcpm_max)
    )
  g + scale_y_discrete(limits=rev) + scale_size_continuous(
    labels=percent, limits=c(0.01,0.99), range=c(1, 10)
  ) + labs(
    x = NULL, y = NULL, color = bquote(log[10]*"(CPM)"), size = "expression"
  ) + theme_cowplot()
}

# Another cluster-gene quantification method: Take the quantile of interest for
# the scaled data.
quantify_quantiles <- function(
  seurat, metadata, assay = "SCT", slot = "scale.data",
  per_million = FALSE
) {
  if (assay == "RNA" && slot == "data")
    seurat <- seurat %>% NormalizeData(assay = assay)
  assay_data <- as.matrix(GetAssayData(seurat, assay = assay, layer = slot))
  assay_cells <- data.frame(rowname = colnames(assay_data))
  metadata <- metadata %>% read.csv(row.names = 1) %>%
    rownames_to_column %>%
    inner_join(assay_cells, by="rowname")
  cells <- sce.clusters$cluster %>%
    sapply(
      \(n) assay_data[
        ,
        metadata %>%
          subset(
            as.character(ident) == n & nCount_RNA_filter == "nCount_RNA_pass"
          ) %>%
          pull(rowname)
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
  sapply(
    sce.clusters$cluster,
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