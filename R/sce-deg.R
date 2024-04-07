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

# Differentially expressed genes.
analyze_sce_to_csv <- function(seurats, output_path) {
  stopifnot(
    length(
      purrr::reduce(sapply(seurats, Cells, simplify=FALSE), intersect)
    ) == 0
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
    list(counts = do.call(cbind, seurats %>% sapply(\(seurat) GetAssayData(seurat, assay='RNA', layer='counts'), simplify=F))),
    colData = do.call(rbind, seurats %>% setNames(NULL) %>% sapply(\(seurat) seurat@meta.data, simplify=F))
  )
  colData(sce) = cbind(
    ident = do.call(
      c,
      sapply(seurats, Idents, simplify=F)
    ) %>% setNames(colnames(sce)),
    colData(sce)
  )
  # Bulk analysis: We will output cell barcodes, and these are the cells that we
  # will analyze 10X transcript-level quantification. We are also going to
  # compute % cells with each transcript present, so we are chiefly interested
  # in applying the same total UMI filter that we are going to apply later.
  # % Cells would be a more reproducible statistic if cells are more similar in
  # sequencing depth.
  colData(sce)$nCount_RNA_filter <- (
    between(sce$nCount_RNA, 2200, 7500)
  ) %>%
    factor %>%
    recode(`FALSE`="nCount_RNA_fail", `TRUE`="nCount_RNA_pass")

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

build_model_matrix <- function(metadata) {
  metadata <- metadata %>%
    read.csv(row.names = 1)
  # Factor levels should come from the filter_integrate_data function.
  metadata$ident <- metadata$ident %>%
    fct_relevel(c("germline", "somatic", "spermatocyte", "muscle"))
  metadata$batch <- metadata$batch %>% factor
  mm <- model.matrix(~ ident + batch, metadata)

  # Let Beta be a row vector:
  # [germline germlinecyst LFCtid LFCmuscle], where germline is the intercept
  # coefficient, and germlinecyst, LFCtid, and LFCmuscle are all deltas relative
  # to germline.
  # We want to add germlinecyst/2 to the first entry (intercept: germ/cyst avg)
  # and subtract it from the third and fourth entries (again, using
  # germ/cyst avg as the baseline for both outlier clusters).
  # X (in model matrix) is also a row vector, and the model is fitted such that:
  #   E[count] = Offset * beta * t(X).
  # If we want a different beta vector, then solve for the column operations
  # that we are going to apply: X -> X * solve(t(M)).
  # In columns: Update to Intercept, germlinecyst unchanged, subtract the new
  # Intercept from the values for the smaller clusters.
  solve(t(matrix(c(1,0.5,0,0, 0,1,0,0, 0,-0.5,1,0, 0,-0.5,0,1), nrow=4)))
  colnames(mm)[1:4] <- c("Mean", "CySCoverGSC", "TidOverMean", "MuscleOverMean")
  mm[, 1:4] <- mm[, 1:4] %*% (
    matrix(c(1,0.5,0,0, 0,1,0,0, 0,-0.5,1,0, 0,-0.5,0,1), nrow=4) %>%
      t %>%
      solve
  )
  mm[
    ,
    # Filtering of doublet clusters comes from the filter_integrate_data func.
    !grepl("doublet", colnames(mm))
  ]
}

fit_glm <- function(Upd_sc, model_matrix, metadata) {
  metadata <- metadata %>% read.csv(row.names = 1)
  Upd_glm <- glm_gp(
    GetAssayData(Upd_sc, assay="RNA", layer="counts"),
    model_matrix[colnames(Upd_sc), ],
    size_factors = metadata[colnames(Upd_sc), "size_factor"],
    on_disk = F,
    verbose = T
  )
  # Quite useful to keep the metadata in the GLM object.
  colData(Upd_glm$data)[, colnames(metadata)] <- metadata[
    colnames(Upd_glm$data),
  ]
  # We did not specify an offset per gene. Grab the offset per cell (column).
  # We can reconstruct a full-size Offset matrix later.
  Upd_glm$Offset = as.matrix(Upd_glm$Offset[1,, drop=F])
  # Mu is quite large and uncompressible, and can be recalculated using
  # glmGamPoi:::calculate_mu.
  Upd_glm$Mu = NULL
  Upd_glm
}

apeglm_coef_table <- function(g, coef=2) {
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

pseudobulk_cpm <- function(tx_files, size_factors, gtf_file) {
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
    length = abs(gtf$end - gtf$start) + 1
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

  # We have per-transcript values which are either "cpm" (bulk) or each
  # pseudobulk cluster name. TODO: Get something more dynamic here using
  # dplyr functions to aggregate transcript_id with each column.
  if ("cpm" %in% colnames(cpm_table))
    perform_summarize_transcripts <- \(df) df %>%
      summarize(
        cpm = transcript_id[which.max(cpm)]
      )
  else
    perform_summarize_transcripts <- \(df) df %>%
      summarize(
        germline = transcript_id[which.max(germline)],
        somatic = transcript_id[which.max(somatic)],
        spermatocyte = transcript_id[which.max(spermatocyte)],
        muscle = transcript_id[which.max(muscle)]
      )

  cluster_transcript <- apply(
    cpm_table %>% subset(select=-c(gene_id, length)),
    2,
    \(v) v / cpm_table$length
  ) %>%
    as.data.frame %>%
    rownames_to_column("transcript_id") %>%
    cbind(gene_id = cpm_table$gene_id) %>%
    group_by(gene_id) %>%
    perform_summarize_transcripts %>%
    right_join(
      data.frame(
        gene_id = metafeatures$flybase, rowname = rownames(metafeatures)
      ), by = "gene_id"
    ) %>%
    column_to_rownames

  colnames(fpkm_table)[1] <- "flybase"
  fpkm_table <- metafeatures %>%
    rownames_to_column %>%
    subset(select=c(rowname,flybase)) %>%
    left_join(fpkm_table, "flybase") %>%
    column_to_rownames %>%
    subset(select = -flybase) %>%
    as.matrix
  attr(fpkm_table, "transcript_id") <- cluster_transcript
  fpkm_table
}

glm_make_cpm_table <- function(Upd_glm) {
  cpm_table <- matrix(
    nrow = nrow(Upd_glm$Beta),
    ncol = nrow(sce.clusters),
    dimnames = list(rownames(Upd_glm$Beta), sce.clusters$cluster)
  )
  for (i in seq(nrow(sce.clusters))) {
    contrast <- sce.clusters[i, "contrast"] %>% unlist
    contrast_cols <- grep(
      "batch",
      colnames(Upd_glm$model_matrix),
      invert = TRUE
    )
    col_data = colData(Upd_glm$data)[
      colAlls(t(Upd_glm$model_matrix[contrast_cols, ]))
      == contrast[contrast_cols],
    ]
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

join_cpm_data <- function(
    Upd_cpm_transcripts, Upd_fpkm, Upd_cpm_regression) {
  Upd_fpkm_transcript_id <- (
    attr(Upd_fpkm, "transcript_id")[rownames(Upd_fpkm), ]
  )
  Upd_cpm_by_fpkm <- mapply(
    \(n, transcripts) Upd_cpm_transcripts[transcripts, ] %>%
      cbind(symbol = rownames(Upd_fpkm)) %>%
      pull(n, symbol),
    sce.clusters$cluster,
    Upd_fpkm_transcript_id %>%
      subset(select = c(germline, somatic, spermatocyte, muscle))
  )

  # Find genes that didn't have a transcript_id in the 10X BAM file, but which
  # were present in enough cells so that we included the gene when computing
  # Upd_cpm_regression.
  replace_genes <- is.na(Upd_cpm_by_fpkm) %>%
    rowAnys(useNames = TRUE) %>%
    which %>%
    names
  replace_genes <- intersect(replace_genes, rownames(Upd_cpm_regression))
  Upd_cpm_by_fpkm[replace_genes,] <- Upd_cpm_regression[replace_genes,]

  Upd_cpm_by_fpkm
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