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
# the nos and tj-driven samples. This is saved in the csv metadata and is
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

pooled_size_factors_seurat_ident_autosome <- function(seurat, assay = "RNA") {
  sce <- GetAssayData(seurat, assay = assay, layer = "counts") %>%
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

build_model_matrix <- function(data_frame, decontXcontam) {
  contrasts(data_frame$ident) <- Upd_celltype_model_matrix[, -1]
  idents <- model.matrix(~ ident + batch, data_frame)[
    ,
    c("(Intercept)", "identCySCoverGSC", "batchnos.2", "batchtj.1", "batchtj.2")
  ]
  colnames(idents)[1] <- ""
  decontXcontam <- rowSums2(decontXcontam)
  model.matrix(~ ident + batch + idents : decontXcontam, data_frame)
  # model.matrix(~ ident + batch, data_frame)
}

fit_glm <- function(Upd_assay, Upd_model_matrix, Upd_metadata) {
  Upd_metadata <- Upd_metadata %>% read.csv(row.names = 1)
  # Upd_assay <- Upd_assay %>% subset(cells = Cells(.)[Upd_subset])
  # Upd_model_matrix <- Upd_model_matrix[Upd_subset, ]
  # Upd_metadata <- Upd_metadata[Upd_subset, ]
  Upd_glm <- glm_gp(
    GetAssayData(Upd_assay, layer="counts"),
    Upd_model_matrix,
    size_factors = Upd_metadata[Cells(Upd_assay), "size_factors"],
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

apeglm_coef_timebound <- function(g, coef=2, prior_var) {
  message('Optimize Approximate Posterior')
  assay <- as.matrix(assay(g$data))
  Offset <- g$Offset[1,, drop=T]
  basemean <- rowMeans(assay)
  # From apeglm source code.
  prior.control <- list(
    no.shrink = setdiff(seq(ncol(g$Beta)), coef),
    prior.mean = 0,
    prior.scale = if (is.null(prior_var)) apeglm:::priorVar(lfc) else prior_var,
    prior.df = 1,
    prior.no.shrink.mean = 0,
    prior.no.shrink.scale = 15
  )
  # apeglm sets the prior scale parameter to a min of 1 (coefficient represents 
  # a fold-change of e) or empirical value. Wow such interpretable, very cool.
  prior.control$prior.scale <- prior.control$prior.scale %>% min(1)
  apeglm_list <- future_lapply(
    seq(nrow(g$data)),
    \(i) tryCatch(
      withTimeout(
        apeglm:::apeglm.single(
          y = assay[i, ],
          x = g$model_matrix,
          log.lik = logLikNB,
          param = g$overdispersion_shrinkage_list$dispersion_trend[i],
          coef = coef,
          method = "general",
          offset = Offset,
          prefit.beta = NULL,
          basemean = basemean[i],
          optim.method = "BFGS",
          weights = NULL,
          prior.control = prior.control,
          prefit.conv = NULL,
          interval.type = "laplace",
          param.sd = NULL,
          interval.level = 0.95,
          threshold = NULL
        ),
        timeout = 30
      ),
      TimeoutException = \(ex) NULL
    )
  )
  names(apeglm_list) <- rownames(g$data)
  apeglm_list <- apeglm_list
}

apeglm_coef_table_sample <- function(g, coef=2, test_mle=TRUE, prior_var = NULL, shrinkage_cutoff = 1) {
  message('Construct dense matrices in memory for regression')
  g$Offset = matrix(
    g$Offset[1,, drop=T],
    byrow=T,
    nrow=nrow(g$data),
    ncol=ncol(g$data),
    dimnames=dimnames(g$data)
  )
  if (is.null(prior_var) || test_mle) {
    g$Mu = glmGamPoi:::calculate_mu(
      Beta=g$Beta,
      model_matrix=g$model_matrix,
      offset_matrix=g$Offset
    )
  }
  if (test_mle) {
    message('F-test of regression without shrinkage')
    mle <- test_de(g, colnames(g$Beta)[coef])
  } else {
    mle <- NULL
  }
  message('Calculate regression s.e. using QR decomposition')
  if (is.null(prior_var)) {
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
  }
  # From apeglm source code.
  prior.control <- list(
    no.shrink = setdiff(seq(ncol(g$Beta)), coef),
    prior.mean = 0,
    prior.scale = if (is.null(prior_var)) apeglm:::priorVar(lfc) else prior_var,
    prior.df = 1,
    prior.no.shrink.mean = 0,
    prior.no.shrink.scale = 15
  )
  # apeglm sets the prior scale parameter to a min of 1 (coefficient represents 
  # a fold-change of e) or empirical value. Wow such interpretable, very cool.
  prior.control$prior.scale <- prior.control$prior.scale %>% min(1)
  if (shrinkage_cutoff != 0)
    genes.test <- (
      rowMins(abs(lfc[,1] + qnorm(0.975) * lfc[,2] %*% t(c(-1, 1)))) / log(2) > shrinkage_cutoff
    ) %>%
      which
  else
    genes.test <- rownames(g$data)
  # Garbage-collect Mu.
  g$Mu = NULL
  message('Optimize Approximate Posterior')
  lfc = apeglm(
    assay(g$data)[genes.test, ],
    g$model_matrix,
    log.lik=NULL,
    param=g$overdispersion_shrinkage_list$dispersion_trend[genes.test],
    coef=coef,
    method='nbinomCR',
    offset=g$Offset[genes.test, ],
    prior.control = prior.control
  )
  lfc$map <- lfc$map[match(seq(nrow(g$data)), genes.test), ]
  lfc$sd <- lfc$sd[match(seq(nrow(g$data)), genes.test), ]
  lfc$fsr <- lfc$fsr[match(seq(nrow(g$data)), genes.test),, drop=F]
  if (test_mle) lfc$mle.test <- mle
  lfc
}

apeglm_coef_table_prior_control <- function(sce, g, coef=2, test_mle=TRUE) {
  offset_matrix <- (
    log(assay(sce, "counts"))
    - log(assay(sce, "decontXcounts"))
    # + rep(log(sce$size_factors), each = nrow(assay(sce)))
  ) %>%
    replace(!is.finite(.), 0)
  offset_prior_control <- -solve(
    crossprod(g$model_matrix, g$model_matrix),
    t(offset_matrix %*% g$model_matrix)
  ) %>% t

  offset_size_factors_matrix <- matrix(
    log(sce$size_factors),
    nrow = nrow(sce),
    ncol = ncol(sce),
    byrow = TRUE,
    dimnames = dimnames(sce)
  )

  if (test_mle) {
    message('F-test of regression without shrinkage')
    mle <- test_de(g, colnames(g$Beta)[coef], pval_adj="holm")
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

  message('Calculate Prior Distribution')
  prior.control <- tibble(
    no.shrink = list(setdiff(seq(ncol(g$Beta)), coef)),
    prior.mean = sapply(
      offset_prior_control[,coef],
      \(n) rep(0, ncol(g$Beta)) %>% replace(coef, n),
      simplify=FALSE
    ),
    prior.scale = apeglm:::priorVar(lfc),
    prior.df = 1,
    prior.no.shrink.mean = 0,
    prior.no.shrink.scale = 15,
    my_gene_name = rownames(offset_prior_control)
  )

  message('Optimize Approximate Posterior')
  apeglm_results = apply(
    prior.control[1:18,, drop=F],
    1,
    \(prior.control.lst) apeglm(
      assay(g$data)[prior.control.lst$my_gene_name,, drop=F],
      g$model_matrix,
      log.lik=logLikNB,
      param=g$overdispersion_shrinkage_list$dispersion_trend,
      coef=coef,
      offset=offset_size_factors_matrix[prior.control.lst$my_gene_name,, drop=F],
      prior.control = prior.control.lst
    ),
    simplify=F
  )
  lfc = apeglm(
    assay(g$data),
    g$model_matrix,
    log.lik=logLikNB,
    param=g$overdispersion_shrinkage_list$dispersion_trend,
    coef=coef,
    mle=lfc,
    offset=g$Offset
  )
  if (test_mle) lfc$mle.test <- mle
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
      `/`(if (correct_for_exon_length) Upd_cpm_transcripts$exon_length else 1) %>%
      as.data.frame,
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
  replace_genes <- replace_na(Upd_cpm == 0, replace=T) %>%
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

sce_select_transcript_to_use <- function(
  Upd_sc, assay.data.sc,
  Upd_sc_pseudobulk_transcripts, Upd_cpm_transcript_to_use
) {
  assay.data.sc <- assay.data.sc %>% read.csv(row.names = 1)

  # For genes such as H3-GFP,
  # use the Cell Ranger include-introns counts without
  # decontamination.
  counts <- GetAssayData(Upd_sc, "RNA", "counts")
  counts <- pseudobulk(
    SingleCellExperiment(list(counts=counts)), vars(ident, batch),
    col_data = FetchData(Upd_sc, c("ident", "batch"))
  )
  stopifnot(all.equal(counts$ident, Upd_sc_pseudobulk_transcripts$ident))
  stopifnot(all.equal(as.character(counts$batch), as.character(Upd_sc_pseudobulk_transcripts$batch)))
  counts <- counts %>% assay
  decontXcounts <- counts

  dest_row <- rep(match(rownames(Upd_cpm_transcript_to_use), assay.data.sc$flybase), 20)
  dest_col <- rep(1:20, each = length(dest_row) / 20)
  source_row <- rep(match(unlist(Upd_cpm_transcript_to_use), rownames(Upd_sc_pseudobulk_transcripts)), 4)
  source_col <- dest_col

  dest_row <- match(rownames(assay.data.sc)[dest_row], rownames(Upd_sc))
  mask_lookups <- !is.na(dest_row)
  dest_row <- dest_row[mask_lookups]
  dest_col <- dest_col[mask_lookups]
  source_row <- source_row[mask_lookups]
  source_col <- source_col[mask_lookups]

  counts[cbind(dest_row, dest_col)] <- assay(Upd_sc_pseudobulk_transcripts, "counts")[cbind(source_row, source_col)]
  decontXcounts[cbind(dest_row, dest_col)] <- assay(Upd_sc_pseudobulk_transcripts, "decontXcounts")[cbind(source_row, source_col)]
  SingleCellExperiment(
    list(counts=counts, decontXcounts=decontXcounts),
    colData = colData(Upd_sc_pseudobulk_transcripts)
  )
}

table_to_tpm <- function(quant_table) {
  quant_table %>%
    t %>%
    `*`(1000 * 1000 / rowSums(., na.rm=T)) %>%
    t
}
