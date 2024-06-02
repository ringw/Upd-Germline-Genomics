excel_tables = list(
  deepgreen='TableStyleMedium18',
  cyan='TableStyleMedium6',
  orange='TableStyleMedium14',
  red='TableStyleMedium3',
  mono='TableStyleMedium1'
)

publish_excel_results <- function(
  Upd_regression_somatic,
  Upd_cpm_transcripts, Upd_cpm, Upd_tpm_do_not_use_for_quantification,
  Upd_isoform_exonic_length,
  sctransform_quantile, supplemental_bulk_pct_expressed, metafeatures, gtf_path,
  batch_data,
  target_path
) {
  wb = createWorkbook()

  metafeatures = metafeatures %>% read.csv(row.names=1)
  write_abundance_table(
    wb,
    'Gene Quantification',
    metafeatures$flybase,
    Upd_cpm,
    sctransform_quantile,
    supplemental_bulk_pct_expressed,
    excel_tables$mono
  )

  write_quant_stats(
    wb, "Gene Quantification Summary", Upd_cpm_transcripts, Upd_cpm,
    Upd_tpm_do_not_use_for_quantification, Upd_isoform_exonic_length,
    Upd_regression_somatic, sctransform_quantile, gtf_path, metafeatures
  )

  baseMeanCPM <- exp(Upd_regression_somatic$map[, "(Intercept)"]) %>% `*`(
    1000 * 1000 / sum(., na.rm=T)
  ) %>% signif(digits=2)
  flybase <- metafeatures %>%
    rownames_to_column %>%
    subset(pass.min.cells) %>%
    pull(flybase, rowname)
  write_regression_table(wb, 'GSC&CySC Regression', 'log(CySC/GSC)', baseMeanCPM, Upd_regression_somatic, 2, flybase, excel_tables$cyan)

  write_batch_idents(wb, "Cluster&Batch Statistics", batch_data)

  dir.create('scRNA-seq-Regression', showW=F)
  saveWorkbook(wb, target_path, overwrite = T)
  target_path
}

write_abundance_table <- function(
    wb, title, flybase, Upd_cpm, sctransform_quantile, supplemental_bulk_pct_expressed, table_style) {
  rename = c(germline='GSC', somatic="CySC", spermatocyte="tid", somaticprecursor="SP", muscle="muscle")
  data = data.frame(
    symbol = rownames(Upd_cpm),
    flybase = flybase
  )

  for (n in names(rename)) {
    name = rename[n]
    data[str_glue("CPM_{name}")] <- Upd_cpm[,n]
    if (n %in% names(sctransform_quantile)) {
      sct_ranking <- sctransform_quantile[[n]][, '90%']
      sct_ranking <- sct_ranking %>% round(2)
      sct_ranking_data <- data.frame(
        symbol = names(sct_ranking), v = sct_ranking
      )
      colnames(sct_ranking_data)[2] <- paste('SCT90%', name, sep='_')
      data <- data %>%
        left_join(
          sct_ranking_data,
          by = "symbol"
        )
    }
  }
  data = data %>%
    left_join(
      data.frame(
        symbol = rownames(supplemental_bulk_pct_expressed),
        pctAllCells = supplemental_bulk_pct_expressed$percent_expressed
      ),
      by = "symbol"
    )
  data$pctAllCells <- data$pctAllCells %>% replace(is.na(.), 0)
  class(data$pctAllCells) <- "percentage"
  addWorksheet(wb, title)
  writeDataTable(wb, title,
                 data,
                 startCol = 1, startRow = 1,
                 withFilter = T,
                 tableStyle = table_style)
}

write_quant_stats <- function(
  wb, title, Upd_cpm_transcripts, Upd_cpm,
  Upd_tpm_do_not_use_for_quantification, Upd_isoform_exonic_length,
  Upd_regression_somatic, sctransform_quantile, gtf_path, metafeatures
) {
  addWorksheet(wb, title)

  n_genes <- tribble(
    ~ cluster, ~ CPM, ~ SCT, ~ L2FC,
    "GSC",
    sum((Upd_cpm[, "germline"] > 0) %>% replace(is.na(.), FALSE)),
    sum(is.finite(sctransform_quantile[["germline"]][, "90%"])),
    sum(is.finite(Upd_regression_somatic$map[, 2])),
    "CySC",
    sum((Upd_cpm[, "somatic"] > 0) %>% replace(is.na(.), FALSE)),
    sum(is.finite(sctransform_quantile[["somatic"]][, "90%"])),
    sum(is.finite(Upd_regression_somatic$map[, 2]))
  )

  writeData(wb, title, paste0("# Genes Quantified by Method (Genes in Reference: ", nrow(Upd_cpm), ")"), colNames = F, startRow = 1, startCol = 1)
  writeDataTable(wb, title,
                 n_genes, withFilter = FALSE,
                 startCol = 1, startRow = 2)

  start_row <- 6
  group_names <- c(germline="GSC", somatic="CySC")
  Upd_cpm_transcripts$gene_id <- Upd_cpm_transcripts$gene_id %>% factor

  analyze_transcripts <- (as.matrix(Upd_cpm_transcripts[,-1]) > 0) %>%
    rowAlls(useNames=T) %>%
    which %>%
    names
  gene_id_transcripts <- Upd_cpm_transcripts[analyze_transcripts, "gene_id"] %>%
    droplevels
  exons <- log(Upd_cpm_transcripts[analyze_transcripts, "exon_length"])
  for (n in names(group_names)) {
    writeData(
      wb,
      title,
      matrix(
        c(
          str_glue("{group_names[n]}: CPM or TPM for Isoform Calling (avoid isoform corr. with exon_length)"),
          "Beta (coef. for exon_length) similar to Pearson's R is shown."
        ),
        ncol = 1
      ),
      colNames = F,
      startRow = start_row, startCol = 1
    )
    start_row <- start_row + 2
    logCPM <- log(Upd_cpm_transcripts[analyze_transcripts, n])
    logTPM <- log(Upd_cpm_transcripts[analyze_transcripts, n]) - exons
    # Fit slope and intercept to isoform length and CPM:
    # germline ~ gene_id + exon_length
    lme.cpm <- lmer(logCPM ~ exons + (1 | gene_id_transcripts))
    lme.tpm <- lmer(logTPM ~ exons + (1 | gene_id_transcripts))
    # Calculate Beta (standardize the exons coefficient) which is similar to a
    # Pearson correlation estimate between log(exon_length) and log(CPM).
    beta <- data.frame(
      CPM=coef(lme.cpm)[[1]][1, "exons"] * sd(exons) / sd(logCPM),
      TPM=coef(lme.tpm)[[1]][1, "exons"] * sd(exons) / sd(logTPM)
    )
    colnames(beta) <- str_glue(
      "log({group_names[n]}_{colnames(beta)}) ~ gene + log(exon_length)"
    )
    writeDataTable(
      wb,
      title,
      beta,
      startRow = start_row, startCol = 1
    )
    start_row <- start_row + 3
  }

  group_names <- c(germline="GSC", somatic="CySC")
  analyze_genes <- rownames(Upd_tpm_do_not_use_for_quantification) %>%
    setdiff(rownames(Upd_tpm_do_not_use_for_quantification) %>% subset(is.na(Upd_tpm_do_not_use_for_quantification[, "germline"]))) %>%
    setdiff(rownames(Upd_tpm_do_not_use_for_quantification) %>% subset(is.na(Upd_tpm_do_not_use_for_quantification[, "somatic"]))) %>%
    setdiff(names(which(rowAnys(Upd_tpm_do_not_use_for_quantification == 0, na.rm=T, useNames=T)))) %>%
    setdiff(names(which(rowAnys(is.na(Upd_isoform_exonic_length), useNames=T)))) %>%
    setdiff(names(which(rowAnys((Upd_isoform_exonic_length == 0), na.rm=T, useNames=T)))) %>%
    intersect(rownames(Upd_regression_somatic$map) %>% subset(is.finite(Upd_regression_somatic$map[, 2])))
  gtf <- read.table(
    gtf_path,
    sep = '\t',
    col.names = c('chr', 'source', 'type', 'start', 'end', 'sc', 'strand', 'fr', 'annotation'),
    header = F,
    quote = ''
  ) %>% subset(grepl('RNA', type)) # include "transcript" types in the reference
  gtf$tx <- gtf$annotation %>% str_extract(
    'transcript_id "([^"]+)"',
    group = 1
  )
  gtf$length = abs(gtf$end - gtf$start) + 1

  for (n in c("germline", "somatic")) {
    writeData(
      wb,
      title,
      paste0(group_names[n], ": Gene CPM or TPM (exon-length-normalized) to avoid corr. with Exon Length"),
      colNames = F,
      startRow = start_row, startCol = 1
    )
    start_row <- start_row + 1
    CPM <- Upd_cpm[analyze_genes, n]
    TPM <- Upd_tpm_do_not_use_for_quantification[analyze_genes, n]
    SCT = sctransform_quantile[[n]][analyze_genes, "90%"]
    lengths <- Upd_isoform_exonic_length[analyze_genes, n]

    analysis <- data.frame(
      `cor(log(CPM), log(exon_length))` = cor(log(CPM), log(lengths)),
      `cor(log(TPM), log(exon_length))` = cor(log(TPM), log(lengths)),
      `cor(SCT90%, log(length))` = cor(SCT, log(lengths)),
      check.names = FALSE
    )
    writeDataTable(
      wb,
      title,
      analysis,
      startRow = start_row, startCol = 1
    )
    start_row <- start_row + 3
  }
}

write_regression_table <- function(
    wb, title, header, baseMeanCPM, apeglm, coef, flybase, table_style) {
  addWorksheet(wb, title)
  writeData(wb, title,
            matrix(
              rep(c('CPM(GSC&CySC)', header), c(1,3)),
              nrow=1
            ),
            colNames = F,
            startCol = 3, startRow = 1)
  data = tibble(
    symbol = names(flybase),
    flybase = flybase,
    baseMeanCPM = baseMeanCPM,
    p95minus = apeglm$map[,coef] - qnorm(0.975) * apeglm$sd[,coef],
    map = apeglm$map[,coef],
    p95plus = apeglm$map[,coef] + qnorm(0.975) * apeglm$sd[,coef],
    fsr = apeglm$fsr[,1] %>% signif(digits=2),
    padj = if('mle.test' %in% names(apeglm)) apeglm$mle.test$adj_pval[match(names(flybase), apeglm$mle.test$name)] %>% signif(digits=2) else 0.01
  )
  data[, 4:6] = round(as.matrix(data[, 4:6]) / log(2), digits=2)
  data = data %>% arrange(
    desc(pmax(p95minus, 0)),
    desc(pmin(p95plus, 0)),
    desc(map)
  )
  colnames(data)[4:6] = c(
    'l2FC(-)95%', 'l2FC', 'l2FC(+)95%'
  )
  writeDataTable(wb, title,
                 data,
                 startCol = 1, startRow = 2,
                 withFilter = T,
                 tableStyle = table_style)
}

write_batch_idents <- function(wb, title, batch_data) {
  addWorksheet(wb, title)
  # Add percent style for 5 of the columns.
  addStyle(wb, title, createStyle(numFmt = "0%"), rows = 1:100, cols = 2:6, gridExpand=TRUE)
  batch_data$batch <- batch_data$batch %>% factor(sce.data$batch)
  batch_data$ident <- batch_data$ident %>% factor(c(sce.clusters$cluster, "doublet"))
  tbl <- with(batch_data, table(batch, ident)) %>%
    matrix(nrow = nrow(.), dimnames = dimnames(.)) %>%
    as.data.frame
  ncells <- with(
    tbl,
    germline + somatic + spermatocyte + somaticprecursor + muscle
  )
  nice_experiment_name <- str_extract(sce.data$tenx_path, "/([^/]+)/", group=1)
  data_table <- tibble(
    Sample = nice_experiment_name,
    Germline = tbl$germline / ncells,
    Somatic = tbl$somatic / ncells,
    `Germline Differentiated` = tbl$spermatocyte / ncells,
    `Somatic Differentiated` = tbl$somaticprecursor / ncells,
    Muscle = tbl$muscle / ncells,
    Cells = ncells,
    `Doublets (Removed)` = tbl$doublet
  )
  data_table[2:6] <- data_table[2:6] %>% round(4)
  writeDataTable(
    wb,
    title,
    data_table,
    startCol = 1, startRow = 1
  )
}

write_excel_tables_list_percentages <- function(lst, output_path) {
  wb <- createWorkbook()
  for (n in names(lst)) {
    addWorksheet(wb, n)
    addStyle(wb, n, createStyle(numFmt = "0%"), rows = seq(1+nrow(lst[[n]])), cols = seq(ncol(lst[[n]])), gridExpand=TRUE)
    writeDataTable(
      wb,
      n,
      lst[[n]]
    )
  }
  saveWorkbook(wb, output_path, overwrite = TRUE)
}

publish_heatmap_named_cuts <- function(named_dendros, assay_data_sc, target_path) {
  wb = createWorkbook()
  sheet <- "Sheet1"
  addWorksheet(wb, sheet)
  assay_data_sc <- assay_data_sc %>% read.csv(row.names = 1) %>% rownames_to_column

  startCol <- 1
  for (n in names(named_dendros)) {
    writeData(wb, sheet, matrix(str_glue("Group {n}")), startCol, 1, colNames=F)
    writeDataTable(
      wb,
      sheet,
      left_join(
        tibble(rowname = labels(named_dendros[[n]])),
        assay_data_sc,
        "rowname"
      ) %>%
        reframe(`gene name` = rowname, flybase),
      startCol,
      2
    )

    startCol <- startCol + 3
  }

  dir.create(dirname(target_path), showW=F, rec=T)
  saveWorkbook(wb, target_path, overwrite=T)
  target_path
}