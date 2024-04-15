excel_tables = list(
  deepgreen='TableStyleMedium18',
  cyan='TableStyleMedium6',
  orange='TableStyleMedium14',
  red='TableStyleMedium3',
  mono='TableStyleMedium1'
)

publish_excel_results <- function(
  Upd_regression_somatic, Upd_regression_tid,
  Upd_regression_sompre, Upd_regression_mscl,
  Upd_cpm_transcripts, Upd_cpm, Upd_fpkm, sctransform_quantile,
  supplemental_bulk_fpkm, metafeatures, gtf_path, target_path
) {
  wb = createWorkbook()

  metafeatures = metafeatures %>% read.csv(row.names=1)
  write_abundance_table(
    wb,
    'Gene Quantification',
    metafeatures$flybase,
    Upd_cpm,
    Upd_fpkm,
    sctransform_quantile,
    supplemental_bulk_fpkm,
    excel_tables$mono
  )

  write_quant_stats(
    wb, "Gene Quantification Summary", Upd_cpm_transcripts, Upd_cpm, Upd_fpkm,
    Upd_regression_somatic, sctransform_quantile, gtf_path, metafeatures
  )

  baseMeanCPM <- exp(Upd_regression_somatic$map[, "(Intercept)"]) %>% `*`(
    1000 * 1000 / sum(.)
  ) %>% signif(digits=2)
  flybase = metafeatures$flybase[match(rownames(Upd_regression_somatic$map), rownames(metafeatures))]
  write_regression_table(wb, 'GSC&CySC Regression', 'log(CySC/GSC)', baseMeanCPM, Upd_regression_somatic, 2, flybase, excel_tables$cyan)
  write_regression_table(wb, 'S-Cyte&Tid-Like Regression', 'log(cyte&tid/GSC&CySC)', baseMeanCPM, Upd_regression_tid, 3, flybase, excel_tables$orange)
  write_regression_table(wb, 'Somatic Precursor Regression', 'log(SomPre/GSC&CySC)', baseMeanCPM, Upd_regression_sompre, 4, flybase, 'TableStyleMedium13')
  write_regression_table(wb, 'Muscle Regression', 'log(mscl/GSC&CySC)', baseMeanCPM, Upd_regression_mscl, 5, flybase, excel_tables$red)
  dir.create('scRNA-seq-Regression', showW=F)
  saveWorkbook(wb, target_path, overwrite = T)
  target_path
}

write_abundance_table <- function(
    wb, title, flybase, Upd_cpm, Upd_fpkm, sctransform_quantile, supplemental_bulk_fpkm, table_style) {
  rename = c(germline='GSC', somatic="CySC", spermatocyte="tid", somaticprecursor="SC", muscle="muscle")
  data = data.frame(
    symbol = rownames(Upd_fpkm),
    flybase = flybase
  )

  for (n in names(rename)) {
    name = rename[n]
    data = cbind(
      data,
      setNames(
        list(Upd_cpm[,n], Upd_fpkm[,n]),
        paste(c('CPM','FPKM'), name, sep='_')
      )
    )
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
        symbol = rownames(supplemental_bulk_fpkm),
        pctAllCells = supplemental_bulk_fpkm$percent_expressed
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
  wb, title, Upd_cpm_transcripts, Upd_cpm, Upd_fpkm, Upd_regression_somatic,
  sctransform_quantile, gtf_path, metafeatures
) {
  addWorksheet(wb, title)

  n_genes <- tribble(
    ~ cluster, ~ CPM, ~ FPKM, ~ SCT, ~ L2FC,
    "GSC",
    sum((Upd_cpm[, "germline"] > 0) %>% replace(is.na(.), FALSE)),
    sum((Upd_fpkm[, "germline"] > 0) %>% replace(is.na(.), FALSE)),
    sum(is.finite(sctransform_quantile[["germline"]][, "90%"])),
    sum(is.finite(Upd_regression_somatic$map[, 2])),
    "CySC",
    sum((Upd_cpm[, "somatic"] > 0) %>% replace(is.na(.), FALSE)),
    sum((Upd_fpkm[, "somatic"] > 0) %>% replace(is.na(.), FALSE)),
    sum(is.finite(sctransform_quantile[["somatic"]][, "90%"])),
    sum(is.finite(Upd_regression_somatic$map[, 2]))
  )

  writeData(wb, title, paste0("# Genes Quantified by Method (Genes in Reference: ", nrow(Upd_fpkm), ")"), colNames = F, startRow = 1, startCol = 1)
  writeDataTable(wb, title,
                 n_genes, withFilter = FALSE,
                 startCol = 1, startRow = 2)
  
  start_row <- 6
  group_names <- c(germline="GSC", somatic="CySC")
  analyze_genes <- rownames(Upd_fpkm) %>%
    setdiff(rownames(Upd_fpkm) %>% subset(is.na(Upd_fpkm[, "germline"]))) %>%
    setdiff(rownames(Upd_fpkm) %>% subset(is.na(Upd_fpkm[, "somatic"]))) %>%
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
      paste0(group_names[n], ": Transcript Length as a Confounding Factor"),
      colNames = F,
      startRow = start_row, startCol = 1
    )
    start_row <- start_row + 1
    transcripts <- attr(Upd_fpkm, "transcript_id")[analyze_genes, n]
    CPM_naive = (
      Upd_cpm_transcripts[, n] %>% split(Upd_cpm_transcripts$gene_id) %>% sapply(max)
    )[metafeatures[analyze_genes, "flybase"]]
    CPM = Upd_cpm[analyze_genes, n]
    FPKM = Upd_fpkm[analyze_genes, n]
    SCT = sctransform_quantile[[n]][analyze_genes, "90%"]
    lengths <- gtf$length[match(transcripts, gtf$tx)]

    analysis <- data.frame(
      `cor(log(CPM[naive]), log(length))` = cor(log(CPM_naive) %>% subset(CPM != 0), log(lengths) %>% subset(CPM != 0)),
      `cor(log(CPM), log(length))` = cor(log(CPM) %>% subset(CPM != 0), log(lengths) %>% subset(CPM != 0)),
      `cor(log(FPKM), log(length))` = cor(log(FPKM) %>% subset(CPM != 0), log(lengths) %>% subset(CPM != 0)),
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
  data = data.frame(
    symbol = rownames(apeglm$map),
    flybase = flybase,
    baseMeanCPM = baseMeanCPM,
    p95minus = apeglm$map[,coef] - qnorm(0.975) * apeglm$sd[,coef],
    map = apeglm$map[,coef],
    p95plus = apeglm$map[,coef] + qnorm(0.975) * apeglm$sd[,coef],
    fsr = apeglm$fsr[,1] %>% signif(digits=2),
    padj = if('mle.test' %in% names(apeglm)) apeglm$mle.test$adj_pval[match(rownames(apeglm$map), apeglm$mle.test$name)] %>% signif(digits=2) else 0.01
  )
  data[, 4:6] = round(data[, 4:6] / log(2), digits=2)
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

write_excel_tables_list <- function(lst, output_path) {
  wb <- createWorkbook()
  for (n in names(lst)) {
    addWorksheet(wb, n)
    writeDataTable(
      wb,
      n,
      lst[[n]],
      rowNames = TRUE
    )
  }
  saveWorkbook(wb, output_path, overwrite = TRUE)
}