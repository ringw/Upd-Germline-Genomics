excel_tables = list(
  deepgreen='TableStyleMedium18',
  cyan='TableStyleMedium6',
  orange='TableStyleMedium14',
  red='TableStyleMedium3',
  mono='TableStyleMedium1'
)

# Write abundance table (SupplementalTable1) and regression (SupplementalTable2) table
publish_excel_results <- function(
  Upd_regression_somatic,
  Upd_cpm,
  metafeatures,
  target_path
) {
  wb = createWorkbook()

  metafeatures = metafeatures %>% read.csv(row.names=1)
  write_abundance_table(
    wb,
    'Gene Quantification',
    metafeatures$flybase,
    Upd_cpm,
    excel_tables$mono
  )

  baseMeanCPM <- exp(Upd_regression_somatic$map[, "(Intercept)"]) %>% `*`(
    1000 * 1000 / sum(., na.rm=T)
  ) %>% signif(digits=2)
  flybase <- metafeatures %>%
    rownames_to_column %>%
    subset(pass.min.cells) %>%
    pull(flybase, rowname)
  write_regression_table(wb, 'GSC&CySC Regression', 'log(CySC/GSC)', baseMeanCPM, Upd_regression_somatic, 2, flybase, excel_tables$cyan)

  dir.create('scRNA-seq-Regression', showW=F)
  saveWorkbook(wb, target_path, overwrite = T)
  target_path
}

write_abundance_table <- function(wb, title, flybase, Upd_cpm, table_style) {
  rename = c(germline='GSC', somatic="CySC", spermatocyte="tid", somaticprecursor="SP", muscle="muscle")
  data = data.frame(
    symbol = rownames(Upd_cpm),
    flybase = flybase
  )

  for (n in names(rename)) {
    name = rename[n]
    data[str_glue("CPM_{name}")] <- Upd_cpm[,n]
  }
  addWorksheet(wb, title)
  writeDataTable(wb, title,
                 data,
                 startCol = 1, startRow = 1,
                 withFilter = T,
                 tableStyle = table_style)
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
    s.value = apeglm$svalue[
      rownames(apeglm$map) %>% match(rownames(apeglm$svalue)),
      1
    ] %>%
      signif(digits=2),
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
