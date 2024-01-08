excel_tables = list(
  deepgreen='TableStyleMedium18',
  cyan='TableStyleMedium6',
  orange='TableStyleMedium14',
  red='TableStyleMedium3',
  mono='TableStyleMedium1'
)

publish_excel_results <- function(glm, apeglm_somatic, apeglm_spermatocyte, apeglm_muscle, metafeatures, target_path) {
  wb = createWorkbook()
  metafeatures = metafeatures %>% read.csv(row.names=1)
  flybase = metafeatures$flybase[match(rownames(apeglm_somatic$map), rownames(metafeatures))]
  write_regression_table(wb, 'Germline and Somatic Tumor', 'log(som/germ)', glm, apeglm_somatic, flybase, excel_tables$cyan)
  write_regression_table(wb, 'Spermatocyte', 'log(spmtc/others)', glm, apeglm_spermatocyte, flybase, excel_tables$orange)
  write_regression_table(wb, 'Muscle', 'log(mscl/others)', glm, apeglm_muscle, flybase, excel_tables$red)
  flybase = metafeatures$flybase[match(rownames(glm$Beta), rownames(metafeatures))]
  tx_length = metafeatures$transcript.length[match(rownames(glm$Beta), rownames(metafeatures))]
  write_abundance_table(wb, 'CPM and FPKM', glm, flybase, tx_length, excel_tables$mono)
  dir.create('scRNA-seq-Regression', showW=F)
  saveWorkbook(wb, target_path, overwrite = T)
  target_path
}

write_regression_table <- function(
    wb, title, header, glm, apeglm, flybase, table_style) {
  addWorksheet(wb, title)
  writeData(wb, title,
            matrix(
              rep(c('avg(germ)', header), c(1,3)),
              nrow=1
            ),
            colNames = F,
            startCol = 2, startRow = 1)
  data = data.frame(
    symbol = rownames(apeglm$map),
    flybase = flybase,
    baseMean = exp(
      mean(glm$data$size_factor %>% subset(glm$data$ident == 'germline'))
      + apeglm$map[,1]
    ) %>% signif(digits=2),
    p95minus = apeglm$map[,2] - qnorm(0.975) * apeglm$sd[,2],
    map = apeglm$map[,2],
    p95plus = apeglm$map[,2] + qnorm(0.975) * apeglm$sd[,2],
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

write_abundance_table <- function(wb, title, glm, flybase, tx_length, table_style) {
  contrasts = cbind(
    germline = c(1,0,0,0,0,0,0),
    somatic = c(1,1,0,0,0,0,0),
    spermatocyte = c(1,0,1,0,0,0,0),
    muscle = c(1,0,0,1,0,0,0)
  )
  rename = c(germline='germ')
  data = data.frame(
    symbol = rownames(glm$Beta),
    flybase = flybase
  )
  for (n in colnames(contrasts)) {
    umi_count = glm$data$nCount_RNA %>% subset(
      as.character(glm$data$ident) == n
    ) %>% sum
    size_factor = glm$data$size_factor %>% subset(
      as.character(glm$data$ident) == n
    ) %>% sum
    cpm = exp(glm$Beta %*% contrasts[,n]) * (
      size_factor * 1000 * 1000 / umi_count
    )
    fpkm = cpm * 1000 / tx_length
    name = if (n %in% names(rename)) rename[n] else n
    data = cbind(
      data,
      setNames(list(cpm,fpkm), paste(c('CPM','FPKM'), name, sep='_'))
    )
  }
  addWorksheet(wb, title)
  writeDataTable(wb, title,
                 data,
                 startCol = 1, startRow = 1,
                 withFilter = T,
                 tableStyle = table_style)
}
