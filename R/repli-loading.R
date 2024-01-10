library(dplyr)
library(glmGamPoi)
library(Matrix)
library(purrr)
library(SingleCellExperiment)
library(tibble)
library(zoo)
# source('C:/Users/dringwa1/OneDrive - Johns Hopkins/Documents/scTJ/Rollutils.R')

REPLI_WINDOW_SIZE = 1000

colData = rbind(
  c(filename='Tj_Early-tech1_S6_L001', ident='Tj', frac=1, replicate=1),
  c(filename='Tj_Early-tech1_S6_L002', ident='Tj', frac=1, replicate=1),
  c(filename='Tj_Early-tech2_S9_L001', ident='Tj', frac=1, replicate=2),
  c(filename='Tj_Early-tech2_S9_L002', ident='Tj', frac=1, replicate=2),
  c(filename='Tj_Early-Mid-tech_1_S12_L001', ident='Tj', frac=2, replicate=1),
  c(filename='Tj_Early-Mid-tech_1_S12_L002', ident='Tj', frac=2, replicate=1),
  c(filename='Tj_Early-Mid-tech_2_S3_L001', ident='Tj', frac=2, replicate=2),
  c(filename='Tj_Early-Mid-tech_2_S3_L002', ident='Tj', frac=2, replicate=2),
  c(filename='Tj_Mid-Late-tech_1_S4_L001', ident='Tj', frac=3, replicate=1),
  c(filename='Tj_Mid-Late-tech_1_S4_L002', ident='Tj', frac=3, replicate=1),
  c(filename='Tj_Mid-Late-tech_2_S7_L001', ident='Tj', frac=3, replicate=2),
  c(filename='Tj_Mid-Late-tech_2_S7_L002', ident='Tj', frac=3, replicate=2),
  c(filename='Tj_Late-tech_1_S5_L001', ident='Tj', frac=4, replicate=1),
  c(filename='Tj_Late-tech_1_S5_L002', ident='Tj', frac=4, replicate=1),
  c(filename='Tj_Late-tech_2_S2_L001', ident='Tj', frac=4, replicate=2),
  c(filename='Tj_Late-tech_2_S2_L002', ident='Tj', frac=4, replicate=2),
  c(filename='GSC_Early_S10_L001', ident='GSC', frac=1, replicate=1),
  c(filename='GSC_Early_S10_L002', ident='GSC', frac=1, replicate=1),
  c(filename='GSC_Early-Mid_S8_L001', ident='GSC', frac=2, replicate=1),
  c(filename='GSC_Early-Mid_S8_L002', ident='GSC', frac=2, replicate=1),
  c(filename='GSC_Mid-Late_S11_L001', ident='GSC', frac=3, replicate=1),
  c(filename='GSC_Mid-Late_S11_L002', ident='GSC', frac=3, replicate=1),
  c(filename='GSC_Late_S1_L001', ident='GSC', frac=4, replicate=1),
  c(filename='GSC_Late_S1_L002', ident='GSC', frac=4, replicate=1)
) %>% as.data.frame %>% column_to_rownames('filename')
colData$frac = colData$frac %>% factor

colData_concat = rbind(
  c(filename='Tj_Early', ident='Tj', frac=1),
  c(filename='Tj_Early-Mid', ident='Tj', frac=2),
  c(filename='Tj_Mid-Late', ident='Tj', frac=3),
  c(filename='Tj_Late', ident='Tj', frac=4)
) %>% as.data.frame %>% column_to_rownames('filename')

colData_kc = rbind(
  c(filename='SRR586004', ident='Kc167', frac=1),
  c(filename='SRR586005', ident='Kc167', frac=2),
  c(filename='SRR586006', ident='Kc167', frac=3),
  c(filename='SRR586007', ident='Kc167', frac=4),
  c(filename='SRR586008', ident='Kc167', frac=1),
  c(filename='SRR586009', ident='Kc167', frac=2),
  c(filename='SRR586010', ident='Kc167', frac=3),
  c(filename='SRR586011', ident='Kc167', frac=4)
) %>% as.data.frame %>% column_to_rownames('filename')

chr_correction = function(data) {
  bin_sizes = data$end - data$start
  bin_size = min(bin_sizes[-length(bin_sizes)])
  data.frame(
    chr = data[1, 'chr'],
    start = seq(0, max(data$end), by=bin_size),
    end = c(seq(bin_size, max(data$end), by=bin_size), max(data$end)),
    value = rep(data$value, ceiling((data$end - data$start) / bin_size))
  )
}
bedgraph_correction = function(bedgraph) {
  bedgraph %>% as.data.frame %>% subset(
    chr %in% c('2L','2R','3L','3R','4','X','Y')
  ) %>% split(.$chr) %>% sapply(chr_correction, simplify=F) %>% do.call(rbind,.)
}

load_repli_colData <- function(repli_dir, colData, suffix) {
  repli = mapply(
    \(f) read.table(
      paste0(repli_dir, '/', f, suffix),
      header=F,
      col.names=c('chr', 'start', 'end', 'value')
    ) %>% bedgraph_correction,
    rownames(colData),
    SIMPLIFY=F
  )
  for (n in names(repli)) colnames(repli[[n]])[4] = n
  dim(repli[[1]])
  repli = purrr::reduce(
    repli,
    \(a,b) full_join(a,b, c('chr','start','end'))
  )
  dim(repli)
  repli = SingleCellExperiment(
    list(rpkm=repli[, -(1:3)] %>% as.matrix %>% replace(!is.finite(.), 0)),
    colData=colData,
    rowData=repli[, 1:3] %>% as.data.frame
  )
  repli
}

load_repli <- function(repli_dir, suffix='.q20.10kb.bedgraph') {
  repli = mapply(
    \(f) read.table(
      paste0(repli_dir, '/', f, suffix),
      header=F,
      col.names=c('chr', 'start', 'end', 'value')
    ) %>% bedgraph_correction,
    rownames(colData),
    SIMPLIFY=F
  )
  for (n in names(repli)) colnames(repli[[n]])[4] = n
  dim(repli[[1]])
  repli = purrr::reduce(
    repli,
    \(a,b) full_join(a,b, c('chr','start','end'))
  )
  dim(repli)
  repli = SingleCellExperiment(
    list(rpkm=repli[, -(1:3)] %>% as.matrix %>% replace(!is.finite(.), 0)),
    colData=colData,
    rowData=repli[, 1:3] %>% as.data.frame
  )
  rowData(repli)$chr = rowData(repli)$chr %>% factor(names(chr.lengths))
  repli
}

pseudobulk_repli <- function(repli) {
  repli = pseudobulk(repli, group_by=vars(ident,frac))
  # Normalize chromosome content summing up all conditions (fractions). The
  # assays will each have a sum of approximately 1MM / (window size kb) (RPKM
  # normalizes each observation, does not normalize the total). The scaled
  # assays will have a sum of exactly this quantity.
  chr_norm_factor = as.data.frame(assay(repli)) %>% split(rowData(repli)$chr) %>% sapply(\(df) mean(as.matrix(df)))
  chr_vec = (1/chr_norm_factor)[rowData(repli)$chr]
  assay(repli) = assay(repli) * chr_vec
  assay(repli, withDimnames=F) = assay(repli) %*% diag(
    1000 * 1000 * 1000 / REPLI_WINDOW_SIZE / colSums(assay(repli))
  )
  # Now we cap every track at the 99.9%ile.
  assay(repli) = pmin(
    assay(repli),
    rep(
      colQuantiles(assay(repli), probs=0.999),
      each=nrow(repli)
    )
  )
  colData(repli)
  # colnames(repli) = paste(
  #   rep(c('GSC','Tj'), each=4),
  #   rep(c('Early','Early-Mid','Mid-Late','Late'), 2),
  #   sep='_'
  # )
  colnames(repli) = paste(
    'Tj',
    c('Early','Early-Mid','Mid-Late','Late'),
    sep='_'
  )

  # Compute the pct after pseudobulk.
  assay(repli,'pct') = (
    assay(repli,'rpkm') / rowSums(assay(repli,'rpkm'))
  )
  assay(repli,'pct') = (
    assay(repli,'pct')
    * ifelse(
      rowSums(assay(repli,'rpkm')[, c('Tj_Early','Tj_Late')]) < 10,
      NA,
      1
    )
  )

  repli
}

write_repli_bedgraph_out <- function(repli, out_dir) {
  options(scipen=100)
  for (name in colnames(repli)) {
    write.table(
      cbind(
        rowData(repli),
        assay(repli,'rpkm')[, name]
      ),
      paste0(out_dir, '/RPKM_', name, '_window_5kb.bedgraph'),
      sep='\t',
      quote=F,
      row.names=F,
      col.names=F
    )
  }
  options(scipen=0)
}

repli_pct <- function(repli, maxgap=Inf) {
  rpkm_sum = rowSums(assay(repli)[,c('Tj_Early','Tj_Late')])
  pct_data = data.frame(pct = assay(repli)[,'Tj_Early'] / rpkm_sum)
  pct_data$pct = pct_data$pct %>% replace(rpkm_sum < 10, NA)
  pct_data = na.approx(pct_data, na.rm=F, maxgap=maxgap) %>% as.data.frame
  pct_data$pct
}

write_repli_pct <- function(repli, out_dir) {
  options(scipen=100)
  pct_data = data.frame(pct=repli_pct(repli))
  write.table(
    cbind(
      rowData(repli),
      pct_data
    ) %>% subset(!is.na(pct_data$pct)),
    paste0(out_dir, '/Early_Score_Tj.bedgraph'),
    sep='\t', quote=F, row.names=F, col.names=F
  )
}
