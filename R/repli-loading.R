library(dplyr)
library(glmGamPoi)
library(Matrix)
library(purrr)
library(SingleCellExperiment)
library(tibble)
library(zoo)
source('C:/Users/dringwa1/OneDrive - Johns Hopkins/Documents/scTJ/Rollutils.R')

colData = rbind(
  c(filename='GSC_Early_S10_L001', ident='GSC', frac=1, replicate=1),
  c(filename='GSC_Early_S10_L002', ident='GSC', frac=1, replicate=1),
  c(filename='GSC_Early-Mid_S8_L001', ident='GSC', frac=2, replicate=1),
  c(filename='GSC_Early-Mid_S8_L002', ident='GSC', frac=2, replicate=1),
  c(filename='GSC_Mid-Late_S11_L001', ident='GSC', frac=3, replicate=1),
  c(filename='GSC_Mid-Late_S11_L002', ident='GSC', frac=3, replicate=1),
  c(filename='GSC_Late_S1_L001', ident='GSC', frac=4, replicate=1),
  c(filename='GSC_Late_S1_L002', ident='GSC', frac=4, replicate=1),
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
  c(filename='Tj_Late-tech_2_S2_L002', ident='Tj', frac=4, replicate=2)
) %>% as.data.frame %>% column_to_rownames('filename')
suffixes = data.frame(
  suffix=c('.bedgraph', '.complete.bedgraph'),
  coef=list(p1=1 - exp(-1/10*log(10))) %>% with(c(1-p1, p1))
)

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

load_repli <- function(repli_dir) {
  repli = mapply(
    \(f, suffix) read.table(
      paste0(repli_dir, '/', f, suffix),
      header=F,
      col.names=c('chr', 'start', 'end', 'value')
    ) %>% bedgraph_correction,
    merge(rownames_to_column(colData),suffixes)$rowname,
    merge(colData,suffixes)$suffix,
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
    list(counts=repli[, -(1:3)] %>% as.matrix %>% replace(!is.finite(.), 0)),
    colData=merge(rownames_to_column(colData),suffixes),
    rowData=repli[, 1:3] %>% as.data.frame
  )
  repli
}

normalize_repli <- function(repli) {
  assay(repli) = assay(repli) %>% t %>% `*`(repli$coef) %>% t
  repli = pseudobulk(repli, group_by=vars(rowname,ident,frac,replicate))
  colnames(repli) = repli$rowname
  assay(repli, 'rpkm') = (
    assay(repli, 'counts') %>% t %>% `*`(1000 * 1000 / rowSums(.)) %>% t
  ) * ifelse(
    rowData(repli)$chr %in% c('X','Y'),
    2,
    1
  )
  repli = pseudobulk(repli, group_by=vars(ident,frac))
  colData(repli)
  colnames(repli) = paste(
    rep(c('GSC','Tj'), each=4),
    rep(c('Early','Early-Mid','Mid-Late','Late'), 2),
    sep='_'
  )
  # pseudo-log2-RPKM, tolerating sparse data
  assay(repli, 'pl2rpkm') = log1p(exp(log(assay(repli,'rpkm'))/log(2)))
  assay(repli, 'l2rpkm') = (log(assay(repli,'rpkm'))/log(2)) %>% pmax(-10)
  assay(repli, 'pct', withDimnames=F) = matrix(NA, nrow=nrow(repli), ncol=ncol(repli))
  for (ident in unique(repli$ident)) {
    repli.subset = assay(repli, 'rpkm')[, repli$ident == ident]
    repli.subset = repli.subset %>% assayRollapply(10, mean, rowData(repli))
    repli.pct = repli.subset / rowSums(repli.subset)
    assay(repli, 'pct')[, colnames(repli.pct)] = repli.pct
  }
  repli
}
