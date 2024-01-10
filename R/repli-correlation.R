repli_qc_make_figures <- function(objs, output_path) {
  dir.create(output_path, show=F)
  process_repli <- function(prefix, repli) {
    rbind(
      cor(
        assay(repli,'rpkm')[, repli$ident=='Tj'],
        method='spearman'
      ) %>% melt %>% within(samples <- 'Tj')
    ) %>% cbind(
      matrix(as.matrix(list(c('Early','Early-Mid','Mid-Late','Late')) %>% rep(2) %>% do.call(expand.grid, .)), ncol=2, dimnames=list(NULL, c('sample1', 'sample2')))
      %>% as.data.frame
      %>% within(sample1 <- factor(sample1, c('Early','Early-Mid','Mid-Late','Late')))
      %>% within(sample2 <- factor(sample2, c('Early','Early-Mid','Mid-Late','Late')))
    ) %>% ggplot(
      aes(x=sample1, y=sample2, fill=value)
    ) + facet_wrap(vars(samples)) + geom_tile() + scale_fill_viridis_c(
      option='magma', limits=c(0.25,1), begin=0, end=0.7
    ) + scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0), limits=rev) + theme(
      strip.text.x = element_text(size=12)
    )
    ggsave(
      paste0(output_path, '/', prefix, '_Spearman.png'),
      width = 2.5, height = 2, scale = 2
    )
  }
  mapply(process_repli, names(objs), objs)
  output_path
}