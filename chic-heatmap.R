heatmap_by_repli <- function(repli, chic.1kb) {
  tj_data = data.frame(
    early_rep = repli_pct(repli, maxgap=10),
    H3K4 = assay(chic.1kb)[,'tj_H3K4'],
    H3K27 = assay(chic.1kb)[,'tj_H3K27'],
    H3K9 = assay(chic.1kb)[,'tj_H3K9']
  ) %>% subset(!is.na(early_rep))
  rasterise(
    ggplot(
      melt(tj_data, id.vars='early_rep'),
      aes(rank(early_rep), variable, fill=log(value) / log(2))
    ) + geom_tile() + scale_fill_viridis_c(
      option='magma', limits=c(-2.5,2.5), oob=squish
    ) + scale_y_discrete(limits=rev, expand=c(0,0))
    + scale_x_discrete(expand=c(0,0))
    + labs(x='Replication: Late to Early', y = 'ChIC enrichment')
  )
  rasterise(
    ggplot(
      melt(tj_data, id.vars='early_rep'),
      aes(rank(early_rep), variable, fill=value)
    ) + geom_tile() + scale_fill_viridis_c(
      option='magma', limits=c(0,3), oob=squish
    ) + scale_y_discrete(limits=rev, expand=c(0,0))
    + scale_x_discrete(expand=c(0,0))
    + labs(x='Replication: Late to Early', y = 'ChIC enrichment'),
    dpi=50
  )
}