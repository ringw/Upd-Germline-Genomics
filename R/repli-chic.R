# Analysis with replication timing on x-axis, ChIC tracks stacked.

analyze_repli_chic <- function(
    repli.hdf5,
    chic_path,
    gtf_path,
    chic_driver,
    output_path,
    num_sorted_bases=10000
) {
  dir.create(output_path, showW = F)

  chic_marks = c('H3K4','H3K27','H3K9')
  chic_coverage <- chic_marks %>%
    sapply(
      \(n) {
        import(paste0(chic_path, '/', chic_driver, '_', n, '.q5.bw'), 'bw') %>%
          coverage(weight='score')
      },
      simplify=F
    )
  chic_coverage
  proj_matrices <- repli.hdf5 %>%
    HDF5ArraySeed(name = '/rank_split') %>%
    DelayedArray %>%
    repli_rank_projection(num_sorted_bases=num_sorted_bases)
  chic_enr <- chic_coverage %>%
    sapply(
      \(rle) apply_repli_bins(proj_matrices, rle[names(chr.lengths)]) %>% rowSums
    ) %>%
    log %>%
    `/`(log(2))
  names(dimnames(chic_enr)) <- c('replication', 'mark')

  timing_heatmap <- rasterise(
    ggplot(
      melt(chic_enr, value.name='log2(mark/input)')
      %>% within(
        mark <- mark %>% factor(c('Timing Ratio','H3K4','H3K27','H3K9'))
      )
      %>% rbind(
        data.frame(
          mark=factor('Timing Ratio', c('Timing Ratio','H3K4','H3K27','H3K9')),
          replication=seq(nrow(chic_enr)),
          `log2(mark/input)`=NA,
          check.names=F
        )
      ),
      aes(replication,mark,fill=`log2(mark/input)`)
    )
    +
    geom_tile()
    + scale_fill_viridis_c(option='magma', limits=c(-1,1), oob=squish)
    + scale_y_discrete(limits=rev, expand=c(0,0))
    + scale_x_continuous(expand=c(0,0))
    + new_scale_fill()
    + geom_tile(aes(replication,mark, fill=replication), \(df) df %>% subset(mark == 'Timing Ratio'))
    + scale_fill_gradient(limits=c(1,nrow(chic_enr)), breaks=c(0.05,0.95)*nrow(chic_enr), labels=c('Late','Early'))
  )
  print(timing_heatmap)
  ggsave(
    paste0(output_path, '/', 'Timing-BP-Heatmap.png'),
    timing_heatmap,
    width=8,
    height=6,
    dpi=120
  )

  mid <- 3900
  veryearly <- 12450
  timing_heatmap_binname <- cut(
    seq(nrow(chic_enr)),
    c(0,mid,veryearly,Inf)
  )
  levels(timing_heatmap_binname) = c('late', 'mid', 'veryearly')
  bin_colors <- list(late=hcl(247,23,29), mid=hcl(235,50,50), veryearly=hcl(226,74,74))
  timing_heatmap_bin <- rasterise(
    ggplot(
      melt(chic_enr, value.name='log2(mark/input)')
      %>% within(
        mark <- mark %>% factor(c('Timing Ratio','Cls(threshold)','H3K4','H3K27','H3K9'))
      )
      %>% rbind(
        data.frame(
          mark=factor('Timing Ratio', c('Timing Ratio','Cls(threshold)','H3K4','H3K27','H3K9')),
          replication=seq(nrow(chic_enr)),
          `log2(mark/input)`=NA,
          check.names=F
        )) %>%
      rbind(
        data.frame(
          mark=factor('Cls(threshold)', c('Timing Ratio','Cls(threshold)','H3K4','H3K27','H3K9')),
          replication=seq(nrow(chic_enr)),
          `log2(mark/input)`=NA,
          check.names=F
        )
      ),
      aes(replication,mark)
    )
    +
    geom_tile(aes(fill=`log2(mark/input)`))
    + scale_fill_viridis_c(option='magma', limits=c(-1,1), oob=squish)
    + scale_y_discrete(limits=rev, expand=c(0,0))
    + scale_x_continuous(expand=c(0,0))
    + new_scale_fill()
    + geom_tile(aes(replication,mark, fill=replication), \(df) df %>% subset(mark == 'Timing Ratio'))
    + scale_fill_gradient(limits=c(1,nrow(chic_enr)), breaks=c(0.05,0.95)*nrow(chic_enr), labels=c('Late','Early'))
    + new_scale_fill()
    + geom_tile(
        aes(replication,mark, fill=class),
        data.frame(
            mark='Cls(threshold)',
            replication=seq(nrow(chic_enr)),
            class=timing_heatmap_binname,
          `log2(mark/input)`=NA,
          check.names=F
        )
    )
    + scale_fill_manual(
        values=unlist(bin_colors)
    )
  )
  print(timing_heatmap_bin)
  ggsave(
    paste0(output_path, '/', 'Timing-BP-Heatmap-Eyeball.png'),
    timing_heatmap_bin,
    width=8,
    height=6,
    dpi=120
  )

  repli_data = repli.hdf5 %>%
    HDF5ArraySeed(name = '/array') %>%
    DelayedArray
  repli_rank <- rank(repli_data[,"score"], na.last="keep", ties.method="first")

  gtf <- read.table(
    gtf_path,
    sep = '\t',
    col.names = c('chr', 'source', 'type', 'start', 'end', 'sc', 'strand', 'fr', 'annotation'),
    header = F,
    quote = ''
  ) %>% subset(grepl('RNA', type)) # include "transcript" types in the reference
  chr.starts = setNames(
    c(0, cumsum(chr.lengths[-length(chr.lengths)])),
    names(chr.lengths)
  )
  tss.data <- gtf %>%
    with(
      data.frame(
        chr = chr,
        start = ifelse(strand == '+', start, end)
      )
    ) %>% within(
        tss_index <- chr.starts[chr] + start
    )
  tss.data$repli_rank = repli_rank[tss.data$tss] / num_sorted_bases
  for (mark in c('H3K4','H3K27'))
    tss.data[, mark] = mapply(
        \(chr, start) chic_coverage[[mark]][[chr]][start] %>% as.numeric,
        tss.data$chr,
        tss.data$start
    ) %>%
    log %>%
    `/`(log(2))
  tss.data = tss.data %>%
    subset(!is.na(repli_rank)) %>%
    split(cut(.$repli_rank, c(0, mid, veryearly, Inf)))
  names(tss.data) = c('late', 'mid', 'veryearly')

  uniform.data <- data.frame(sce.index = sample(length(repli_rank), 100000)) %>%
    subset(!is.na(repli_rank[.$sce.index])) %>%
    cbind(
      repli_rank = repli_rank[.$sce.index] / num_sorted_bases,
      chr = names(chr.starts)[findInterval(.$sce.index, chr.starts)],
      start = .$sce.index - chr.starts[findInterval(.$sce.index, chr.starts)]
    )
  for (mark in c('H3K4','H3K27'))
    uniform.data[, mark] = mapply(
        \(chr, start) chic_coverage[[mark]][[chr]][start] %>% as.numeric,
        uniform.data$chr,
        uniform.data$start
    ) %>%
    log %>%
    `/`(log(2))
  uniform.data = uniform.data %>%
    subset(!is.na(repli_rank)) %>%
    split(cut(.$repli_rank, c(0, mid, veryearly, Inf)))
  names(uniform.data) = c('late', 'mid', 'veryearly')

  uniform.data$veryearly %>% ggplot(aes(H3K4, H3K27)) + stat_density_2d(aes(fill = ..level..), geom='polygon') + scale_fill_viridis_c(option='magma') + theme_bw() + theme(panel.background = element_rect(fill = NA), panel.ontop=T)
  tss.data$veryearly %>% ggplot(aes(H3K4, H3K27)) + stat_density_2d(aes(fill = ..level..), geom='polygon') + scale_fill_viridis_c(option='magma') + theme_bw() + theme(panel.background = element_rect(fill = NA), panel.ontop=T)

  contour_plot <- bind_rows(
    list(TSS=bind_rows(tss.data, .id='repli_rank'),
         allBP=bind_rows(uniform.data, .id='repli_rank')),
    .id='selection'
  ) %>% rev %>%
    ggplot(aes(H3K4, H3K27)) + facet_grid(
        repli_rank ~ selection
    ) + stat_density_2d(
        aes(fill = ..level..), geom='polygon'
    ) + scale_fill_viridis_c(option='magma') + scale_x_continuous(
        breaks=c(-1,0,1), minor_breaks=c(-1,-0.5,0,0.5,1)
    ) + scale_y_continuous(
        breaks=c(-1,0,1), minor_breaks=c(-1,-0.5,0,0.5,1)
    ) + coord_cartesian(
        x=c(-1,1), y=c(-1,1)
    ) + theme_bw() + theme(
        panel.background = element_rect(fill = NA), panel.ontop=T
    ) + labs(
        fill = 'P(K4&K27)',
        x = 'log2(H3K4/input)',
        y = 'log2(H3K27/input)'
    )
    print(contour_plot)
    ggsave(
    paste0(output_path, '/', 'Timing-Cut-Contour-K4-K27.png'),
    contour_plot,
    width=8,
    height=6,
    dpi=120
    )
    contour_plot.veryearly <- contour_plot
    contour_plot.veryearly$data <- contour_plot.veryearly$data %>% subset(repli_rank == 'veryearly')
    print(contour_plot.veryearly)
    for (i in seq_along(names(uniform.data))) {
        contour_name <- names(uniform.data)[i]
        contour_plot.subset <- contour_plot
        contour_plot.subset$data <- contour_plot.subset$data %>% subset(repli_rank == contour_name)
    ggsave(
    paste0(output_path, '/', 'Timing-Cut-Contour-', i, '-', contour_name, '-K4-K27.png'),
    contour_plot.subset + labs(
        title = contour_name
    ) + theme(
        title = element_text(color = bin_colors[[contour_name]])
    ),
    width=8,
    height=2,
    dpi=120
    )
    }
}
