chic_colData = rbind(
  c(name='Nos_H3K4', driver='Nos', mod='K4'),
  c(name='Nos_H3K27', driver='Nos', mod='K27'),
  c(name='Nos_H3K9', driver='Nos', mod='K9'),
  c(name='tj_H3K4', driver='tj', mod='K4'),
  c(name='tj_H3K27', driver='tj', mod='K27'),
  c(name='tj_H3K9', driver='tj', mod='K9')
) %>% as.data.frame %>% column_to_rownames('name')
chic_colData$driver = chic_colData$driver %>% factor
chic_colData$mod = chic_colData$mod %>% factor(c('K4','K27','K9'))

window_chic <- function(dirname, window_size=1000) {
  if (F)
  chic_ranges <- sapply(
    rownames(chic_colData),
    \(name) {
      filename = paste0(dirname, '/', name, '.q5.bw')
      bigwig = import(filename, 'bw')
      bigwig_score = coverage(bigwig, weight='score')
      chr.tiles = tileGenome(chr.lengths, tilewidth=1000, cut.last=T)
      seqlevels(chr.tiles, pruning.mode='coarse') = names(bigwig_score)
      granges = binnedAverage(chr.tiles, bigwig_score, 'score')
      data.frame(
        chr = granges@seqnames,
        start = granges@ranges@start-1,
        end = granges@ranges@start + granges@ranges@width - 1,
        score = granges$score
      )
    },
    simplify=F
  )
  # Median computation
  chic_ranges <- sapply(
    rownames(chic_colData),
    \(name) {
      filename = paste0(dirname, '/', name, '.q5.bw')
      bigwig = import(filename, 'bw')
      bigwig_score = coverage(bigwig, weight='score')
      chr.tiles = tileGenome(chr.lengths, tilewidth=1000, cut.last=T)
      seqlevels(chr.tiles, pruning.mode='coarse') = names(bigwig_score)
      df = data.frame(
        chr = chr.tiles@seqnames,
        start = chr.tiles@ranges@start-1,
        end = chr.tiles@ranges@start + chr.tiles@ranges@width - 1
      )
      df$score = NA
      for (i in seq(nrow(df))) {
        df[i, 'score'] = median(
          bigwig_score[[df[i, 'chr']]][
            (df[i, 'start']+1):(df[i, 'end'])
          ]
        )
      }
      df
    },
    simplify=F
  )
  for (n in names(chic_ranges)) colnames(chic_ranges[[n]])[4] = n
  chic_ranges = chic_ranges %>% purrr::reduce(
    \(a,b) full_join(a,b, c('chr','start','end'))
  )
  chic_ranges$chr = chic_ranges$chr %>% factor(names(chr.lengths))
  SingleCellExperiment(
    list(ratio=chic_ranges[, -(1:3)] %>% as.matrix),
    colData=chic_colData,
    rowData=chic_ranges[, 1:3] %>% as.data.frame
  )
}

make_tj_bp_early_score <- function(bedgraph_in) {
tj_colnames = c(
    'Tj_Early-tech1_S6_L001',
    'Tj_Early-tech1_S6_L002',
    'Tj_Early-tech2_S9_L001',
    'Tj_Early-tech2_S9_L002',
    'Tj_Late-tech_1_S5_L001',
    'Tj_Late-tech_1_S5_L002',
    'Tj_Late-tech_2_S2_L001',
    'Tj_Late-tech_2_S2_L002'
)

tj_data = matrix(
    # i=1, j=1, x=0,
    # dims=c(sum(chr.lengths), length(tj_colnames)),
    nrow=sum(chr.lengths),
    ncol=length(tj_colnames),
    dimnames=list(NULL, tj_colnames)
)

for (n in tj_colnames) {
    bigwig = import(
        paste0('repli/bedgraph.in/', n, '.q20.near10.coverage10k.bw'),
        'bw'
    )
    bigwig_score = coverage(bigwig, weight='score')
    start_index = 1
    for (chr_name in names(chr.lengths)) {
        chr_data = bigwig_score[[chr_name]]
        tj_data[
            start_index:(start_index + length(chr_data) - 1),
            n
        ] <- as.numeric(chr_data)
        start_index <- start_index + length(chr_data)
    }
}

tj_data = tj_data / (rowSums2(tj_data) %>% replace(. < 4, NA))
tj_data = tj_data %*% c(1,1,1,1,0,0,0,0)

# Now rank each base pair, and create a linear transformation from the
# chromosomes to the bucketed base pair statistics.
tj_rank = rank(tj_data, na.last='keep', ties.method='first')
bucket_enr_size = 1000
bucket_rank = cut(
  tj_rank,
  c(
    seq(
      0,
      max(tj_rank, na.rm=T)-bucket_enr_size,
      by=bucket_enr_size),
    Inf)
)
num_ranks = length(levels(bucket_rank))
proj_matrices = list()
start_index = 1
for (chr_name in names(chr.lengths)) {
    i = as.numeric(bucket_rank[
        start_index:(start_index + chr.lengths[chr_name] - 1)
    ])
    proj_matrices[[chr_name]] = sparseMatrix(
      i = i %>% subset(!is.na(i)),
      j = seq(chr.lengths[chr_name])[!is.na(i)],
      x = 1/bucket_enr_size,
      dims=c(num_ranks, chr.lengths[chr_name])
    )
    start_index <- start_index + chr.lengths[chr_name]
}
list(
  tj.pct=tj_data,
  tj.matrices=proj_matrices
)
}

apply_repli_bins <- function(matrices, rles) {
  mapply(
    \(m, v) {
      result = rep(0, nrow(m))
      block_size = 1024 * 1024
      start = 1
      end = block_size
      while (start <= length(v)) {
        end = pmin(end, length(v))
        result = result + (
          m[, start:end]
          %*%
          as.numeric(v)[start:end]
        )
        start <- end + 1
        end <- end + block_size
      }
      as.numeric(result)
    },
    matrices, rles
  )
}

plot_repli_mat_chic <- function(repli.pct.bp) {
  chic_enr <- sapply(
    c('tj_H3K4', 'tj_H3K27', 'tj_H3K9'),
    \(n) {
      bigwig = import(paste0('chic/',n,'.q5.bw'), 'bw')
      bigwig_score = coverage(bigwig, weight='score')
      apply_repli_bins(
        repli.pct.bp$tj.matrices,
        bigwig_score[names(chr.lengths)]
      ) %>% rowSums
    }
  )
  chic_enr = log(chic_enr) / log(2)
  names(dimnames(chic_enr)) = c('replication', 'mark')
  ggsave('Timing-Ratio.png',
  rasterise(
    ggplot(
      melt(chic_enr, value.name='log2(mark/input)')
      %>% within(
        mark <- mark %>% factor(c('Timing Ratio','tj_H3K4','tj_H3K27','tj_H3K9'))
      )
      %>% rbind(
        data.frame(
          mark=factor('Timing Ratio', c('Timing Ratio','tj_H3K4','tj_H3K27','tj_H3K9')),
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
  , dpi=600), width=8, height=6, dpi=120)
}