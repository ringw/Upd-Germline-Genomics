# Base pair-level coverage of replication timing. The intermediate result (the
# rank of the replication timing estimate for every base pair) is large.

na.fill.sparse <- function(vec) sparseVector(
  vec %>% subset(!is.na(.)),
  which(!is.na(vec)),
  length(vec)
)

repli_coverage_contrast <- function(
    repli_dir,
    repli_hdf5,
    ident,
    contrast,
    suffix='.q20.near10.coverage10k.bw',
    min_rpkm=200,
    num_sorted_bases=10000
) {
    names(contrast) <- levels(colData$frac)
    contrast = contrast[contrast != 0]
    rpkm_data = matrix(
        0,
        nrow=sum(chr.lengths),
        ncol=length(contrast),
        dimnames=list(NULL, names(contrast))
    )
    rpkm_data = matrix(
        0,
        nrow=sum(chr.lengths),
        ncol = length(contrast) + 3,
        dimnames=list(
            NULL,
            c(names(contrast), 'score', 'pearson', 'spearman')
        )
    )
    num_samples <- unique(
      colData[colData$ident == ident, ] %>%
        group_by(frac) %>%
        summarise(samples = length(ident)) %>%
        pull(samples)
    )
    stopifnot(length(num_samples) == 1)
    for (filename in rownames(colData)) {
        if (colData[filename,'ident'] == ident && as.character(colData[filename,'frac']) %in% names(contrast)) {
            bigwig <- BiocIO::import(
                paste0(repli_dir, '/', filename, suffix),
                'bw'
            )
            bigwig_score = coverage(bigwig, weight='score')

            fraction <- as.character(colData[filename,'frac'])
            start_index <- 1
            for (chr_name in names(chr.lengths)) {
                rpkm_data[
                    start_index:(start_index + chr.lengths[chr_name] - 1),
                    fraction
                ] <- (
                    rpkm_data[
                        start_index:(start_index + chr.lengths[chr_name] - 1),
                        fraction
                    ]
                    + as.numeric(bigwig_score[[chr_name]])
                )
                start_index <- start_index + length(bigwig_score[[chr_name]])
            }
            stopifnot(start_index == sum(chr.lengths) + 1)
        }
    }
    rpkm_data <- rpkm_data / num_samples
    repli_vec <- as.numeric(rpkm_data %*% c(contrast, 0,0,0)) %>%
      `/`(
        rowSums(rpkm_data) %>% replace(. < min_rpkm, NA)
      )
    rpkm_data[, "score"] <- repli_vec
    # Fix the Early/Late score to be between 0 and 1.
    if (isTRUE(all.equal(min(contrast), -1)) && isTRUE(all.equal(max(contrast), 1)))
      repli_vec <- 0.5 + 0.5 * repli_vec

    rpkm_data[, "pearson"] <- scale(repli_vec) %>%
      as.numeric %>%
      replace(is.na(.), 0)

    repli_rank <- rank(repli_vec, na.last="keep", ties.method="first")
    rpkm_data[, "spearman"] <- scale(repli_rank) %>%
      as.numeric %>%
      replace(is.na(.), 0)
    dir.create(dirname(repli_hdf5), showW = F)
    rpkm_data <- rpkm_data %>%
      writeHDF5Array(filepath = repli_hdf5, name = "/array", with.dimnames = T)
    bucket_rank = repli_rank %>%
      cut(
        c(
            seq(
            0,
            max(repli_rank, na.rm=T)-num_sorted_bases,
            by=num_sorted_bases),
            Inf
        )
      )
  num_ranks <- length(levels(bucket_rank))
  rank_chrs <- matrix(
    NA,
    nrow = max(chr.lengths),
    ncol = length(chr.lengths),
    dimnames = list(NULL, names(chr.lengths))
  )
  start_index <- 1
  for (chr_name in names(chr.lengths)) {
    # A factor is nonzero; replace base pair NA with 0
    rank_chrs[
      seq(chr.lengths[chr_name]),
      chr_name
    ] <- as.numeric(bucket_rank[
      start_index:(start_index + chr.lengths[chr_name] - 1)
    ]) %>%
      replace(is.na(.), 0)
    start_index <- start_index + chr.lengths[chr_name]
  }
  stopifnot(start_index == sum(chr.lengths) + 1)
  writeHDF5Array(
    rank_chrs,
    repli_hdf5,
    name = "/rank_split",
    with.dimnames = T
  )
  return(repli_hdf5)
  proj_matrices <- rank_split %>%
    mapply(
      \(v, chr.length) list(i = v %>% as.numeric) %>%
        with(
          sparseMatrix(
            i = i %>% subset(!is.na(.)),
            j = which(!is.na(i)),
            x = 1/num_sorted_bases,
            dims = c(num_ranks, chr.length),
            repr="R"
          )
        ),
      .,
      chr.lengths,
      SIMPLIFY=F
    )
}

repli_rank_projection <- function(rank_split, num_sorted_bases=10000) {
  num_ranks <- max(rank_split, na.rm=T)
  mapply(
    \(n, chr.length) sparseMatrix(
      i = rank_split[1:chr.length, n] %>% subset(. != 0),
      j = which(rank_split[1:chr.length, n] != 0),
      x = 1/num_sorted_bases,
      dims = c(num_ranks, chr.length),
      repr = "R"
    ),
    names(chr.lengths),
    chr.lengths,
    SIMPLIFY=F
  )
}

analyze_repli_quarters <- function(
  repli.hdf5,
  metafeatures_path,
  gtf_path,
  nquarters=4
) {
  rank_repli_coverage <- repli.hdf5 %>%
    HDF5ArraySeed(name = '/rank_split') %>%
    DelayedArray

  gtf <- read.table(
    gtf_path,
    sep = '\t',
    col.names = c('chr', 'source', 'type', 'start', 'end', 'sc', 'strand', 'fr', 'annotation'),
    header = F,
    quote = ''
  ) %>%
    subset(type == 'gene') %>%
    subset(chr %in% names(chr.lengths)) %>%
    arrange(chr)
  gtf$flybase <- gtf$annotation %>% str_extract(
    'gene_id "([^"]+)"',
    group = 1
  )
  gtf <- gtf %>%
    left_join(
      gtf %>%
        split(gtf$chr) %>%
        sapply(
          \(df) df %>%
            within(
              repli_rank <- rank_repli_coverage[
                ifelse(strand == '+', start, end),
                df[1, 'chr']
              ] %>%
                as.numeric %>%
                replace(. == 0, NA)
            ) %>%
            subset(select=c(flybase, repli_rank)),
          simplify=F
        ) %>%
        do.call(rbind, .),
      by = "flybase"
    )

  num_buckets <- max(rank_repli_coverage, na.rm = T)
  gtf$quarter <- gtf$repli_rank %>%
    cut(
      seq(0, num_buckets, length.out=nquarters+1) %>%
        round
    ) %>%
    as.numeric %>%
    factor %>%
    fct_relabel(\(n) paste0("Q", n))

  metafeatures <- read.csv(metafeatures_path, row.names=1)
  features <- metafeatures %>%
    rownames_to_column %>%
    left_join(
      gtf %>% subset(select=c(flybase, repli_rank, quarter)),
      "flybase"
    ) %>%
    column_to_rownames
  features %>% subset(select=c(flybase, repli_rank, quarter))
}