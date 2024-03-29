# Read bam files into a glm to generate FPKM values.
read_replicated_coverage <- function(
  bam_files,
  ref_lengths,
  feature_lengths,
  tile_width = 1000
) {
  seq_granges <- GRanges(with(enframe(ref_lengths), paste0(name, ":1-", value)))
  seq_lookup <- setNames(
    sapply(seq_along(seq_granges), \(ind) seq_granges[ind], simplify=F),
    names(ref_lengths)
  )
  ref_windows <- slidingWindows(
    seq_granges[ref_lengths == feature_lengths], tile_width, tile_width) %>%
    setNames(names(ref_lengths)[ref_lengths == feature_lengths])
  ref_windows <- sapply(
    names(ref_lengths),
    \(n) if (ref_lengths[n] == feature_lengths[n])
      ref_windows[[n]]
    else
      seq_lookup[[n]]
  ) %>%
    GRangesList
  ref_windows <- unlist(ref_windows)
  
  read_bam <- \(bam) bam %>%
    readGAlignments %>%
    countOverlaps(ref_windows, .)
  # Simplify to matrix - elements as columns, each genomic range as a row
  counts_matrix <- sapply(
    pull(bam_files, source_file, name=name),
    read_bam
  )
  # Size factor for CPM of any genomic region of interest
  size_factor = colSums(counts_matrix) / 1000 / 1000
  # Offset the rows for window size in kb
  offset = matrix(
    log(ref_windows@ranges@width) - log(1000),
    nrow = nrow(counts_matrix),
    ncol = ncol(counts_matrix)
  )
  fit = glm_gp(
    counts_matrix,
    ~ 0 + condition,
    bam_files,
    size_factors = size_factor,
    offset = offset,
    overdispersion = FALSE,
    overdispersion_shrinkage = FALSE,
    verbose = TRUE
  )
  for (c in colnames(fit$Beta)) {
    new_col <- str_replace(c, "condition", "")
    elementMetadata(ref_windows)[[new_col]] <- exp(fit$Beta[, c])
  }
}