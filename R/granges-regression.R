as_bulk_summarized_experiment <- function(granges, colData, sample_name = "name") {
  # For extracting row data about the samples.
  anyGranges <- granges[[1]]
  # Counts assay.
  samples <- granges %>%
    sapply(\(gr) gr$score, simplify=FALSE) %>%
    setNames(pull(colData, sample_name))
  mat <- do.call(cbind, samples)

  # ln-offset in order to get from RPKM parameters to counts.
  omat <- matrix(1, nrow=nrow(mat), ncol=ncol(mat), dimnames=dimnames(mat))
  offsets <- as.matrix(
    Diagonal(x = log(anyGranges@ranges@width) - log(1000))
    %*%
    omat
    +
    omat
    %*%
    Diagonal(x = log(colSums2(mat)) - log(1000) * 2, names=T)
  )

  rowData <- tibble(seqnames = as.factor(anyGranges@seqnames), start = anyGranges@ranges@start, width = anyGranges@ranges@width)
  SummarizedExperiment(
    list(counts = mat),
    colData = colData,
    rowData = rowData,
    metadata = list(
      # Create a GRanges without the elementMetadata (track value) which is now
      # stored in the assay.
      granges = GRanges(
        granges[[1]]@seqnames, granges[[1]]@ranges, seqlengths=seqlengths(granges[[1]])
      ),
      offset = offsets
    )
  )
}
