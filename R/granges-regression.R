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

tiles_to_fseq <- function(glm, fct, fctLevels, grangesList, bw) {
  stopifnot(all(str_glue("{fct}{fctLevels}") %in% colnames(glm$Beta)))
  Mu <- split(
    exp(glm$Beta)[, str_glue("{fct}{fctLevels}")],
    grangesList %>% sapply(length) %>% enframe %>% with(Rle(factor(name, names(grangesList)), value))
  )
  for (n in names(Mu)) colnames(Mu[[n]]) <- fctLevels
  GRangesList(
    mapply(
      \(gr, Mu) list(gr) %>%
        append(
          if (nrow(Mu) > 1)
            apply(
              Mu,
              2,
              \(v) ksmooth_sliding_windows(GRanges(gr, score=v), bw=bw)$score,
              simplify=FALSE
            )
          else as.list(as.data.frame(Mu))
        ) %>%
        do.call(GRanges, .),
      grangesList,
      Mu
    )
  )
}