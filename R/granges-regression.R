as_bulk_summarized_experiment <- function(granges, colData, sample_name = "name") {
  # For extracting row data about the samples.
  unl <- unlist(granges[[1]])
  # Counts assay.
  samples <- granges %>%
    sapply(\(gr) unlist(gr)$score, simplify=FALSE) %>%
    setNames(pull(colData, sample_name))
  mat <- do.call(cbind, samples)

  # ln-offset in order to get from RPKM parameters to counts.
  omat <- matrix(1, nrow=nrow(mat), ncol=ncol(mat), dimnames=dimnames(mat))
  offsets <- as.matrix(
    Diagonal(x = log(unl@ranges@width) - log(1000))
    %*%
    omat
    +
    omat
    %*%
    Diagonal(x = log(colSums2(mat)) - log(1000) * 2, names=T)
  )

  rowData <- tibble(seqnames = as.factor(unl@seqnames), start = unl@ranges@start, width = unl@ranges@width)
  SummarizedExperiment(
    list(counts = mat),
    colData = colData,
    rowData = rowData,
    metadata = list(
      grangesList = GRangesList(
        # GRanges constructor creates a new GRanges from the first sample
        # without including the track metadata on the windows (which is now
        # stored in assay).
        sapply(granges[[1]], \(gr) GRanges(gr@seqnames, gr@ranges, seqlengths=seqlengths(gr)))
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