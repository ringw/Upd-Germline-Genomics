Read10XBamExons <- function(bam, transcript_df, cells_to_use) {
  coords_tmp <- tempfile(fileext = ".coord.txt")
  mm_tmp <- tempfile(fileext = ".mtx")
  mydimnames <- list(transcript_df$tx, cells_to_use)
  granges <- bam %>%
    read_tenx_barcoded_bam_exons
  mat <- granges %>%
    match_umis_bam(transcript_df, cells_to_use, coords_tmp) %>%
    coords_to_matrix_market(mydimnames, length(granges[[1]]$tag$CB), mm_tmp) %>%
    readMM
  dimnames(mat) <- mydimnames
  mat
}

read_tenx_barcoded_bam_exons <- function(bam) {
  scanBam(
    BamFile(bam),
    param=ScanBamParam(
      what=c("rname", "pos"), tag=c("CB", "TX"), tagFilter=list(xf=25L, RE="E")
    )
  )
}

read_transcript_ids <- function(assay_data, gtf_path) {
  assay_data <- assay_data %>% read.csv(row.names = 1)
  gtf <- read.table(
    gtf_path,
    sep = '\t',
    col.names = c('chr', 'source', 'type', 'start', 'end', 'sc', 'strand', 'fr', 'annotation'),
    header = F,
    quote = ''
  ) %>% subset(grepl('RNA', type)) # include "transcript" types in the reference
  gtf$flybase <- gtf$annotation %>% str_extract(
    'gene_id "([^"]+)"',
    group = 1
  )
  gtf$tx <- gtf$annotation %>% str_extract(
    'transcript_id "([^"]+)"',
    group = 1
  )
  left_join(
    data.frame(flybase = assay_data$flybase %>% setdiff(c("", NA))),
    gtf %>% subset(select=c(flybase, tx)),
    "flybase"
  )
}

match_umis_bam_do_not_use <- function(bam, transcript_df, cells) {
  stopifnot(typeof(bam) == "list" && length(bam) == 1)
  
  all_transcripts_str <- bam[[1]]$tag$TX
  cell_column <- match(bam[[1]]$tag$CB, cells)
  all_transcripts_str <- all_transcripts_str[!is.na(cell_column)]
  cell_column <- cell_column[!is.na(cell_column)]
  reframe(
    tibble(transcripts = all_transcripts_str, column = cell_column) %>%
      rowwise,
    row = strsplit(
      strsplit(transcripts, ";")[[1]],
      ","
    ) %>%
      sapply(\(v) v[1]),
    column
  )
}

match_umis_bam <- function(bam, transcript_df, cells, output_coords) {
  stopifnot(typeof(bam) == "list" && length(bam) == 1)
  
  all_transcripts_str <- bam[[1]]$tag$TX
  cell_column <- match(bam[[1]]$tag$CB, cells)
  all_transcripts_str <- all_transcripts_str[!is.na(cell_column)]
  cell_column <- cell_column[!is.na(cell_column)]

  write.table(
    data.frame(comment = character(0)),
    output_coords,
    row.names=F,
    col.names=F,
    quote=F
  )
  for (i in seq_along(cell_column)) {
    data <- data.frame(
      row = strsplit(
        strsplit(all_transcripts_str[i], ";")[[1]],
        ","
      ) %>%
        sapply(\(v) v[1]) %>%
        match(transcript_df$tx) %>%
        subset(!is.na(.))
    )
    if (length(data$row))
      data$column <- cell_column[i]
    write.table(
      data,
      output_coords,
      row.names=F,
      col.names=F,
      quote=F,
      append=T
    )
  }
  return(output_coords)
}

coords_to_matrix_market <- function(coords_file, dimnames, n_counts, mm_output) {
  write.table(
    data.frame(comment = "%%MatrixMarket matrix coordinate integer general"),
    mm_output,
    row.names=F,
    col.names=F,
    quote=F
  )
  with_options(
    list(scipen = 100),
    write.table(
      data.frame(
        nrow = length(dimnames[[1]]),
        ncol = length(dimnames[[2]]),
        nent = n_counts
      ),
      mm_output,
      row.names=F,
      col.names=F,
      sep="\t",
      append=T
    )
  )
  run(
    "bash",
    c(
      "-c",
      str_glue("sort {coords_file} | uniq -c | awk '{{print $2 \"\\t\" $3 \"\\t\" $1}}' >> {mm_output}")
    )
  )
  return(mm_output)
}

assemble_decont_to_glm_offset <- function(
  counts_mats, deconts_objects, genes_to_use, metadata
) {
  for (i in seq_along(counts_mats)) {
    # Correct our issue with colnames(counts) not having the add'l prefix.
    colnames(counts_mats[[i]]) <- colnames(deconts_objects[[i]]$decontXcounts)
  }

  metadata <- metadata %>% read.csv(row.names = 1)
  make_offset <- function(counts, deconts) {
    # If the pattern of zeroes does not match, then it's unclear how we will
    # compute the offset which has log link function.
    stopifnot(isTRUE(all.equal(counts != 0, deconts$decontXcounts != 0)))
    # Subset the genes to use within this function.
    counts <- counts[genes_to_use, ]
    deconts$decontXcounts <- deconts$decontXcounts[genes_to_use, ]
    # Given the decontaminated value which can be estimated by our regression
    # model, then the count that is actually observed comes from "counts". Apply
    # this correction within the offset, which has log (natural) link function.
    offset <- counts
    offset@x <- log(offset@x) - log(deconts$decontXcounts@x)
    return(offset)
  }
  offset_mats <- mapply(
    make_offset, counts_mats, deconts_objects, SIMPLIFY=FALSE
  )
  offset <- do.call(cbind, offset_mats)[, rownames(metadata)]
  offset + rep(log(metadata$size_factors), each=nrow(offset))
}

combine_counts_mats <- function(named_counts_mats) {
  for (n in names(named_counts_mats))
    colnames(named_counts_mats[[n]]) <- colnames(named_counts_mats[[n]]) %>%
      paste0(n, "_", .)
  do.call(cbind, named_counts_mats)
}