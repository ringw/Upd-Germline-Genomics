# Scan the 10X Genomics BAM file - to verify the faster perl-based workflow.
read_tenx_barcoded_bam_exons <- function(bam) {
  scanBam(
    BamFile(bam),
    param=ScanBamParam(
      what=c("rname", "pos"), tag=c("CB", "TX"), tagFilter=list(xf=25L, RE="E")
    )
  )
}

# Produce list of FBtr transcripts, for looking up Matrix Market row index in Perl.
read_transcript_ids <- function(assay_data, gtf_path) {
  assay_data <- assay_data %>% read.csv(row.names = 1)
  # Use exons so that we can count length.
  gtf <- read.table(
    gtf_path,
    sep = '\t',
    col.names = c('chr', 'source', 'type', 'start', 'end', 'sc', 'strand', 'fr', 'annotation'),
    header = F,
    quote = ''
  ) %>% subset(grepl('exon|miRNA', type))
  gtf$flybase <- gtf$annotation %>% str_extract(
    'gene_id "([^"]+)"',
    group = 1
  )
  gtf$tx <- gtf$annotation %>% str_extract(
    'transcript_id "([^"]+)"',
    group = 1
  )
  gtf_transcripts <- gtf %>%
    group_by(flybase, tx) %>%
    summarize(exon_length = sum(abs(end - start) + 1), .groups = "drop")
  left_join(
    data.frame(flybase = assay_data$flybase %>% setdiff(c("", NA))),
    gtf_transcripts,
    "flybase"
  )
}

# Calls BAM to feature matrix coords Perl script for exonic fragments -> UMI
# count matrix.
cell_ranger_bam_feature_exonic_umi <- function(
  bam_path, cell.barcodes.txt, cell_ranger_bam_feature_matrix_coords,
  feature_ids_path, feature_tag_name = c("GX", "TX"), output_path
) {
  feature_tag_name <- match.arg(feature_tag_name)
  run(
    "sh",
    c(
      "-c",
      str_glue(
        "rclone cat 'sharepoint:Chen Lab Backups/Upd_sc/'",
        str_replace(bam_path, "scRNA-seq/", ""),
        " | ",
        "samtools view -@ 4 -D CB:{cell.barcodes.txt} | ",
        "perl {cell_ranger_bam_feature_matrix_coords} ",
        "{cell.barcodes.txt} {feature_ids_path} {feature_tag_name}",
        " > {output_path}"
      )
    )
  )
  output_path
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

cbind_pseudobulk_decontaminated <- function(
  named_counts, named_decontx_counts, metadata
) {
  metadata <- metadata %>% read.csv(row.names = 1)
  metadata$ident <- metadata$ident %>% factor(sce.clusters$cluster)
  # Use nos.2 as the reference level for batch. It has the greatest number of
  # cells in the rare clusters.
  metadata$batch <- metadata$batch %>%
    factor(c("nos.2", "nos.1", "tj.1", "tj.2"))
  addCellIdCols <- function(name, mat) {
    colnames(mat) <- paste0(name, "_", colnames(mat))
    mat
  }
  which_rows <- (
    !duplicated(rownames(named_counts[[1]]))
    & !is.na(rownames(named_counts[[1]]))
  )
  SingleCellExperiment(
    list(
      counts = do.call(
        cbind,
        mapply(addCellIdCols, names(named_counts), named_counts, SIMPLIFY=F)
      )[which_rows, rownames(metadata)],
      decontXcounts = do.call(
        cbind,
        named_decontx_counts
      )[which_rows, rownames(metadata)]
    ),
    colData = metadata
  ) %>%
    pseudobulk(
      vars(ident, batch),
      aggregation_functions = list(
        .default = "rowSums2"
      ),
      nCount_RNA = sum(nCount_RNA),
      nFeature_RNA = sum(nFeature_RNA),
      pct.mito = mean(pct.mito),
      pct.ribo = mean(pct.ribo),
      somatic.score = mean(somatic.score),
      size_factors = sum(size_factors)
    )
}

create_decontx_assay_fbgn_objects <- function(decontx_objects, Upd_sc_metadata, assay_data) {
  Upd_sc_metadata <- Upd_sc_metadata %>% read.csv(row.names = 1)
  assay_data <- assay_data %>% read.csv(row.names = 1)
  # Just need to subset the RNA assay features (from assay_data) because we
  # didn't check for the "GX" tag for genes such as H3-GFP.
  genes_df <- rownames_to_column(assay_data) %>%
    subset(
      flybase %in% rownames(decontx_objects[[1]]$decontXcounts),
      select=c(rowname, flybase)
    )
  counts <- do.call(
    cbind,
    sapply(
      decontx_objects,
      \(l) l$decontXcounts[genes_df$flybase, ],
      simplify=F
    )
  )[, rownames(Upd_sc_metadata)]
  rownames(counts) <- genes_df$rowname
  assay_obj <- CreateAssayObject(
    counts = counts,
    min.cells = 5,
    key = "decontx_"
  )
  assay_obj@misc$decontX <- sapply(
    decontx_objects,
    \(obj) obj[-match("decontXcounts", names(obj))],
    simplify=F
  )
  assay_obj %>% NormalizeData
}

seurat_decontx_using_batches <- function(Upd_sc) {
  dcx <- decontX(GetAssayData(Upd_sc, assay="RNA", layer="counts"), batch=Upd_sc$batch)
  Upd_sc[["DECONTX"]] <- CreateAssayObject(
    counts = dcx$decontXcounts,
    key = "decontx_"
  )
  Upd_sc[["DECONTX"]]@misc$decontX <- dcx[-match("decontXcounts", names(dcx))]
  Upd_sc[["DECONTX"]] %>% NormalizeData
}