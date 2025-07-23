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
  named_counts, named_decontx_counts, Upd_sc
) {
  # must not depend on Metadata.csv internally when building the
  # single-cell assays targets. Use single-cell Seurat object instead.
  metadata <- FetchData(Upd_sc, c("ident", "batch", "nCount_RNA", "nFeature_RNA", "pct.mito", "pct.ribo", "somatic.score"))
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
      somatic.score = mean(somatic.score) # ,
      # size_factors = sum(size_factors)
    )
}

seurat_decontx_using_batches <- function(Upd_sc, assay="RNA") {
  dcx <- decontX(GetAssayData(Upd_sc, assay=assay, layer="counts"), batch=Upd_sc$batch)
  Upd_sc[["DECONTX"]] <- CreateAssay5Object(
    counts = dcx$decontXcounts,
    key = "decontx_"
  )
  Upd_sc[["DECONTX"]]@meta.data <- Upd_sc[["RNA"]]@meta.data
  Upd_sc[["DECONTX"]]@misc$decontX <- dcx[-match("decontXcounts", names(dcx))]
  Upd_sc[["DECONTX"]] %>% NormalizeData
}

seurat_assay_isoforms <- function(Upd_sc, assay.data.sc, Upd_cpm_transcript_to_use, isoforms_named_list) {
  counts <- GetAssayData(Upd_sc, "RNA", "counts")
  meta.data <- FetchData(Upd_sc, c("ident", "batch"))
  assay.data.sc <- assay.data.sc %>% read.csv(row.names = 1)
  for (n in names(isoforms_named_list)) {
    colnames(isoforms_named_list[[n]]) <- colnames(isoforms_named_list[[n]]) %>% paste(n, ., sep="_")
  }

  batch_data <- sapply(
    isoforms_named_list,
    \(data) {
      data <- as.matrix(data)
      r1 <- assay.data.sc[rownames(counts), "flybase"] %>% replace(!(. %in% rownames(Upd_cpm_transcript_to_use)), NA)
      c1 <- as.character(meta.data[colnames(data), "ident"])
      read_isoforms <- Upd_cpm_transcript_to_use[
        cbind(
          rep(
            assay.data.sc[rownames(counts), "flybase"] %>% replace(!(. %in% rownames(Upd_cpm_transcript_to_use)), NA),
            ncol(data)
          ),
          rep(as.character(meta.data[colnames(data), "ident"]), each=nrow(counts))
        )
      ]
      repl_numeric <- data[
        cbind(
          read_isoforms,
          rep(colnames(data), each=nrow(counts))
        )
      ]
      matrix(
        repl_numeric,
        nrow = nrow(counts),
        dimnames = list(rownames(counts), colnames(data))
      ) %>%
        as("sparseMatrix")
    }
  )
  repl_rows <- names(which(rowAnys(as.matrix(is.na(batch_data[[1]])))))
  for (n in names(isoforms_named_list)) {
    batch_data[[n]][repl_rows, ] <- counts[repl_rows, colnames(batch_data[[n]])]
  }
  exon_counts <- do.call(cbind, setNames(batch_data, NULL))
  stopifnot(isTRUE(all.equal(dimnames(counts), dimnames(exon_counts))))

  assay <- CreateAssay5Object(counts = exon_counts, key = "exon_")
  assay@meta.data <- Upd_sc[["RNA"]]@meta.data
  assay
}
