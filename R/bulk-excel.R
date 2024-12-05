publish_chic_fragments <- function(
    frag_size_tables_publish,
    mapq_tables_publish,
    peaks_publish,
    output_path) {
  wb <- createWorkbook()

  title <- "H3 ChIC"
  addWorksheet(wb, title)
  start_row <- 1
  writeData(wb, title, "Fragment Sizes", colNames = F, startRow = start_row, startCol = 1)
  start_row <- start_row + 1
  for (n in names(frag_size_tables_publish)) {
    writeData(wb, title, n, colNames = F, startRow = start_row, startCol = 1)
    frag_size_tables_publish[[n]]$fragment_sizes <- frag_size_tables_publish[[n]]$fragment_sizes %>%
      factor(unique(.))
    tbl <- dcast(frag_size_tables_publish[[n]], chr ~ fragment_sizes, fun = sum, value.var = "n")
    writeDataTable(
      wb,
      title,
      tbl,
      startRow = start_row + 1, startCol = 1
    )
    for (iloc in seq(nrow(tbl))) {
      tbl[iloc, -1] <- tbl[iloc, -1] %>%
        `/`(sum(.)) %>%
        round(3)
    }
    writeDataTable(
      wb,
      title,
      tbl,
      startRow = start_row + 1, startCol = 1 + ncol(tbl) + 1
    )
    addStyle(
      wb, title,
      createStyle(numFmt = "0.0%"),
      rows = start_row + 1 + seq(nrow(tbl)),
      cols = seq(1 + ncol(tbl) + 1 + 1, 1 + ncol(tbl) + ncol(tbl)),
      gridExpand = TRUE
    )
    start_row <- start_row + 1 + 1 + nrow(tbl) + 1
  }

  writeData(wb, title, "MAPQ Values", colNames = F, startRow = start_row, startCol = 1)
  start_row <- start_row + 1
  for (n in names(mapq_tables_publish)) {
    writeData(wb, title, n, colNames = F, startRow = start_row, startCol = 1)
    tbl <- dcast(mapq_tables_publish[[n]], chr ~ mapq_range, fun = sum, value.var = "n")
    writeDataTable(
      wb,
      title,
      tbl,
      startRow = start_row + 1, startCol = 1
    )
    for (iloc in seq(nrow(tbl))) {
      tbl[iloc, -1] <- tbl[iloc, -1] %>%
        `/`(sum(.)) %>%
        round(3)
    }
    writeDataTable(
      wb,
      title,
      tbl,
      startRow = start_row + 1, startCol = 1 + ncol(tbl) + 1
    )
    addStyle(
      wb, title,
      createStyle(numFmt = "0.0%"),
      rows = start_row + 1 + seq(nrow(tbl)),
      cols = seq(1 + ncol(tbl) + 1 + 1, 1 + ncol(tbl) + ncol(tbl)),
      gridExpand = TRUE
    )
    start_row <- start_row + 1 + 1 + nrow(tbl) + 1
  }

  publish_chic_fragments_additional_data(
    wb,
    "ChIC Peak Calling",
    peaks_publish
  )

  publish_peaks_contingency_tables(
    wb,
    "ChIC Mono&bivalency",
    peaks_publish
  )

  dir.create(dirname(output_path), showW = F, rec = F)
  saveWorkbook(wb, output_path, overwrite = T)
  output_path <- output_path
}

publish_chic_fragments_additional_data <- function(
    wb,
    title,
    df,
    percent_columns = 13:18) {
  addWorksheet(wb, title)
  writeData(
    wb,
    title,
    matrix(
      "One-Tailed p-value"
    ),
    startCol = 7, startRow = 1,
    colNames = F
  )
  writeData(
    wb,
    title,
    matrix(
      "Gene Body % Enriched (% of sliding windows p < 0.001)"
    ),
    startCol = 13, startRow = 1,
    colNames = F
  )
  writeData(
    wb,
    title,
    matrix(
      "CDS Head Enriched (bp from TSS along strand p < 0.001)"
    ),
    startCol = 19, startRow = 1,
    colNames = F
  )
  writeDataTable(
    wb,
    title,
    df,
    startCol = 1, startRow = 2,
    withFilter = T,
    tableStyle = excel_tables$deepgreen
  )
  # Write % number format.
  addStyle(
    wb,
    title,
    createStyle(numFmt = "0%"),
    rows = seq(3, 2 + nrow(df)),
    cols = percent_columns,
    gridExpand = TRUE
  )
}

publish_peaks_contingency_tables <- function(
    wb,
    title,
    df) {
  addWorksheet(wb, title)
  writeData(
    wb,
    title,
    matrix(
      c(
        "Germline Valency (p < 0.001)",
        "no H3K4",
        "enr H3K4"
      )
    ),
    startCol = 1, startRow = 1,
    colNames = F
  )
  tbl <- tribble(
    ~`no H3K27`, ~`enr H3K27`,
    sum(
      df$H3K4_Germline >= 0.001 &
        df$H3K27_Germline >= 0.001,
      na.rm = T
    ),
    sum(
      df$H3K4_Germline >= 0.001 &
        df$H3K27_Germline < 0.001,
      na.rm = T
    ),
    sum(
      df$H3K4_Germline < 0.001 &
        df$H3K27_Germline >= 0.001,
      na.rm = T
    ),
    sum(
      df$H3K4_Germline < 0.001 &
        df$H3K27_Germline < 0.001,
      na.rm = T
    )
  )
  writeData(
    wb,
    title,
    tbl,
    startCol = 2, startRow = 2,
    colNames = F
  )

  writeData(
    wb,
    title,
    matrix(
      c(
        "Somatic Valency (p < 0.001)",
        "no H3K4",
        "enr H3K4"
      )
    ),
    startCol = 5, startRow = 1,
    colNames = F
  )
  tbl <- tribble(
    ~`no H3K27`, ~`enr H3K27`,
    sum(
      df$H3K4_Somatic >= 0.001 &
        df$H3K27_Somatic >= 0.001,
      na.rm = T
    ),
    sum(
      df$H3K4_Somatic >= 0.001 &
        df$H3K27_Somatic < 0.001,
      na.rm = T
    ),
    sum(
      df$H3K4_Somatic < 0.001 &
        df$H3K27_Somatic >= 0.001,
      na.rm = T
    ),
    sum(
      df$H3K4_Somatic < 0.001 &
        df$H3K27_Somatic < 0.001,
      na.rm = T
    )
  )
  writeData(
    wb,
    title,
    tbl,
    startCol = 6, startRow = 2,
    colNames = F
  )
}

detail_repli_analysis <- function(
    assay.data.sc,
    flybase.sequence.ontology,
    repli_Germline, repli_Somatic, repli_Bayes_Factor,
    CPM_Germline, CPM_Somatic,
    Upd_regression_somatic,
    chic.gene.enrichment,
    chromosome_pericetromere_label, transcriptome_integrate = 25000
) {
  assay.data.sc <- assay.data.sc %>%
    read.csv() %>%
    dplyr::rename(symbol = "X")
  flybase.sequence.ontology <- flybase.sequence.ontology %>%
    read.table(col.names = c("gene_primary_id", "gene_symbol", "so_term_name", "so_term_id")) %>%
    as_tibble()
  genes <- with(
    assay.data.sc,
    GRanges(chr %>% replace_na("*"), IRanges(ifelse(strand == "+", start, end) %>% replace(is.na(chr), 1), width = 1), names = symbol)
  )
  gene_lookup <- genes %>%
    findOverlaps(repli_Germline) %>%
    sapply(\(v) if (length(v)) v[1] else NA)
  feature <- flybase.sequence.ontology %>%
    mutate(
      so_term_name = so_term_name %>%
        factor(
          c(
            "pseudogene", "retrogene", "protein_coding_gene",
            "antisense_lncRNA_gene",
            "lncRNA_gene", "miRNA_gene", "ncRNA_gene",
            "snoRNA_gene", "snRNA_gene", "tRNA_gene",
            "C_D_box_scaRNA_gene", "C_D_box_snoRNA_gene",
            "gene_array_member"
          )
        ) %>%
        dplyr::recode(
          C_D_box_scaRNA_gene="CDbox_gene",
          C_D_box_snoRNA_gene="CDbox_gene"
        )
    ) %>%
    group_by(gene_primary_id) %>%
    summarise(
      so_term_name = factor(
        levels(so_term_name)[min(c(Inf, as.numeric(so_term_name)), na.rm=T)],
        levels(so_term_name)
      )
    )
  df <- tibble(
    assay.data.sc[c(1, 3, 6)],
    region = chr %>%
      mapply(
        \(n, matches) paste0(
          n,
          if (length(matches)) "C"
        ),
        .,
        findOverlaps(genes, chromosome_pericetromere_label) %>%
          as("List")
      ) %>%
      factor(
        c(
          "2L", "2LC", "2RC", "2R", "3L", "3LC", "3RC", "3R", "4",
          "X", "Y",
          "rDNA"
        )
      ),
    TSS = ifelse(assay.data.sc$strand == "+", assay.data.sc$start, assay.data.sc$end),
    feature = deframe(feature)[flybase],
    timing = cbind(
      GSC = repli_Germline$score[gene_lookup],
      CySC = repli_Somatic$score[gene_lookup]
    ),
    bayes_factor = repli_Bayes_Factor$score[gene_lookup],
  )
  gene_neighbor <- findOverlaps(
    GenomicRanges::resize(genes, transcriptome_integrate, "center"),
    genes
  ) %>%
    as("List") %>%
    as.list() %>%
    replace(
      which(!(assay.data.sc$chr %in% names(chr.lengths))),
      as.list(which(!(assay.data.sc$chr %in% names(chr.lengths))))
    )
  quant_raw <- cbind(
    GSC = (log(CPM_Germline) / log(10)) %>% replace(!is.finite(.), NA),
    CySC = (log(CPM_Somatic) / log(10)) %>% replace(!is.finite(.), NA)
  )
  quant_region <- sapply(
    gene_neighbor,
    \(v) colMeans(quant_raw[v,, drop=F], na.rm=T)
  ) %>%
    t()
  df <- tibble(
    df,
    quant_raw,
    quant_region,
    L2FC = -Upd_regression_somatic$map[
      match(symbol, rownames(Upd_regression_somatic$map)), 2
    ] / log(2),
    logp = -log(as.matrix(chic.gene.enrichment[7:12])) / log(10),
    logp_region = sapply(
      gene_neighbor,
      \(v) colMeans(logp[v,, drop=F], na.rm=T)
    ) %>%
      t()
  )
}

publish_repli_analysis <- function(
    assay.data.sc, repli_Germline, repli_Somatic, repli_Bayes_Factor,
    window_name, CPM_Germline, CPM_Somatic,
    chromosome_pericetromere_label,
    output_path) {
  wb <- createWorkbook()

  title <- "Repliseq"
  addWorksheet(wb, title)

  assay.data.sc <- assay.data.sc %>%
    read.csv() %>%
    dplyr::rename(symbol = "X")
  seqlevels(repli_Germline) <- seqlevels(repli_Germline) %>% c("*")
  TSS <- with(
    assay.data.sc,
    GRanges(
      chr %>% replace_na("*"),
      IRanges(
        ifelse(strand == "+", start, end) %>% replace(is.na(chr), 1), width = 1
      ),
      names = symbol
    )
  )
  gene_lookup <- findOverlaps(TSS, repli_Germline) %>%
    sapply(\(v) if (length(v)) v[1] else NA)
  df <- tibble(
    assay.data.sc,
    region = rep("", nrow(assay.data.sc)) %>%
      replace(from(findOverlaps(TSS, chromosome_pericetromere_label)), "Pericentromere"),
    GSC = repli_Germline$score[gene_lookup],
    CySC = repli_Somatic$score[gene_lookup],
    `Bayes Factor` = repli_Bayes_Factor$score[gene_lookup],
    `GSC CPM` = CPM_Germline,
    `CySC CPM` = CPM_Somatic,
    OnOff_CPM = structure(
      as.numeric(interaction(CPM_Germline >= 5, CPM_Somatic >= 5)),
      levels = c("GSC Off CySC Off", "GSC On CySC Off", "GSC Off CySC On", "GSC On CySC On"),
      class = "factor"
    ),
    `Diff Rep Assay` = window_name,
  )
  writeData(
    wb,
    title,
    matrix(
      c(
        "-1 Late to +1 Early",
        NA,
        "P(Alt. Hypot.) / P(Null Hypot.)",
        "Transcriptome Values",
        NA,
        NA,
        "Nested Peak Calling (** signif) TSS in named region"
      ),
      nrow = 1
    ),
    startCol = 12, startRow = 1,
    colNames = F
  )
  writeDataTable(
    wb,
    title,
    df,
    startCol = 1, startRow = 2,
    withFilter = T,
    tableStyle = excel_tables$deepgreen
  )
  # Shrink number of decimals.
  addStyle(
    wb,
    title,
    createStyle(numFmt = "0.00"),
    rows = seq(3, 2 + nrow(df)),
    cols = 10:14,
    gridExpand = TRUE
  )

  dir.create(dirname(output_path), showW = F, rec = F)
  saveWorkbook(wb, output_path, overwrite = T)
  output_path <- output_path
}
