targets.sce <- list(
  tar_target(
    Upd_cells_dendrogram_raw,
    sc_cells_dendrogram(Upd_sc)
  ),
  tar_target(
    Upd_cells_dendrogram,
    sc_sort_cells_dendrogram(Upd_sc, Upd_cells_dendrogram_raw),
    packages = tar_option_get("packages") %>% c("dendextend")
  ),
  # Report genes with "medium" or "high" quantification in any cluster.
  tar_target(
    Upd_genes_dendrogram_loose_criteria,
    sc_genes_dendrogram(Upd_sc, gene_ = rowAnys(Upd_cpm > 10)[rownames(Upd_sc[["RNA"]])]) %>%
      sc_sort_genes_dendrogram(Upd_cpm, .),
    packages = tar_option_get("packages") %>% c("dendextend")
  ),
  tar_file(
    Upd_gene_clustering_report,
    Upd_genes_dendrogram_loose_criteria %>%
      cut(h = 165) %>%
      with(lower) %>%
      sapply(labels) %>%
      setNames(str_glue("dendrogram_{seq_along(.)}")) %>%
      report_cpm_in_gene_sets(Upd_cpm, ., "scRNA-seq-Regression/Report-Gene-CPM.html", "scRNA-seq-Regression/Report-Gene-CPM.pdf"),
    packages = tar_option_get("packages") %>% c("gt")
  ),
  tar_target(
    Upd_genes_dendrogram_raw,
    sc_genes_dendrogram(Upd_sc, gene_ = rowAnys(Upd_cpm > 50)[rownames(Upd_sc[["RNA"]])])
  ),
  tar_target(
    Upd_genes_dendrogram,
    sc_sort_genes_dendrogram(Upd_cpm, Upd_genes_dendrogram_raw),
    packages = tar_option_get("packages") %>% c("dendextend")
  ),
  tar_target(
    Upd_genes_dendrogram_cpm_raw,
    Upd_cpm %>%
      subset(rowAlls((. > 0) %>% replace_na(FALSE))) %>%
      log %>%
      `%*%`(
        diag(
          table(
            read.csv(metadata, row.names=1)$ident
          )[colnames(.)]
        )
      ) %>%
      dist %>%
      hclust(method = "average") %>%
      as.dendrogram,
    packages = tar_option_get("packages") %>% c("tidyr")
  ),
  tar_target(
    Upd_genes_dendrogram_cpm,
    sc_sort_genes_dendrogram(Upd_cpm, Upd_genes_dendrogram_cpm_raw),
    packages = tar_option_get("packages") %>% c("dendextend")
  ),
  tar_file(
    Upd_genes_heatmap,
    save_figures(
      "figure/Integrated-scRNAseq",
      ".pdf",
      tibble(
        "Heatmap-With-Clusters",
        figure = plot_single_cell_heatmap(
          Upd_sc, Upd_genes_dendrogram, Upd_cells_dendrogram
        ) %>%
          list,
        width=8,
        height=6
      )
    ),
    packages = tar_option_get("packages") %>% c("dendextend", "ggdendro")
  ),
  tar_target(
    Upd_genes_newcut, cut(Upd_genes_dendrogram, 158)
  ),
  tar_target(
    Upd_genes_newcut_minor_breaks,
    c(1, cumsum(sapply(Upd_genes_newcut$lower, \(dd) length(labels(dd)))))
  ),
  tar_target(
    Upd_genes_newcut_breaks,
    (
      Upd_genes_newcut_minor_breaks[-length(Upd_genes_newcut_minor_breaks)]
      + Upd_genes_newcut_minor_breaks[-1]
    ) / 2
  ),
  tar_file(
    Upd_genes_logCPM_heatmap,
    save_figures(
      "figure/Integrated-scRNAseq",
      ".pdf",
      tibble(
        c("Heatmap-logCPM", "Heatmap-logCPM-Annotated"),
        figure = list(
          plot_logcpm_heatmap(
            Upd_cpm, Upd_genes_dendrogram
          )
          + theme(aspect.ratio = 4/3),
          plot_logcpm_heatmap(
            Upd_cpm, Upd_genes_dendrogram, genes_cut=158
          )
          + theme(aspect.ratio = 4/3)
          + scale_y_continuous(
            pos = "right",
            breaks = sort(c(Upd_genes_newcut_minor_breaks, Upd_genes_newcut_breaks)),
            labels = as.character(rbind("", head(LETTERS, length(Upd_genes_newcut_breaks))) %>% c(""))
          )
        ),
        width=4.5,
        height=6
      )
    ),
    packages = tar_option_get("packages") %>%
      c("dendextend", "ggdendro", "tidyr")
  ),
  tar_file(
    Upd_genes_heatmap_excel,
    publish_heatmap_named_cuts(
      Upd_genes_newcut$lower %>% setNames(head(LETTERS, length(.))),
      assay.data.sc,
      "scRNA-seq-Regression/Dendrogram-Gene-Names.xlsx"
    )
  ),
  tar_map(
    tibble(extension = c(".pdf", ".png")),
    tar_target(
      sc_track_figure,
      save_figures(
        "figure/Integrated-scRNAseq", extension,
        tribble(
          ~name, ~figure, ~width, ~height,
          "RNAseq-log2FC-Track",
          plot_genes_on_chrs(Upd_regression_somatic$map, 2, assay.data.sc, do_smooth=FALSE)
          + labs(x = NULL, y = "log2FC(CySC / GSC)"),
          8,
          2.5,
          "RNAseq-log2FC-Track-LOESS",
          plot_genes_on_chrs(Upd_regression_somatic$map, 2, assay.data.sc, do_smooth=TRUE)
          + labs(x = NULL, y = "log2FC(CySC / GSC)"),
          8,
          2.5
        )
      ),
      format = "file"
    ),
    tar_target(
      sc_plot_genes_on_chrs,
      save_figures(
        "figure/Integrated-scRNAseq", extension,
        tribble(
          ~name, ~figure, ~width, ~height,
          "RNAseq-AvgFC-X-vs-AA",
          plot_chr_ratio_on_clusters(Upd_sc),
          6,
          8,
          "RNAseq-AvgFC-X-vs-AA-Demo-Batch-Variance",
          plot_chr_ratio_on_clusters(Upd_sc) + facet_wrap(vars(batch)),
          12,
          16
        )
      ),
      packages = c(tar_option_get("packages"), "tidyr"),
      format = "file"
    )
  ),
  tar_target(
    manova_batch_effect,
    analyze_pcasubset_batch_effect(Upd_sc),
    cue = tar_cue("never")
  ),
  tar_target(
    manova_ident,
    analyze_pcasubset_ident(Upd_sc, cell_cycle_drosophila, assay.data.sc),
    cue = tar_cue("never")
  ),
  tar_target(
    manova_gt,
    gt_group(
      manova_batch_effect %>%
        analyze_manova(c("genotype", "batch_effect"), 3) %>%
        dcast(response ~ term, value.var="R2") %>%
        gt("response", cap="Batch R-squared") %>%
        fmt_number(decimals = 2),
      manova_ident %>%
        analyze_manova(c("ident", "Phase"), 3) %>%
        dcast(response ~ term, value.var="R2") %>%
        gt("response", cap="Cell identity R-squared") %>%
        fmt_number(decimals = 2)
    ),
    packages = c(tar_option_get("packages"), "gt")
  ),

  # Unintegrated clusters - SD02. Which clusters (from sce.R RenameIdents) are
  # removed as doublets? This is justified based on expressing all of the GSC &
  # CySC markers including H3-GFP.
  tar_target(
    unintegrated_marker_genes,
    c("nos", "vas", "tj", "Egfr", "lncRNA:roX2", "Mst77F", "wb", "Act57B")
  ),
  tar_target(
    unintegrated_clusters_doublets,
    tribble(
      ~ obj, ~ cluster, ~ label,
      "nos.1", "3", "doublet",
      "nos.2", "8", "doublet",
      "nos.2", "9", "doublet",
      "tj.1", "2", "doublet",
      "tj.2", "1", "doublet"
    )
  ),
  tar_target(
    unintegrated_clusters_report,
    mapply(
      unintegrated_report_cluster_expression,
      list(
        seurat_qc_nos.1 %>% FindNeighbors(dims=1:8, verb=F) %>% FindClusters(res = 0.1, verb=F),
        seurat_qc_nos.2 %>% FindNeighbors(dims=1:9, verb=F) %>% FindClusters(res = 0.5, verb=F),
        seurat_qc_tj.1 %>% FindNeighbors(dims=1:8, verb=F) %>% FindClusters(res = 0.15, verb=F),
        seurat_qc_tj.2 %>% FindNeighbors(dims=1:8, verb=F) %>% FindClusters(res = 0.1, verb=F)
      ),
      unintegrated_clusters_doublets$cluster %>%
        split(factor(unintegrated_clusters_doublets$obj, c("nos.1", "nos.2", "tj.1", "tj.2"))),
      list(unintegrated_marker_genes),
      SIMPLIFY=FALSE
    )
  )
)