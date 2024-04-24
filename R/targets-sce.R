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
  tar_target(
    Upd_genes_dendrogram_raw,
    sc_genes_dendrogram(Upd_sc, gene_ = rowAnys(Upd_tpm > 100)[rownames(Upd_sc[["RNA"]])])
  ),
  tar_target(
    Upd_genes_dendrogram,
    sc_sort_genes_dendrogram(Upd_tpm, Upd_genes_dendrogram_raw),
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
    analyze_pcasubset_batch_effect(Upd_sc)
  ),
  tar_target(
    manova_ident,
    analyze_pcasubset_ident(Upd_sc, cell_cycle_drosophila, assay.data.sc)
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
  )
)