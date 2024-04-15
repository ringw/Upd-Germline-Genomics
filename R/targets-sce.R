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
    sc_genes_dendrogram(Upd_sc, gene_ = rowAnys(Upd_cpm > 100)[rownames(Upd_sc[["RNA"]])])
  ),
  tar_target(
    Upd_genes_dendrogram,
    sc_sort_genes_dendrogram(Upd_cpm, Upd_genes_dendrogram_raw),
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
  )
)