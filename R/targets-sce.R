targets.sce <- list(
  # Each 10X Cell Ranger -> a Seurat object to filter out doublet clusters later.
  tar_map(
    sce.data %>% mutate(obj = rlang::syms(batch)),
    names = batch,
    tar_file(
      tenx_file, tenx_path,
      cue = tar_cue("never")
    ),
    # Read 10X and immediately apply mt_pct and ribo_pct to filter cells. This
    # target is not going to filter nFeature or nCount. Doublets are low-
    # complexity (they contain both germline and somatic transcripts and are a
    # big cluster somewhere in the middle in terms of transcriptome). Keeping
    # barcodes with excessive transcripts makes the doublet clusters more
    # apparent in the next step.
    tar_target(
      seurat_qc,
      read_seurat_sctransform(
        tenx_file, batch, sce.present.features, assay.data.sc
      )
    ),
    # Plot every batch, for validation.
    tar_target(
      batch_umap,
      run_umap_on_batch(obj, metadata),
      packages = c(tar_option_get("packages"), "tidyr")
    )
  ),
  # Call clusters in the biological replicate. The purpose of this indep
  # clustering step before integration is to identify doublets.
  tar_target(nos.1, call_nos.1(seurat_qc_nos.1)),
  tar_target(nos.2, call_nos.2(seurat_qc_nos.2)),
  tar_target(tj.1, call_tj.1(seurat_qc_tj.1)),
  tar_target(tj.2, call_tj.2(seurat_qc_tj.2)),
  # Filter by nCount/nFeature, integrate, and separate cells into 5 clusters.
  # The Seurat clustering functions using the Upd_sc target's random seed and
  # the parameters in this function are what really matter for exactly
  # reproducing our clustering of the cells.
  tar_target(Upd_sc, filter_integrate_data(list(nos.1,nos.2,tj.1,tj.2))),
  
  # Write the seurat[["RNA"]]@meta.data data frame to a CSV.
  tar_file(h3.gfp.gtf, "scRNA-seq/H3-GFP-transcript-descriptive.gtf"),
  tar_file(
    assay.data.sc,
    create_assay_data_sc(
      tenx_file_nos.1, sce.features,
      flybase.annotations, flybase.gtf, h3.gfp.gtf, sce.present.features,
      'scRNA-seq-Assay-Metadata.csv')
  ),
  # Later, we will need the CDS starts in a GenomicRanges.
  tar_target(
    tss_location,
    read.csv(assay.data.sc) %>%
      with(
        setNames(
          GRanges(
            chr %>% replace(is.na(chr), "Y"),
            IRanges(
              ifelse(strand == "+", start, end) %>%
                replace(is.na(chr), -10000),
              width = 1
            )
          ),
          X
        )
      )
  ),

  # Store Seurat cell cycle feature result, as we didn't put this into the
  # Upd_sc target.
  tar_target(
    Upd_phase,
    Upd_sc %>%
      apply_cell_cycle_score(cell_cycle_drosophila, assay.data.sc) %>%
      FetchData("Phase") %>%
      rownames_to_column() %>%
      pull(Phase, rowname)
  ),
  # Cells which are (for each cluster) 33% G1, 33% S, 33% G2M -classified.
  tar_target(
    Upd_cells_by_phase,
    with(
      read.csv(metadata),
      tibble(
        rowname = X,
        ident,
        batch,
        phase = Upd_phase
      )
    ) %>%
      subsample_ident_normalize_phase() %>%
      pull(rowname)
  ),

  # Fig. S1. UMAP plots.
  tar_target(
    supplemental_elbow_figure,
    data.frame(stdev = Upd_sc[['pcasubset']]@stdev, x = 1:50) %>%
      head(10) %>%
      ggplot(aes(x, stdev^2))
      + geom_point(color = "#06470c", size = 2)
      + scale_x_continuous(breaks = c(1, 5, 10))
      + scale_y_continuous(
        trans = "sqrt", breaks = c(4, 36, 100, 150), limits = c(3, 155)
      )
      + labs(x = "PC (Germline/Somatic)", y = "Explained Variance")
      + theme_cowplot()
  ),
  tar_file(
    sc_idents_supp,
    save_figures(
      "figure/Integrated-scRNAseq", ".pdf",
      tribble(
        ~rowname, ~figure, ~width, ~height,
        "RNAseq-Validation-Subset-PCA",
        Upd_sc_plot_idents_pcasubset(Upd_sc),
        6, 4,
        "RNAseq-Validation-Subset-PCA-Batch",
        Upd_sc_plot_idents_pcasubset_batch(Upd_sc),
        5.75, 3.5,
      )
    ),
    packages = tar_option_get("packages") %>% c("cowplot", "grid", "gtable")
  ),
  tar_map(
    tibble(extension = c(".pdf", ".png")),
    tar_target(
      sc_idents,
      save_figures(
        "figure/Integrated-scRNAseq", extension,
        tribble(
          ~name, ~figure, ~width, ~height,
          "RNAseq-UMAP-Ident", Upd_sc_plot_idents(Upd_sc), 6, 4,
          "RNAseq-UMAP-Germline-Somatic",
          Upd_sc %>% Upd_sc_plot_idents %>% Upd_sc_plot_subset, 6, 3,
          "RNAseq-Quantification-Quarters-CPM",
          fpkm_third_density(log(Upd_cpm) / log(10), ylim = c(-2.75, 4.5))
          + scale_y_continuous(
            breaks = seq(-2, 4)
          )
          + theme(
            # Aspect ratio found in Fig 1 when the graphic is square (4x4 in).
            # Maintain this aspect ratio the same in +
            aspect.ratio = 1.34
          ),
          4, 4,
          "RNAseq-Quantification-Additional-CPM-Violin",
          fpkm_simple_violin(
            as.data.frame(log(Upd_cpm) / log(10))[3:5] %>%
              dplyr::rename(
                Spermatocyte="spermatocyte",
                `Other Soma`="somaticprecursor",
                Muscle="muscle"
              )
          ),
          4, 4,
          "RNAseq-Quantification-Quarters-CPM-Criteria",
          ggarrange(
            plotlist = {
              cpm_gene_lists_extended[
                c("NonExclusiveGermlineAndSomaticGenes", "ExclusiveSomaticGenes", "ExclusiveGermlineGenes", "OffGenes")
              ] %>%
                sapply(
                  \(lst) fpkm_third_density(
                    log(Upd_cpm[lst, ]) / log(10),
                    ylim = c(-2.75, 4.5),
                    inter_cutoffs = rep(list(rep(Inf, 2)), 2)
                  ) + scale_y_continuous(
                    breaks = seq(-2, 4)
                  ) + theme(
                    aspect.ratio = 1.34
                  ),
                  simplify=F
                )
            },
            nrow = 2,
            ncol = 2
          ),
          8, 8,
          "RNAseq-Quantification-Quarters-Percentage",
          gene_group_bar_plot(
            quartile.factor_Germline, quartile.factor_Somatic, Upd_cpm
          ),
          6, 5.25,
          "RNAseq-Scatter",
          ggplot(
            Upd_cpm %>%
              pmax(1e-10) %>%
              as.data.frame,
            aes(germline, somatic)
          )
          + rasterise(
            geom_point(
              stroke=NA, size=0.5,
              # size=0.75, alpha=0.25,
              color=hcl(256, 40, 25)
            ),
            dpi=240
          )
          + geom_segment(
            aes(xend=germline_end, yend=somatic_end),
            tribble(
              ~germline, ~somatic, ~germline_end, ~somatic_end,
              0, 5, Inf, 5,
              5, 0, 5, Inf
            )
          )
          + scale_x_continuous(trans="log", oob=squish, expand=rep(0.01,2), breaks=10^seq(-3,3), labels=partial(round, dig=3), limits=c(1e-3, 1e4))
          + scale_y_continuous(trans="log", oob=squish, expand=rep(0.01,2), breaks=10^seq(-3,3), labels=partial(round, dig=3), limits=c(1e-3, 1e4))
          + labs(x = "Germline CPM", y = "Somatic CPM")
          + theme_bw()
          + theme(aspect.ratio = 1),
          4.5, 4
        )
      ),
      format = "file"
    ),
    tar_target(
      sc_genes,
      save_figures(
        "figure/Integrated-scRNAseq/Genes-of-Interest", extension,
        tribble(
          ~gene,
          "AGO3",
          "RpL22-like",
          "RpL22",
          "vas",
          "nos",
          "bam",
          "tj",
          "Stat92E",
          "zfh1",
          "lncRNA:roX1",
          "lncRNA:roX2",
          "Mst87F",
          "soti",
          "sunz",
          "w-cup",
          "can",
          "Act57B",
          "Dl",
          "E(spl)m3-HLH",
          "Amy-d",
          "scpr-B",
          "wb",
          "Act5C",
          "ey",
          "eya",
          "caps",
          "hth"
        ) %>%
          rowwise %>%
          mutate(
            gene_short_name = gene %>% str_replace("lncRNA:", ""),
            name = paste0("RNAseq-FeaturePlot-", gene_short_name),
            figure = (
              Upd_sc %>% `[[<-`("DECONTX", value = Upd_decontX) %>%
                Upd_sc_feature_plot(gene, shuffle_feature_plot, assay="DECONTX")
              + labs(tag = gene_short_name)
            ) %>%
              list,
            width = 3,
            height = 2
          ) %>%
          subset(select=c(name, figure, width, height)),
        dpi = 480
      ),
      format = "file"
    ),
    tar_target(
      sc_genes_histone_mod,
      save_figures(
        "figure/Integrated-scRNAseq/Genes-Histone-Modifying", extension,
        tribble(
          ~gene,
          "Set1", "trx", "trr", "egg", "G9a", "Su(var)3-9", "ash1", "E(z)",
          "Set2", "NSD", "Su(var)3-3", "Kdm5", "Kdm4A", "Kdm4B", "Jarid2",
          "Utx", "Kdm2", "Kdm4A", "Kdm4B", "HP1b", "HP1c", "rhi", "Su(var)205",
          "HP6", "Su(var)3-7", "Su(var)2-10"
        ) %>%
          rowwise %>%
          mutate(
            gene_short_name = gene %>% str_replace("lncRNA:", ""),
            name = paste0("RNAseq-FeaturePlot-", gene_short_name),
            figure = (
              Upd_sc %>% `[[<-`("DECONTX", value = Upd_decontX) %>%
                Upd_sc_feature_plot(gene, shuffle_feature_plot, assay="DECONTX")
              + labs(tag = gene_short_name)
            ) %>%
              list,
            width = 3,
            height = 2
          ) %>%
          subset(select=c(name, figure, width, height)),
        dpi = 480
      ),
      format = "file"
    ),
    tar_target(
      sc_genes_remodeling,
      save_figures(
        "figure/Integrated-scRNAseq/Genes-Remodeling", extension,
        tribble(
          ~gene,
          "Bap60", "Bap55", "brm", "Iswi", "Acf", "Chrac-16", "Nurf-38",
          "Ino80", "Arp8", "Arp5"
        ) %>%
          rowwise %>%
          mutate(
            gene_short_name = gene %>% str_replace("lncRNA:", ""),
            name = paste0("RNAseq-FeaturePlot-", gene_short_name),
            figure = (
              Upd_sc %>% `[[<-`("DECONTX", value = Upd_decontX) %>%
                Upd_sc_feature_plot(gene, shuffle_feature_plot, assay="DECONTX")
              + labs(tag = gene_short_name)
            ) %>%
              list,
            width = 3,
            height = 2
          ) %>%
          subset(select=c(name, figure, width, height)),
        dpi = 480
      ),
      format = "file"
    ),
    tar_target(
      sc_genes_replication,
      save_figures(
        "figure/Integrated-scRNAseq/Genes-Replication", extension,
        tribble(
          ~gene,
          "PolA1", "PolA2", "Prim1", "Prim2", "PolD1", "PolD2", "PolD3",
          "Chrac-14", "PolE1", "PolE2", "PolE4", "Mcm2", "Mcm3", "dpa", "Mcm5",
          "Mcm6", "Mcm7", "Orc1", "Orc2", "Orc4", "Orc5", "Orc6"
        ) %>%
          rowwise %>%
          mutate(
            gene_short_name = gene %>% str_replace("lncRNA:", ""),
            name = paste0("RNAseq-FeaturePlot-", gene_short_name),
            figure = (
              Upd_sc %>% `[[<-`("DECONTX", value = Upd_decontX) %>%
                Upd_sc_feature_plot(gene, shuffle_feature_plot, assay="DECONTX")
              + labs(tag = gene_short_name)
            ) %>%
              list,
            width = 3,
            height = 2
          ) %>%
          subset(select=c(name, figure, width, height)),
        dpi = 480
      ),
      format = "file"
    )
  ),
  tar_file(
    Upd_volcano,
    save_figures(
      "figure/Integrated-scRNAseq", ".pdf",
      tibble(
        name="RNAseq-Volcano",
        figure=plot_volcano_apeglm(Upd_regression_somatic) %>%
          rasterise(dpi=300) %>%
          list,
        width = 3,
        height = 4
      )
    )
  ),
  tar_target(
    Upd_volcano_bivalent_genes,
    with(
      chic.gene.enrichment,
      subset(symbol, H3K4_Germline < 1e-3 & H3K27_Germline < 1e-3) %>%
        union(subset(symbol, H3K4_Somatic < 1e-3 & H3K27_Somatic < 1e-3)) %>%
        intersect(rownames(Upd_regression_somatic$map))
    )
  ),
  tar_target(
    gg_volcano_bivalent,
    Upd_regression_somatic[
      c("map", "svalue")
    ] %>%
      sapply(\(arr) arr[Upd_volcano_bivalent_genes,, drop=F], simplify=F) %>%
      plot_volcano_apeglm(
        color_column = chic.gene.enrichment %>%
          tibble(., quartile.factor_Germline, quartile.factor_Somatic) %>%
          group_by(symbol) %>%
          reframe(
            color = if (isTRUE(H3K4_Germline < 1e-3) & isTRUE(H3K27_Germline < 1e-3) & quartile.factor_Germline != "Q1") {
              if (isTRUE(H3K4_Somatic < 1e-3) & isTRUE(H3K27_Somatic < 1e-3) & quartile.factor_Somatic != "Q1")
                classification_colors_fig4$both
              else
                classification_colors_fig4$germline
            } else {
              if (isTRUE(H3K4_Somatic < 1e-3) & isTRUE(H3K27_Somatic < 1e-3) & quartile.factor_Somatic != "Q1")
                classification_colors_fig4$somatic
              else
                ""
            }
          ) %>%
          deframe() %>%
          subset(. != "")
      )
  ),
  tar_file(
    Upd_volcano_bivalent,
    save_figures(
      "figure/Integrated-scRNAseq", ".pdf",
      tibble(
        name="Bivalent-Germline-Somatic-Volcano",
        figure=gtable(unit(c(3, 1), "in"), unit(4, "in")) %>%
          gtable_add_grob(
            list(
              gg_volcano_bivalent %>% rasterise(dpi = 300) %>% ggplotGrob(),
              (
                ggplot(
                  tibble(
                    x = 1:3,
                    y = 0,
                    label = str_to_title(names(classification_colors_fig4)) %>%
                      factor(., .),
                  ),
                  aes(x, y, color=label)
                ) +
                  geom_point(size = 3) +
                  scale_color_manual(values = unlist(classification_colors_fig4, use.names=F)) +
                  labs(color = "Bivalent in:")
              ) %>%
                get_legend()
            ),
            1,
            l = 1:2
          ) %>%
          list(),
        width = 4,
        height = 4
      )
    ),
    packages = tar_option_get("packages") %>% c("cowplot", "grid", "gtable")
  ),
  tar_map(
    tibble(extension = c(".pdf", ".png")),
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
        )
      ),
      packages = c(tar_option_get("packages"), "tidyr"),
      format = "file"
    )
  ),

  tar_file(go_gene_sets_input, "Supplemental_Data/Gene_Sets_of_Interest.csv"),
  tar_target(
    go_gene_sets,
    apply(
      read.csv(go_gene_sets_input),
      2,
      \(v) v %>% subset(str_length(.) > 0),
      simplify=FALSE
    ) %>%
      enframe
  ),
  tar_file(
    go_gene_set_figure,
    save_figures(
      "figure/Integrated-scRNAseq",
      ".pdf",
      tibble(
        filename=paste0("Gene_Set_", pull(go_gene_sets, name)),
        figure=(Upd_regression_somatic$map[
          pull(go_gene_sets, value)[[1]],
          2
        ] / log(2)) %>%
          l2fc_bar_plot %>%
          list,
        width=8,
        height=c(
          Regulators_of_Nucleosome_Spacing=4,
          Origin_Binding_Components=8
        )[go_gene_sets$name]
      )
    ),
    pattern = map(go_gene_sets)
  ),

  # Unintegrated clusters - SD02. Which clusters (from sce.R RenameIdents) are
  # removed as doublets? This is justified based on expressing all of the GSC &
  # CySC markers including H3-GFP.
  tar_target(
    unintegrated_marker_genes,
    c("nos", "vas", "tj", "Egfr", "lncRNA:roX2", "Mst77F", "Act57B")
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
        `nos-Upd_H3-GFP_Rep1`=seurat_qc_nos.1 %>% FindNeighbors(dims=1:8, verb=F) %>% FindClusters(res = 0.1, verb=F),
        `nos-Upd_H3-GFP_Rep2`=seurat_qc_nos.2 %>% FindNeighbors(dims=1:9, verb=F) %>% FindClusters(res = 0.5, verb=F),
        `tj-Upd_H3-GFP_Rep1`=seurat_qc_tj.1 %>% FindNeighbors(dims=1:8, verb=F) %>% FindClusters(res = 0.15, verb=F),
        `tj-Upd_H3-GFP_Rep2`=seurat_qc_tj.2 %>% FindNeighbors(dims=1:8, verb=F) %>% FindClusters(res = 0.1, verb=F)
      ),
      unintegrated_clusters_doublets$cluster %>%
        split(factor(unintegrated_clusters_doublets$obj, c("nos.1", "nos.2", "tj.1", "tj.2"))),
      list(unintegrated_marker_genes),
      list(read.csv(metadata, row.names=1)),
      SIMPLIFY=FALSE
    ),
    packages = tar_option_get("packages") %>% c("scDblFinder", "scuttle")
  ),
  tar_file(
    fig.g1.l2fc,
    save_figures(
      "figure/Integrated-scRNAseq",
      ".pdf",
      tribble(
        ~rowname, ~figure, ~width, ~height,
        "L2FC-Scatter-Cell-Phase-Invariant-Model",
        plot_apeglm_async_phase_model(
          Upd_regression_somatic, Upd_regression_somatic_standardize_phase
        ),
        3.5,
        3.5
      )
    )
  )
)