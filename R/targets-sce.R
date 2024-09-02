targets.sce <- list(
  # Fig. 1. UMAP plots.
  tar_file(
    sc_idents_supp,
    save_figures(
      "figure/Integrated-scRNAseq", ".pdf",
      tibble(
        rowname="RNAseq-Validation-Subset-PCA",
        figure=Upd_sc_plot_idents_pcasubset(Upd_sc) %>% list,
        width=6, height=4
      )
    )
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
          "RNAseq-Quantification-Quarters-SCT",
          fpkm_quarter_density(
            with(
              list(
                sct = cbind(
                  germline = sctransform_quantile$germline[, "90%"],
                  somatic = sctransform_quantile$somatic[, "90%"]
                )
              ),
              sct %>% subset(rowAlls(abs(.) >= 0.01))
            ),
            y_label = bquote(Q["90%"]*"(SCT)"),
            ylim = c(-0.5, 5)
          ),
          4, 4,
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

  # Dendrogram for supp. gene expression.
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
        width=5,
        height=4
      )
    ),
    packages = tar_option_get("packages") %>%
      c("dendextend", "ggdendro", "tidyr")
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
        width = 6,
        height = 4
      )
    )
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
  tar_target(
    sd02_xlsx,
    publish_sd02(
      unintegrated_clusters_report,
      "scRNA-seq-Regression/SD02-scRNA-seq-Cell-Level.xlsx"
    ),
    packages = tar_option_get("packages") %>% union("openxlsx")
  )
)