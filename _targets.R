# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# HDF5Array is used for glm and will write to custom location.
# HDF5Array::setHDF5DumpDir('D:/HDF5Array_dump')

# Load packages required to define the pipeline:
library(grDevices)
library(magrittr)
library(reshape2)
library(stringr)
library(targets)
library(tarchetypes)
library(tibble)
library(viridis)

pdfFonts(sans = pdfFonts()$Helvetica)
postscriptFonts(sans = postscriptFonts()$Helvetica)

# Set target options:
tar_option_set(
  packages = c(
    'apeglm',
    'BiocIO',
    'Cairo',
    'cowplot',
    'DESeq2',
    'dplyr',
    'forcats',
    'GenomicRanges',
    'ggnewscale',
    'ggpattern',
    'ggplot2',
    'ggpubr',
    'ggrastr',
    'glmGamPoi',
    'HDF5Array',
    'irlba',
    'magrittr',
    'Matrix',
    'OpenImageR',
    'openxlsx',
    'orthogene',
    'pillar',
    'processx',
    'purrr',
    'readxl',
    'reshape2',
    'reticulate',
    'rtracklayer',
    'S4Vectors',
    'scales',
    'scran',
    'scuttle',
    'Seurat',
    'SingleCellExperiment',
    'stringr',
    'tibble',
    'viridis',
    'withr',
    'zoo'
  )
)

# tar_make_clustermq() is an older (pre-{crew}) way to do distributed computing
# in {targets}, and its configuration for your machine is below.
options(clustermq.scheduler = "multiprocess")

Sys.setenv(BOWTIE_THREADS="12")

future::plan(future::multicore, workers=8)

# All plots will have white background.
ggplot2::theme_set(
  ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "#000000"),
      panel.background = ggplot2::element_rect(fill="transparent"),
      # We want the facet text to appear as a title of the graphic. Not as a
      # detail akin to the axis break labels!
      strip.text = ggplot2::element_text(size = ggplot2::rel(1.2))
    )
)

# tar_make_future() is an older (pre-{crew}) way to do distributed computing
# in {targets}, and its configuration for your machine is below.
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

sce.data = tibble(
  batch=c('nos.1','nos.2','tj.1','tj.2'),
  capitalized=c('Nos','Nos','tj','tj'),
  batch_num=c('1','2','1','2'),
  dir_name=str_glue("scRNA-seq/{capitalized}-Upd_H3-GFP_Rep{batch_num}"),
  tenx_path=str_replace(str_glue("{dir_name}/outs/filtered_feature_bc_matrix"), "tj", "Tj"),
  bam_path=str_glue("{dir_name}/outs/possorted_genome_bam.bam")
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

sce.clusters = data.frame(cluster = c('germline','somatic','spermatocyte','muscle'))
sce.clusters <- tribble(
  ~ cluster, ~ contrast,
  # Contrasts taken from Upd_celltype_model_matrix and hardcoded into the
  # targets here.
  'germline', c(1, CySCoverGSC=-0.5, 0,0,0),
  'somatic', c(1, CySCoverGSC=0.5, 0,0,0),
  'spermatocyte', c(1,-0.5,1,0,0),
  'somaticprecursor', c(1,0.5,0,1,0),
  'muscle', c(1,0.5,0,0,1)
)

list(
  tar_file(
    flybase.annotations.current,
    'fbgn_annotation_ID_fb_2022_04.tsv.gz'
  ),
  tar_file(
    flybase.annotations.previous,
    'fbgn_annotation_ID_fb_2020_06.tsv.gz'
  ),
  tar_target(
    flybase.annotations,
    c(
      `2022_04`=flybase.annotations.current,
      `2020_06`=flybase.annotations.previous
    )
  ),
  tar_target(
    sce.features,
    load_feature_names(tenx_file_nos.1, flybase.annotations)
  ),
  tar_target(
    sce.present.features,
    select_features(tenx_file_nos.1, flybase.annotations)
  ),
  tar_target(
    flybase.lengths,
    bulk_reads_idxstats_chic.bam_GC3772016_S1_L002_chr %>%
      subset(rname != "*") %>%
      pull(rlength, rname)
  ),
  tar_target(
    feature.lengths,
    flybase.lengths %>% replace(!(names(.) %in% names(chr.lengths)), 1)
  ),
  tar_target(
    feature.rle,
    Rle(factor(names(feature.lengths), names(feature.lengths)), as.numeric(feature.lengths))
  ),
 
  tar_target(Upd_model_matrix, build_model_matrix(FetchData(Upd_sc, c("ident", "batch")), Upd_decontX_contamination)),
  tar_target(
    shuffle_feature_plot,
    sample(Cells(Upd_sc))
  ),
  tar_target(
    Upd_sc_size_factors,
    Upd_sc %>%
      `[[<-`("EXONICRNA", value = Upd_exons) %>%
      pooled_size_factors_seurat_ident_autosome(assay = "EXONICRNA")
  ),
  tar_target(
    metadata,
    extract_upd_metadata_to_csv(Upd_sc, Upd_sc_size_factors, "scRNA-seq-Metadata.csv"),
    format='file'
  ),
  # GLM fitting: Prior step required to determine the nuisance parameter and
  # apply the apeglm model for SupplementalTable2.
  tar_target(
    Upd_glm,
    fit_glm(Upd_exons, Upd_model_matrix, metadata),
  ),
  # Fit glmGamPoi to the 33/33/33 G1/S/G2M data for the supplement.
  tar_target(
    Upd_glm_standardize_phase,
    fit_glm(
      CreateSeuratObject(Upd_exons) %>%
        `[`(, Upd_cells_by_phase) %>%
        GetAssay("RNA"),
      Upd_model_matrix[Upd_cells_by_phase, ],
      metadata
    )
  ),
  # SupplementalTable2 data.
  tar_target(
    Upd_regression_somatic,
    # prior_var was determined by fitting Upd_glm with ident + batch and without
    # the decontX terms (which grow the coef of interest when there is a trend
    # in the data where contamination is negatively correlated with the fitted
    # values based on ident). Then apeglm:::priorVar was applied. We want to
    # shrink LFC estimates based on a prior from a highly interpretable model
    # (ident + batch), not weaken the regularization because we purposefully
    # added additional complexity to the model.
    apeglm_coef_table_sample(Upd_glm, shrinkage_cutoff = 0, prior_var = 0.8935249),
    cue = tar_cue("never")
  ),
  # Fit apeglm to the 33/33/33 G1/S/G2M data for the supplement.
  tar_target(
    Upd_regression_somatic_standardize_phase,
    with_options(
      list(future.globals.maxSize = 5 * 1024^3),
      apeglm_coef_timebound(Upd_glm_standardize_phase, prior_var = 0.8027215)
    ),
    packages = tar_option_get("packages") %>% c("future.apply", "R.utils")
  ),
  # For cell cycle scoring.
  tar_download(
    cell_cycle_drosophila,
    "https://github.com/hbc/tinyatlas/raw/add6f25/cell_cycle/Drosophila_melanogaster.csv",
    "cell_cycle_drosophila.csv"
  ),
  tar_map(
    tibble(extension = c(".pdf", ".png")),
    tar_target(
      supplemental_figures,
      save_figures(
        "figure/Integrated-scRNAseq", extension,
        tribble(
          ~name, ~figure, ~width, ~height,
          "RNAseq-Validation-Elbow", supplemental_elbow_figure, 3, 2.5,
          "RNAseq-UMAP-Genotype",
          Upd_sc %>%
            AddMetaData(recode(Upd_sc$batch, nos.1="nos", nos.2="nos", tj.1="tj", tj.2="tj"), "genotype") %>%
            Upd_sc_group_plot("genotype", shuffle_feature_plot)
          + scale_color_manual(
            values = c(
              nos="#afd34d",
              tj="#b85dd4"
            )
          ),
          3, 2,
        )
      ),
      format = "file"
    )
  ),
  # Predict gene mean from glmGamPoi. This is mostly still used for H3-GFP:
  # The 10X BAM scripts will focus on the reference chromosomes and will later
  # test normalizing by transcript length using the original reference. That
  # analysis is not actually needed for the H3-GFP construct.
  tar_target(
    Upd_cpm_regression,
    glm_make_cpm_table(Upd_glm)
  ),
  tar_target(
    Upd_fpkm_regression_do_not_use_for_quantification,
    Upd_cpm_regression %>% cpm_to_fpkm_using_cds(assay.data.sc)
  ),
  tar_target(
    Upd_cpm_transcript_to_use,
    cpm_select_transcript_to_quantify(
      Upd_count_transcripts,
      assay.data.sc
    ),
    packages = tar_option_get("packages") %>% c("tidyr")
  ),
  tar_target(
    Upd_isoform_exonic_length,
    lookup_mat_transcripts_exon_length(
      Upd_count_transcripts, Upd_cpm_transcript_to_use, assay.data.sc
    )
  ),
  # TPM values based on exon length correction, then scale to sum to 1MM.
  tar_target(
    Upd_tpm_do_not_use_for_quantification,
    join_cpm_data(
      assay.data.sc,
      Upd_count_transcripts %>%
        subset(select=c(germline, somatic, spermatocyte, somaticprecursor, muscle)) %>%
        table_to_tpm %>%
        as.data.frame %>%
        cbind(exon_length = Upd_count_transcripts$exon_length, .),
      Upd_cpm_transcript_to_use,
      Upd_fpkm_regression_do_not_use_for_quantification %>%
        table_to_tpm()
    ) %>%
      table_to_tpm()
  ),
  # CPM values based on dominant isoform (by max FPKM or TPM value for the
  # isoform). We will fix the values to sum to 1MM because we are adding back in
  # the H3-GFP value from Upd_cpm_regression, but we do not apply the correction
  # by exonic length before calculating this TPM.
  tar_target(
    Upd_cpm,
    join_cpm_data(
      assay.data.sc,
      Upd_count_transcripts %>%
        subset(select=c(germline, somatic, spermatocyte, somaticprecursor, muscle)) %>%
        table_to_tpm,
      Upd_cpm_transcript_to_use,
      Upd_cpm_regression %>% table_to_tpm,
      corr=F) %>%
      table_to_tpm
  ),
  tar_target(
    # Single Cell Excel supports SupplementalTable1 Quantification_Results.
    sc_excel,
    publish_excel_results(
      Upd_regression_somatic,
      Upd_cpm,
      assay.data.sc,
      'scRNA-seq-Regression/SupplementalTable1_2_outputs.xlsx'
    ),
    format='file'
  ),
  tar_map(
    chic.fpkm.data,
    names = name,
    tar_target(
      quartile.factor,
      quant_quartile_factor(Upd_cpm[, tolower(name)], q1_threshold=5),
      packages = tar_option_get("packages") %>% c("tidyr")
    )
  ),
  tar_target(
    expression.factor,
    interaction(
      quartile.factor_Germline != "Q1",
      quartile.factor_Somatic != "Q1"
    ) %>%
      `levels<-`(value = c("off", "GSC", "CySC", "both")) %>%
      fct_relevel("off", after=4) %>%
      setNames(names(quartile.factor_Germline))
  ),

  # Paths to bowtie2-calling scripts. Used by the chic and repli targets.
  tar_file(align_chic_lightfiltering, "scripts/align_chic_lightfiltering.sh"),
  tar_file(align_repli_lightfiltering, "scripts/align_repli_lightfiltering.sh"),
  tar_file(align_repli, "scripts/align_repli.sh"),

  # Description: Heatmaps of ChIC vs sorted genes. ####
  tar_map(
    tibble(
      chic.experiments,
      chic.heatmap.tss = rlang::syms(
        str_glue("chic.heatmap.tss_{mark}_{name}_CN_chr")
      )
    ),
    names = experiment_name,

    # Molecule/input FPKM ratio track, in bigwig format.
    tar_target(
      tss_mark_heatmap,
      display_tss_tile_matrix(
        chic.heatmap.tss,
        paste0("figure/", name, "/FPKM-", mark, ".pdf"),
        scale_image_filter = rep(1/50, 50),
        fc_max = 4,
        fc_filter = Inf
      ),
      format = "file"
    ),
    tar_target(
      tss_mark_heatmap_png,
      display_tss_tile_matrix(
        chic.heatmap.tss,
        paste0("figure/", name, "/FPKM-", mark, ".png"),
        scale_image_filter = rep(1/50, 50),
        fc_max = 4,
        fc_filter = Inf
      ),
      format = "file"
    ),
    tar_target(
      inv_tss_mark_heatmap,
      display_tss_tile_matrix(
        chic.heatmap.tss,
        paste0("figure/", name, "/INV-FPKM-", mark, ".pdf"),
        scale_image_filter = rep(1/50, 50),
        fc_max = 4,
        fc_filter = Inf,
        direction = -1
      ),
      format = "file"
    ),
    tar_target(
      inv_tss_mark_heatmap_png,
      display_tss_tile_matrix(
        chic.heatmap.tss,
        paste0("figure/", name, "/INV-FPKM-", mark, ".png"),
        scale_image_filter = rep(1/50, 50),
        fc_max = 4,
        fc_filter = Inf,
        direction = -1
      ),
      format = "file"
    )
  ),

  tar_map(
    data.frame(extension=c(".png", ".pdf")),
    tar_target(
      tss_mark_heatmap_legend,
      {
        dir.create("figure/Germline", showW=FALSE, recursive=TRUE)
        dir.create("figure/Somatic", showW=FALSE, recursive=TRUE)
        save_figures(
          "figure",
          extension,
          tribble(
            ~name, ~figure, ~width, ~height,
            "Germline/Heatmap-Legend",
            get_legend(
              ggplot(
                # color from fc_max
                data.frame(x=0, y=0, c=c(0, 4)),
                aes(x, y, color=c)
              ) + geom_point() + scale_color_viridis_c(
                option="magma",
                guide=guide_colorbar(title = "mark/input", barheight = 10)
              )
            ),
            1,
            3,
            "Germline/INV-Heatmap-Legend",
            get_legend(
              ggplot(
                # color from fc_max
                data.frame(x=0, y=0, c=c(0, 4)),
                aes(x, y, fill=c)
              ) + geom_tile() + create_direction_invert_tss_tile_matrix_gradient() +
                guides(
                  fill=guide_colorbar(title = "mark/input", barheight = 10)
                )
            ),
            1,
            3,
            "Somatic/Heatmap-Legend",
            get_legend(
              ggplot(
                # color from fc_max
                data.frame(x=0, y=0, c=c(0, 4)),
                aes(x, y, color=c)
              ) + geom_point() + scale_color_viridis_c(
                option="magma",
                guide=guide_colorbar(title = "mark/input", barheight = 10)
              )
            ),
            1,
            3,
            "Somatic/INV-Heatmap-Legend",
            get_legend(
              ggplot(
                # color from fc_max
                data.frame(x=0, y=0, c=c(0, 4)),
                aes(x, y, fill=c)
              ) + geom_tile() + create_direction_invert_tss_tile_matrix_gradient() +
                guides(
                  fill=guide_colorbar(title = "mark/input", barheight = 10)
                )
            ),
            1,
            3
          )
        )
      },
      format = "file"
    )
  ),

  tar_map(
    mutate(
      chic.fpkm.data,
      sc_chr_quartile_data_TSS = rlang::syms(str_glue("sc_chr_quartile_data_{name}_TSS")),
      sc_chr_active_data_TSS = rlang::syms(str_glue("sc_chr_active_data_{name}_TSS"))
    ),
    names = name,
    tar_file(
      fig.fpkm.chic.facet,
      save_figures(
        paste0("figure/", name),
        ".pdf",
        tribble(
          ~rowname, ~figure, ~width, ~height,
          "CHIC-TSS-Chr-AllMarks-RNAseq-Quartile", 
          chic_plot_average_profiles_facet_grid(
            dplyr::rename(sc_chr_quartile_data_TSS, genes=quant),
            "CPM Quartile",
            setNames(sc_quartile_colors, NULL)
          ) +
            coord_cartesian(NULL, chic_average_profile_limits, ex=F),
          9,
          12,
          "CHIC-TSS-Chr-AllMarks-RNAseq",
          chic_plot_average_profiles_facet_grid(
            dplyr::rename(sc_chr_active_data_TSS, genes=activity),
            "CPM Quartile",
            rep(chic_line_track_colors[[tolower(name)]], 2)
          ) +
            coord_cartesian(NULL, chic_average_profile_limits, ex=F) +
            guides(linewidth = guide_none(), color = guide_none()),
          8,
          10
        )
      )
    )
  ),
  tar_file(
    fig.repli.chic.diff.timing,
    save_figures(
      "figure/Both-Cell-Types",
      ".pdf",
      tribble(
        ~rowname, ~figure, ~width, ~height,
        "CHIC-TSS-Diff-Replication-Program",
        list(
          Germline=repli_diff_timing_data_Germline_TSS,
          Somatic=repli_diff_timing_data_Somatic_TSS
        ) %>%
          bind_rows(.id = "celltype") %>%
          facet_diff_replication_program(),
        8,
        8
      )
    ),
    packages = tar_option_get("packages") %>% c("grid", "gtable")
  ),

  tar_target(
    germline_somatic_line_plot_legend,
    (
      ggplot(tibble(name=c("Germline", "Somatic"), x=0, y=0, xend=1, yend=0), aes(x, y))
      +
      geom_segment(aes(xend=xend, yend=yend, color=name), linewidth=0.66)
      +
      scale_color_manual(values = chic_line_track_colors %>% unlist %>% setNames(NULL))
    ) %>%
      get_legend
  ),
  tar_file(
    fig.fpkm.chic.both,
    save_figures(
      "figure/Both-Cell-Types",
      ".pdf",
      tribble(
        ~rowname, ~figure, ~width, ~height,
        "CHIC-TSS-Chr-AllMarks-RNAseq",
        list(Germline=sc_chr_active_data_Germline_TSS, Somatic=sc_chr_active_data_Somatic_TSS) %>%
          bind_rows(.id = "name") %>%
          arrange(desc(row_number())) %>%
          mutate(genes=interaction(activity, ordered(name))) %>%
          chic_plot_average_profiles_facet_grid(
            "CPM Quartile",
            c(muted(chic_line_track_colors$germline, l=70), chic_line_track_colors$germline, muted(chic_line_track_colors$somatic, l=70), chic_line_track_colors$somatic),
            linewidth = c(0.33, 0.66, 0.33, 0.66)
          ) %>%
          `+`(coord_cartesian(NULL, chic_average_profile_limits, ex=F)) %>%
          replace_legend(germline_somatic_line_plot_legend),
        12,
        8
      )
    )
  ),

  # Plots of ChIC profiles. ####
  tar_map(
    tibble(
      chic.fpkm.data,
      figure_dir = str_glue("figure/{name}"),
      plot.chic.peak.location = rlang::syms(str_glue("plot.chic.peak.location_{name}")),
      sc_quartile_data_TSS = rlang::syms(str_glue("sc_quartile_data_{name}_TSS")),
      sc_quartile_data_Paneled = rlang::syms(str_glue("sc_quartile_data_{name}_Paneled")),
      sc_nucleosome_quartile_data_TSS = rlang::syms(str_glue("sc_nucleosome_quartile_data_{name}_TSS")),
      sc_nucleosome_exclusive_quartile_data_TSS = rlang::syms(str_glue("sc_nucleosome_exclusive_quartile_data_{name}_TSS")),
    ),
    names = name,
    tar_file(
      chic.results,
      save_figures(
        figure_dir,
        ".pdf",
        tribble(
          ~name, ~figure, ~width, ~height,
          "CHIC-TSS-AllMarks-RNAseq-Quartile",
          chic_plot_average_profiles_facet_grid(
            dplyr::rename(sc_quartile_data_TSS, genes="quant"),
            "Quant.",
            setNames(sc_quartile_annotations, NULL),
            seq(0.5, 0.85, length.out=4),
            facet_wrap(vars(mark), scales = "free")
          ) +
            geom_blank(
              aes(x=0, y=value, color=NULL, linewidth=NULL, group=NULL),
              chic_lineplot_limit_data_TSS
            ),
          9, 4,
          "CHIC-AllMarks-RNAseq-Quartile",
          chic_plot_paneled_profiles_facet_grid(
            dplyr::rename(sc_quartile_data_Paneled, genes="quant"),
            "Quant.",
            setNames(sc_quartile_annotations, NULL),
            seq(0.5, 0.85, length.out=4),
            facet_wrap(vars(mark))
          ),
          15, 4,
          "CHIC-TSS-H3-RNAseq-Quartile",
          chic_plot_h3_enrichment(
            tibble(sc_nucleosome_quartile_data_TSS, genes=quant, mark="H3"),
            legend_title = "Quant.",
            quartile_colors = setNames(sc_quartile_annotations, NULL),
            linewidth = seq(0.5, 0.85, length.out=4),
            faceter = facet_wrap(vars(mark))
          ),
          4, 4,
          paste0("CHIC-TSS-H3-RNAseq-", name, "Exclusive-Quartile"),
          chic_plot_h3_enrichment(
            tibble(sc_nucleosome_exclusive_quartile_data_TSS, genes=quant, mark="H3"),
            legend_title = "Quant.",
            quartile_colors = setNames(sc_quartile_annotations, NULL),
            linewidth = seq(0.5, 0.85, length.out=4),
            faceter = facet_wrap(vars(mark))
          ),
          4, 4,
        )
      )
    )
  ),

  targets.bulk.samples,
  targets.chic,
  targets.flybase,
  targets.quantification,
  targets.repli,
  targets.rnaseq,
  targets.sce
)
 