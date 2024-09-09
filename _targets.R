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
    'normr',
    'OpenImageR',
    'openxlsx',
    'orthogene',
    'pillar',
    'processx',
    'propagate',
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

repli.processing = data.frame(
  name=c('.q20.bedgraph', '.q20.markdup.bedgraph', '.q20.p10.markdup.bedgraph', '.q20.near10.bedgraph', '.q20.near25.bedgraph'),
  display=c('OrigPileup', 'DupePosition', 'DownsampleAndDupePosition', 'DeleteNearbyPositions', 'DeleteNearbyPositions25')
)
repli.processing = data.frame(
  name='.q20.near10.bedgraph',
  display='DeleteNearbyPositions'
)

repli.coverage.contrasts = list(
  EarlyLate = c(1,0,0,-1),
  Weighted = c(7/8, 5/8, 3/8, 1/8)
)

repli.coverage <- list(
  ident = c('Tj','Tj','GSC','GSC'),
  contrast = c('EarlyLate', 'Weighted','EarlyLate','Weighted')
)

chic.raw.tracks <- tar_map(
  tibble(
    chic.samples,
    chr_target = rlang::syms(paste0("chic.bam_", sample, "_chr"))
  ),
  names = sample,
  unlist = FALSE,
  tar_file(
    chic.bam,
    chr_target
  ),
  tar_target(
    chic.raw,
    bam_paired_fragment_ends(chic.bam, feature.lengths),
    cue = tar_cue("never")
  )
)

# ChIC lookup: For every driver x mark x input/mod aggregation.
chic.lookup <- chic.fpkm.data %>%
  cross_join(chic.mark.data) %>%
  cross_join(tribble(~input, "input", "mod"))
# Pull the sample names from the sample sheet csv. For each target chic.bam and
# chic.raw, create rlang syms referring to the list of targets.
chic.lookup$files <- mapply(
  \(mark, driver, input) chic.samples[
    chic.samples$driver == driver
      & chic.samples$group == mark
      & if (input == "input") chic.samples$molecule == "H3" else chic.samples$molecule != "H3",
    "sample"
  ] %>%
    paste0("chic.bam_", .) %>%
    rlang::syms() %>%
    append(list("c"), .) %>% do.call(call, ., quote=T),
  chic.lookup$mark,
  chic.lookup$driver,
  chic.lookup$input,
  SIMPLIFY = FALSE
)
chic.lookup$sample_names <- mapply(
  \(mark, driver, input) chic.samples[
    chic.samples$driver == driver
      & chic.samples$group == mark
      & if (input == "input") chic.samples$molecule == "H3" else chic.samples$molecule != "H3",
    "sample"
  ] %>%
    paste0("chic.raw_", .),
  chic.lookup$mark,
  chic.lookup$driver,
  chic.lookup$input,
  SIMPLIFY = FALSE
)
chic.lookup$samples <- sapply(
  chic.lookup$sample_names,
  rlang::syms,
  simplify = FALSE
)
chic.lookup <- chic.lookup %>% within(
  lookup_name <- paste(input, mark, name, sep="_")
)
chic.experiments$chic_input_files <- chic.experiments %>%
  mutate(input = "input") %>%
  left_join(
    chic.lookup %>% mutate(chic_input_files = files),
    by = c("name", "driver", "mark", "input")
  ) %>%
  pull(chic_input_files)
chic.experiments$chic_mod_files <- chic.experiments %>%
  mutate(input = "mod") %>%
  left_join(
    chic.lookup %>% mutate(chic_mod_files = files),
    by = c("name", "driver", "mark", "input")
  ) %>%
  pull(chic_mod_files)

sce_targets <- tar_map(
  unlist = FALSE,
  sce.data %>% mutate(obj = rlang::syms(batch)),
  names = batch,
  tar_file(
    tenx_file, tenx_path,
    cue = tar_cue("never")
  ),
  # Our light QC (before SCTransform) will be based on mitochondrial percent and
  # ribo transcript percent. This is not the full QC. Later we will filter
  # nCount_RNA and nFeature_RNA. At this stage, these criteria did not need to
  # be applied to get clear clusters.
  tar_target(
    seurat_qc,
    read_seurat_sctransform(
      tenx_file, batch, sce.present.features, assay.data.sc
    )
  ),
  # Quantiles of gene expression value per cluster. We will apply one of the
  # quantiles as a nonparametric statistic for calling some genes "off", because
  # their % expressed is too low, or because the normalized expression value at
  # the low end is tiny.
  tar_target(
    sctransform_quantile,
    read_seurat_sctransform(
      tenx_file, batch, sce.features, assay.data.sc,
      return.only.var.genes = FALSE, run_pca = FALSE
    ) %>%
      quantify_quantiles(metadata),
    cue = tar_cue("never")
  ),
  tar_target(
    lognormalize_quantile,
    read_seurat_sctransform(
      tenx_file, batch, sce.present.features, assay.data.sc,
      return.only.var.genes = FALSE, run_pca = FALSE
    ) %>%
      quantify_quantiles(
        metadata, assay = "RNA", slot = "data", per_million = FALSE
      )
  ),
  # Plot every batch, for validation.
  tar_target(
    batch_umap,
    run_umap_on_batch(obj, metadata),
    packages = c(tar_option_get("packages"), "tidyr")
  )
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
 
  tar_file(h3.gfp.gtf, "scRNA-seq/H3-GFP-transcript-descriptive.gtf"),
  tar_file(
    assay.data.sc,
    create_assay_data_sc(
      tenx_file_nos.1, sce.features,
      flybase.annotations, flybase.gtf, h3.gfp.gtf, sce.present.features,
      'scRNA-seq-Assay-Metadata.csv')
  ),
  tar_target(nos.1, call_nos.1(seurat_qc_nos.1)),
  tar_target(nos.2, call_nos.2(seurat_qc_nos.2)),
  tar_target(tj.1, call_tj.1(seurat_qc_tj.1)),
  tar_target(tj.2, call_tj.2(seurat_qc_tj.2)),
  tar_target(Upd_sc, filter_integrate_data(list(nos.1,nos.2,tj.1,tj.2))),
  tar_target(Upd_model_matrix, build_model_matrix(FetchData(Upd_sc, c("ident", "batch")), Upd_decontX_contamination)),
  tar_target(
    Upd_model_matrix_add_phase_covariate,
    build_model_matrix_add_phase_covariate(
      Idents(Upd_sc),
      Upd_sc %>%
        apply_cell_cycle_score(cell_cycle_drosophila, assay.data.sc) %>%
        FetchData("Phase") %>%
        pull(Phase)
    )
  ),
  tar_target(
    Upd_phase,
    Upd_sc %>%
      apply_cell_cycle_score(cell_cycle_drosophila, assay.data.sc) %>%
      FetchData("Phase") %>%
      rownames_to_column() %>%
      pull(Phase, rowname)
  ),
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
  tar_combine(
    sctransform_quantile,
    sce_targets$sctransform_quantile,
    command = combine_gene_quantiles(list(!!!.x))
  ),
  tar_combine(
    lognormalize_quantile,
    sce_targets$lognormalize_quantile,
    command = combine_gene_quantiles(list(!!!.x))
  ),
  tar_target(
    Upd_glm,
    fit_glm(Upd_exons, Upd_model_matrix, metadata),
  ),
  tar_target(
    Upd_glm_g1,
    fit_glm(
      CreateSeuratObject(Upd_exons) %>%
        `[`(, names(Upd_phase)[Upd_phase == "G1"]) %>%
        GetAssay("RNA"),
      Upd_model_matrix[names(Upd_phase)[Upd_phase == "G1"], ],
      metadata
    )
  ),
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
  tar_target(
    Upd_regression_somatic_g1,
    with_options(
      list(future.globals.maxSize = 5 * 1024^3),
      apeglm_coef_timebound(Upd_glm_g1, prior_var = 0.8055104)
    ),
    packages = tar_option_get("packages") %>% c("future.apply", "R.utils")
  ),
  sce_targets,
  # Pseudobulk by the genotype (Nos-GAL4 or tj-GAL4).
  tar_target(
    supplemental_genes,
    "tj,vas" %>% strsplit(",") %>% unlist
  ),
  tar_target(
    supplemental_bulk_cpm,
    gene_pseudobulk_cpm(
      paste0(
        "scRNA-seq-Regression/pseudobulk_",
        c("nos", "nos", "tj", "tj"),
        c(".1", ".2", ".1", ".2"),
        "_filtered.txt"
      ),
      mapply(
        # Load all of the 10X matrices using sce.features (all features instead
        # of our min # cells features).
        \(n, f) read_seurat_sctransform(
          f, n, sce.features, assay.data.sc, sctransform=FALSE
        ),
        c("nos.1", "nos.2", "tj.1", "tj.2"),
        c(tenx_file_nos.1, tenx_file_nos.2, tenx_file_tj.1, tenx_file_tj.2),
        SIMPLIFY = FALSE
      ),
      c("nos.1", "nos.2", "tj.1", "tj.2"),
      flybase.gtf,
      metadata,
      assay.data.sc
    )
  ),
  tar_target(
    supplemental_gene_list,
    c(
      "Act5C", "alphaTub84B", "AGO3", "vas", "tj", "zfh1",
      "soti", "w-cup",
      "lncRNA:roX2",
      "wb",
      "Act57B"
    )
  ),
  tar_map(
    sce.clusters,
    names = cluster,
    tar_target(
      supplemental_cluster_cpm,
      gene_cluster_cpm(
        Upd_cpm,
        list(nos.1=nos.1, nos.2=nos.2, tj.1=tj.1, tj.2=tj.2),
        cluster,
        metadata
      )
    )
  ),
  tar_target(
    supplemental_cluster_dot_plot_figure,
    dot_plot_cpm(
      list(
        GSC=supplemental_cluster_cpm_germline,
        CySC=supplemental_cluster_cpm_somatic,
        `other(germ)`=supplemental_cluster_cpm_spermatocyte,
        `other(soma)`=supplemental_cluster_cpm_somaticprecursor,
        mus=supplemental_cluster_cpm_muscle
      ),
      supplemental_gene_list,
      logcpm_max=3.55
    )
  ),
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
  # For cell cycle scoring.
  tar_download(
    cell_cycle_drosophila,
    "https://github.com/hbc/tinyatlas/raw/add6f25/cell_cycle/Drosophila_melanogaster.csv",
    "cell_cycle_drosophila.csv"
  ),
  tar_target(
    cell_cycle_scoring_excel,
    tibble(
      filename = "scRNA-seq-Regression/Cell-Cycle-Scoring.xlsx",
      do_write = write_Upd_sc_cell_cycle_phases(Upd_sc, cell_cycle_drosophila, assay.data.sc) %>%
        write_excel_tables_list_percentages(filename)
    ) %>%
      pull(filename),
    format = "file"
  ),
  tar_target(
    cpm_gene_lists,
    apply(
      Upd_cpm[, c("germline", "somatic")],
      2,
      \(v) (v >= 5) %>% which %>% names,
      simplify=FALSE
    )
  ),
  tar_target(
    cpm_gene_venn,
    plot_gene_lists(
      cpm_gene_lists
    )
  ),
  tar_map(
    tibble(extension = c(".pdf", ".png")),
    tar_target(
      supplemental_figures,
      save_figures(
        "figure/Integrated-scRNAseq", extension,
        tribble(
          ~name, ~figure, ~width, ~height,
          "RNAseq-Integrated-Cluster-Dots", supplemental_cluster_dot_plot_figure, 3, 4,
          "RNAseq-Validation-Elbow", supplemental_elbow_figure, 3, 2.5,
          "RNAseq-CPM-Gene-List-Venn-Area-Blank", cpm_gene_venn, 3.6, 2.7,
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
          "RNAseq-Batch-UMAP-Doublets",
          plot_multiple_umap_data(
            list(nos.1=batch_umap_nos.1, nos.2=batch_umap_nos.2, tj.1=batch_umap_tj.1, tj.2=batch_umap_tj.2) %>%
              bind_rows(.id = "batch")
          ) + theme(
            aspect.ratio = 1,
            axis.text = element_text(size = 8),
            panel.spacing = unit(1, "lines"),
            panel.background = element_rect(fill = "white")
          ),
          6,
          6
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
    Upd_fpkm_regression,
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
      Upd_fpkm_regression %>% table_to_tpm) %>%
      table_to_tpm
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
    sc_excel,
    publish_excel_results(
      Upd_regression_somatic,
      Upd_count_transcripts, Upd_cpm, Upd_tpm_do_not_use_for_quantification,
      Upd_isoform_exonic_length,
      sctransform_quantile, supplemental_bulk_cpm, assay.data.sc, flybase.gtf,
      list(nos.1=batch_umap_nos.1, nos.2=batch_umap_nos.2, tj.1=batch_umap_tj.1, tj.2=batch_umap_tj.2) %>%
        bind_rows(.id = "batch"),
      'scRNA-seq-Regression/Enriched-Genes.xlsx'
    ),
    format='file',
    packages = tar_option_get("packages") %>% c("lme4")
  ),

  # ChIC paired-end alignment targets.
  tar_file(align_chic_lightfiltering, "scripts/align_chic_lightfiltering.sh"),
  tar_file(align_repli_lightfiltering, "scripts/align_repli_lightfiltering.sh"),
  tar_file(align_repli, "scripts/align_repli.sh"),
  chic.raw.tracks,

  tar_map(
    chic.fpkm.data,
    names = name,
    tar_target(
      bed,
      reference_sort_by_fpkm_table(
        Upd_cpm, tolower(name), assay.data.sc,
        paste0("scRNA-seq-Regression/", name, "-FPKM.bed")
      ),
      format = "file"
    ),
    tar_target(
      quartile.factor,
      quant_quartile_factor(Upd_cpm[, tolower(name)], q1_threshold=5),
      packages = tar_option_get("packages") %>% c("tidyr")
    )
  ),

  # TODO: Try replacing this target (data frame) with feature.rle
  tar_target(
    chic.genome,
    feature.lengths %>% data.frame(chr = names(.), len = .)
  ),
  tar_target(
    genomic_feature_factor,
    factor_genome(flybase.gtf, assay.data.sc, feature.lengths)
  ),

  # Description: Heatmaps of ChIC vs sorted genes.
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

  # tar_target(
  #   plot.chic.peak.location_Germline,
  #   display_peak_location_stats(
  #     list(
  #       H3K4=chic.peak.location.stat_H3K4_Germline,
  #       H3K27=chic.peak.location.stat_H3K27_Germline,
  #       H3K9=chic.peak.location.stat_H3K9_Germline
  #     )
  #   )
  # ),
  # tar_target(
  #   plot.chic.peak.location_Somatic,
  #   display_peak_location_stats(
  #     list(
  #       H3K4=chic.peak.location.stat_H3K4_Somatic,
  #       H3K27=chic.peak.location.stat_H3K27_Somatic,
  #       H3K9=chic.peak.location.stat_H3K9_Somatic
  #     )
  #   )
  # ),
  # tar_target(
  #   limits_plot.chic.peak.location,
  #   c(
  #     0,
  #     sapply(
  #       list(plot.chic.peak.location_Germline, plot.chic.peak.location_Somatic),
  #       \(pl) pl$data$value
  #     ) %>% max
  #   )
  # ),

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

  tar_target(demo.f.distribution, demo_f_distribution()),
  tar_target(plot.scaled.f.distribution, plot_scaled_f()),

  tar_target(
    cpm_gene_lists_extended,
    cpm_gene_lists %>%
      with(
        list(
          AllGermlineGenes = germline,
          AllSomaticGenes = somatic,
          ExclusiveGermlineGenes = setdiff(germline, somatic),
          ExclusiveSomaticGenes = setdiff(somatic, germline),
          NonExclusiveGermlineAndSomaticGenes = intersect(germline, somatic),
          AllGermlineOrSomaticGenes = union(germline, somatic),
          OffGenes = rownames(read.csv(assay.data.sc, row.names=1)) %>% setdiff(germline) %>% setdiff(somatic)
        )
    )
  ),
  tar_map(
    tibble(
      chic.fpkm.data,
      sc_extended_data_TSS = rlang::syms(
        str_glue("sc_extended_data_{name}_TSS")
      ),
      sc_extended_data_Paneled = rlang::syms(
        str_glue("sc_extended_data_{name}_Paneled")
      )
    ),
    names = name,
    tar_target(
      name = fpkm.chic.plots,
      tibble(
        experiment = name,
        gene_list = names(cpm_gene_lists_extended),
        tss_plot = list(
          sc_extended_data_TSS %>%
            subset(gene_list %in% c(names(cpm_gene_lists_extended), "OffGenes")) %>%
            mutate(
              genes = gene_list %>%
                list %>%
                append(setNames(c(names(cpm_gene_lists_extended), "OffGenes"), c("on", "off"))) %>%
                do.call(fct_recode, .) %>%
                fct_relevel("off", "on")
            ) %>%
            arrange(genes) %>%
            chic_plot_average_profiles_facet_grid(
              "",
              chic_line_track_colors[[tolower(name)]] %>%
                c(muted(., c=20, l=80), .),
              c(0.33, 0.66),
              facet_wrap(vars(mark))
            )
        ),
        paneled_plot = list(
          sc_extended_data_Paneled %>%
            subset(gene_list %in% c(names(cpm_gene_lists_extended), "OffGenes")) %>%
            mutate(
              genes = gene_list %>%
                list %>%
                append(setNames(c(names(cpm_gene_lists_extended), "OffGenes"), c("on", "off"))) %>%
                do.call(fct_recode, .) %>%
                fct_relevel("off", "on")
            ) %>%
            arrange(genes) %>%
            chic_plot_paneled_profiles_facet_grid(
              "",
              chic_line_track_colors[[tolower(name)]] %>%
                c(muted(., c=20, l=80), .),
              c(0.5, 1),
              facet_wrap(vars(mark))
            )
        )
      ),
      pattern = map(cpm_gene_lists_extended)
    ),
    tar_file(
      name = fig.fpkm.chic,
      save_figures(
        paste0("figure/", name),
        ".pdf",
        fpkm.chic.plots %>%
          filter(gene_list != "OffGenes") %>%
          rowwise %>%
          reframe(
            experiment,
            gene_list,
            prefix = c("CHIC-TSS-", "CHIC-"),
            width = c(8, 12),
            figure = list(tss_plot, paneled_plot)
          ) %>%
          reframe(
            filename = paste0(prefix, "AllMarks-RNAseq-CPM-", gene_list),
            figure,
            width,
            height = 4
          )
      )
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
          ),
          9,
          12,
          "CHIC-TSS-Chr-AllMarks-RNAseq",
          chic_plot_average_profiles_facet_grid(
            dplyr::rename(sc_chr_active_data_TSS, genes=activity),
            "CPM Quartile",
            rep(chic_line_track_colors[[tolower(name)]], 2)
          )
          + guides(linewidth = guide_none(), color = guide_none()),
          8,
          10
        )
      )
    )
  ),
  tar_map(
    tribble(
      ~celltype, ~sc_chr_nucleosome_data_TSS, ~granges,
      "Germline", rlang::sym("sc_chr_nucleosome_data_Germline_TSS"), rlang::sym("chic.experiment.quantify_H3K4_Germline_CN_chr"),
      "Somatic", rlang::sym("sc_chr_nucleosome_data_Somatic_TSS"), rlang::sym("chic.experiment.quantify_H3K4_Somatic_CN_chr")
    ),
    names = celltype,
    tar_target(
      fig.fpkm.chic.facet.nucleosome,
      save_figures(
        str_glue("figure/", celltype),
        ".pdf",
        tibble(
          rowname="CHIC-TSS-Chr-Nucleosome-Occupancy-RNAseq",
          figure=list(
            chic_plot_average_profiles_facet_grid(
              dplyr::rename(sc_chr_nucleosome_data_TSS, genes=activity),
              "",
              rep(chic_line_track_colors[[tolower(celltype)]], 2),
              faceter = facet_wrap(vars(facet)),
              x_intercept = NA
            ) + scale_y_continuous(
              name = 'Enrichment (vs Auto Monosome Median)',
              trans = 'log',
              limits = c(0.9, 32),
              expand = c(0, 0),
              breaks = c(1, 4, 16)
            ) + guides(linewidth = guide_none(), color = guide_none())
          ),
          width=6,
          height=6
        )
      )
    )
  ),
  tar_target(
    fig.fpkm.chic.facet.nucleosome.both,
    save_figures(
      "figure/Both-Cell-Types",
      ".pdf",
      tibble(
        rowname="CHIC-TSS-Chr-Nucleosome-Occupancy-RNAseq",
        figure=list(
          (
            chic_plot_average_profiles_facet_grid(
              list(Germline=sc_chr_nucleosome_data_Germline_TSS, Somatic=sc_chr_nucleosome_data_Somatic_TSS) %>%
                bind_rows(.id = "name") %>%
                arrange(desc(row_number())) %>%
                mutate(genes=interaction(activity, ordered(name)), facet=facet %>% ordered(c("X","2","3","4"))),
              "",
              c(muted(chic_line_track_colors$germline, l=70), chic_line_track_colors$germline, muted(chic_line_track_colors$somatic, l=70), chic_line_track_colors$somatic),
              linewidth = c(0.33, 0.66, 0.33, 0.66),
              faceter = facet_wrap(vars(facet), ncol=4),
              x_intercept = NA
            ) + scale_y_continuous(
              name = 'Enrichment (vs Auto Monosome Median)',
              trans = 'log',
              limits = c(0.9, 32),
              expand = c(0, 0),
              breaks = c(1, 4, 16)
            ) + theme(
              aspect.ratio = 1
            )
          ) %>%
            replace_legend(germline_somatic_line_plot_legend)
        ),
        width=10,
        height=3
      )
    )
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
          replace_legend(germline_somatic_line_plot_legend),
        12,
        8
      )
    )
  ),
  tar_file(
    fig.fpkm.chic.both.genelists,
    save_figures(
      "figure/Both-Cell-Types",
      ".pdf",
      tibble(
        fpkm.chic.plots_Germline %>% rename_with(\(n) paste0(n, "_Germline")),
        fpkm.chic.plots_Somatic %>% rename_with(\(n) paste0(n, "_Somatic"))
      ) %>%
        filter(gene_list_Germline != "OffGenes") %>%
        rowwise %>%
        reframe(
          gene_list = gene_list_Germline,
          tss_plot_data = bind_rows(list(Germline=tss_plot_Germline$data, Somatic=tss_plot_Somatic$data), .id="name") %>%
            mutate(genes = interaction(genes, name) %>% factor(c("off.Germline", "on.Germline", "off.Somatic", "on.Somatic"), ordered=T)) %>%
            arrange(name == "Germline") %>%
            list,
          paneled_plot_data = bind_rows(list(Germline=paneled_plot_Germline$data, Somatic=paneled_plot_Somatic$data), .id="name") %>%
            mutate(genes = interaction(genes, name) %>% factor(c("off.Germline", "on.Germline", "off.Somatic", "on.Somatic"), ordered=T)) %>%
            arrange(name == "Germline") %>%
            list
        ) %>%
        rowwise %>%
        reframe(
          rowname = str_glue("CHIC-{c('TSS-','')}AllMarks-RNAseq-CPM-{gene_list}") %>% as.character,
          figure = list(
            chic_plot_average_profiles_facet_grid(
              tss_plot_data,
              "",
              sapply(
                chic_line_track_colors,
                \(col) col %>%
                c(muted(., c=20, l=80), .)
              ) %>%
                as.character,
              c(0.33, 0.66, 0.33, 0.66),
              facet_wrap(vars(mark))
            ) %>%
              replace_legend(germline_somatic_line_plot_legend),
            chic_plot_paneled_profiles_facet_grid(
              paneled_plot_data,
              "",
              sapply(
                chic_line_track_colors,
                \(col) col %>%
                c(muted(., c=20, l=80), .)
              ) %>%
                as.character,
              c(0.33, 0.66, 0.33, 0.66),
              facet_wrap(vars(mark))
            ) %>%
              replace_legend(germline_somatic_line_plot_legend)
          ),
          width = c(9, 15),
          height = 4
        )
    )
  ),

  tar_map(
    tibble(
      chic.fpkm.data,
      figure_dir = str_glue("figure/{name}"),
      plot.chic.peak.location = rlang::syms(str_glue("plot.chic.peak.location_{name}")),
      sc_quartile_data_TSS = rlang::syms(str_glue("sc_quartile_data_{name}_TSS")),
      sc_quartile_data_Paneled = rlang::syms(str_glue("sc_quartile_data_{name}_Paneled"))
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
            facet_wrap(vars(mark))
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
          15, 4
          # "CHIC-AllMarks-Peak-Annotation",
          # plot.chic.peak.location
          # + scale_y_continuous(
          #   limits = limits_plot.chic.peak.location,
          #   expand = expansion(mult = c(0, 0.05))
          # ),
          # 6, 3
        )
      )
    )
  ),

  # Did not make one of these for the completed replication timing values and
  # their cuts yet. Do so!
  # tar_map(
  #   data.frame(extension = c(".pdf", ".png")),
  #   tar_target(
  #     name = repli.chic.quarter_Somatic,
  #     chic_average_profiles(
  #       repli.quarters_Tj_Weighted %>%
  #         rownames_to_column %>%
  #         pull(quarter, rowname) %>%
  #         # Some genes were not quantified in 10X Genomics BAM TX. The reason why
  #         # those genomic features were unsuitable for quantifying at the isoform
  #         # level by 10X Cell Ranger is unclear. To make Repli and FPKM somewhat
  #         # more comparable (although there may still be genomic features with low
  #         # read count in the Repli window and we will remove those features
  #         # before any Repli computation), we will filter Repli features using 10X
  #         # FPKM feature criteria.
  #         subset(names(.) %in% names(quartile.factor_Somatic)),
  #       dirname(
  #         c(
  #           chic.bw_H3K4_Somatic,
  #           chic.bw_H3K27_Somatic,
  #           chic.bw_H3K9_Somatic
  #         )[1]
  #       ),
  #       assay.data.sc,
  #       'tj',
  #       'Repli Quartile',
  #       setNames(repli_quartile_fills, NULL)
  #     ) %>%
  #       list %>%
  #       tibble(
  #         name = "Somatic_marks",
  #         figure = .,
  #         width = 9,
  #         height = 4
  #       ) %>%
  #       save_figures("repli/profile", extension, ., dpi = 300),
  #     format = "file"
  #   )
  # ),

  #,

  # tar_target(repli_dir, 'Repli-Seq/Validation', format='file'),
  # tar_target(name = repli.raw, command = load_repli(repli_dir), cue=tar_cue('never')),
  # tar_target(name = repli, command = normalize_repli(repli.raw), cue=tar_cue('never')),
  # tar_target(name = repli.hmm, command = repli_fit_hmm(repli), cue=tar_cue('never'))

  # Repli-Seq
  tar_file(run_fastqc_sh, "scripts/run_fastqc.sh"),

  targets.bulk.samples,
  targets.chic,
  targets.flybase,
  targets.quantification,
  targets.repli,
  targets.sce
)
 