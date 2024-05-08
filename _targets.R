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

# tar_make_future() is an older (pre-{crew}) way to do distributed computing
# in {targets}, and its configuration for your machine is below.
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

sce.data = tibble(
  batch=c('nos.1','nos.2','tj.1','tj.2'),
  capitalized=c('Nos','Nos','tj','tj'),
  batch_num=c('1','2','1','2'),
  dir_name=str_glue("scRNA-seq/{capitalized}-Upd_H3-GFP_Rep{batch_num}"),
  tenx_path=str_glue("{dir_name}/outs/filtered_feature_bc_matrix"),
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
  chic.samples,
  names = sample,
  unlist = FALSE,
  tar_file(
    chic.bam,
    tibble(
      tmp_name = tempfile(pattern = sample, fileext = ".bam"),
      filtered_name = paste0("chic/", group, "/", sample, ".bam"),
      align_lightfiltering = run(
        "bash",
        c(
          "-i",
          align_chic_lightfiltering,
          flybase.bowtie %>% paste("chic_bowtie2", sep="/"),
          paste0(batch, "/", sample, "_R1_001.fastq.gz"),
          paste0(batch, "/", sample, "_R2_001.fastq.gz"),
          tmp_name
        )
      ) %>%
        list,
      chromatin_specific = run(
        "bash",
        c(
          "-i",
          align_chic_chromatin_specific_filtering,
          tmp_name,
          filtered_name
        )
      ) %>%
        list,
      move_log = file.copy(
        str_replace(tmp_name, ".bam", ".log"),
        str_replace(filtered_name, ".bam", ".log")
      ),
      move_markdup_log = file.copy(
        str_replace(tmp_name, ".bam", ".markdup.log"),
        str_replace(filtered_name, ".bam", ".markdup.log")
      ),
      clean_up_large_tmp_file = file.remove(tmp_name)
    ) %>%
      pull(filtered_name),
    cue = tar_cue("never")
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

repli.samples = read.csv('repli/repli_samples.csv')
repli.samples$replication_value = repli.samples$replication_value %>% factor(unique(.))
repli.samples$full = repli.samples$replication_value
levels(repli.samples$full) <- c("Early", "Early-Mid", "Mid-Late", "Late")
repli.samples$abbrev = repli.samples$replication_value
levels(repli.samples$abbrev) <- c("E", "G", "J", "L")
repli.samples <- repli.samples %>%
  mutate(name = paste0(genotype, "_", abbrev, rep))

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
      quantify_quantiles(metadata)
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
  ),

  # Consider that all of these targets can be replaced with one call to a new
  # perl script per batch. We are intentionally streaming down the bam file as
  # that may be faster than the disk, and we only need to do that once and split
  # the outputs. Or instead of a single-threaded script, we can generate a tee
  # command-line that uses process substitution to read the SAM output in
  # several processes.
  tar_map(
    unlist = FALSE,
    sce.clusters,
    names = cluster,
    tar_target(
      pseudobulk.tx,
      download_bam_transcripts(
        download_bam_tx_sh,
        metadata,
        batch,
        cluster,
        paste0('scRNA-seq-Regression/pseudobulk_', batch, '_', cluster, '.txt')
      ),
      format = 'file',
      # Remove tar_cue if the Upd_sc object changes. WHY do these targets keep
      # re-running so often when we expect no real changes to the Upd_sc and its
      # metadata csv?
      cue = tar_cue("never")
    ),
    tar_target(
      pseudobulk.library.size,
      system2(
        c(
          estimate_library_size_bam_tx_sh,
          metadata,
          batch,
          cluster
        ),
        stdout=TRUE
      ) %>%
        as.numeric,
      # Remove tar_cue if the Upd_sc object changes. WHY do these targets keep
      # re-running so often when we expect no real changes to the Upd_sc and its
      # metadata csv?
      cue = tar_cue("never")
    )
  )
) %>%
  append(
    list(
      # Map over clusters (below, from metadata) and extract pseudobulk counts.
      tar_file(download_bam_tx_sh, "scripts/download_bam_tx.sh"),
      tar_file(estimate_library_size_bam_tx_sh, "scripts/estimate_library_size_bam_tx.sh")
  )
)

repli_targets <- tar_map(
  repli.samples,
  names = name,
  tar_file(
    repli.fastqc,
    tibble(
      input_path = paste0("Upd_Tumor/Repli/", filename),
      output_dir = "repli/fastqc",
      output_dir_create = dir.create(output_dir, rec=TRUE, showW=FALSE),
      output_path = paste0(
        output_dir, "/", str_replace(filename, ".fastq.gz", "_fastqc.zip")
      ),
      run_fastqc = processx::run(run_fastqc_sh, c(input_path, output_path)) %>% list
    ) %>%
      pull(output_path)
  ),
  tar_file(
    repli.bam,
    tibble(
      output_path = paste0("repli/bam/", str_replace(filename, ".fastq.gz", ".bam")),
      align = run(
        "bash",
        c(
          "-i",
          align_repli,
          flybase.bowtie %>% paste("chic_bowtie2", sep="/"),
          paste0("Upd_Tumor/Repli/", filename),
          output_path
        )
      ) %>%
        list
    ) %>%
      pull(output_path)
  )
) %>%
  list(
    sapply(
      c("nos", "tj"),
      \(suffix) tar_target_raw(
        paste0("repli.bams_", suffix),
        subset(repli.samples, genotype == suffix) %>%
          with(
            call(
              "tibble",
              condition = as.character(abbrev),
              rep = rep,
              name = name,
              source_file = "repli.bam_" %>%
                paste0(name) %>%
                rlang::syms() %>%
                append("c", .) %>%
                do.call(call, ., quote=T)
            )
          )
      )
    ),
    tar_map(
      tribble(~suffix, ~my_table, "nos", quote(repli.bams_nos), "tj", quote(repli.bams_tj)),
      names = suffix,
      tar_target(
        repli.coverage,
        read_replicated_coverage(
          my_table,
          flybase.lengths,
          feature.lengths
        ),
        packages = "GenomicAlignments"
      )
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
    read.table(flybase.genome.index, header=F) %>%
      reframe(name=V1, value=V2) %>%
      deframe
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
  tar_map(
    tibble(seurat = rlang::syms(c("nos.1", "nos.2", "tj.1", "tj.2"))),
    tar_target(cluster.names, levels(Idents(seurat))),
    tar_target(
      cluster.deg.analysis,
      tibble(
        cluster = cluster.names,
        glm_fit = if (grepl("doublet|[0-9]", cluster.names))
          Upd_glm %>%
            small_cluster_test_apeglm(
              metadata, seurat, as.character(quote(seurat)), cluster.names,
              test_mle=FALSE
            ) %>%
            list
        else list(NA),
        marker_genes = if (is.list(glm_fit[[1]]))
          tibble(
            name = rownames(glm_fit[[1]]$map),
            `l2FC95-` = (
              glm_fit[[1]]$map[, "mmOtherOverMean"] - qnorm(0.975) * glm_fit[[1]]$sd[, "mmOtherOverMean"]
            ) / log(2),
            `l2FC95+` = (
              glm_fit[[1]]$map[, "mmOtherOverMean"] + qnorm(0.975) * glm_fit[[1]]$sd[, "mmOtherOverMean"]
            ) / log(2),
            lfsr = glm_fit[[1]]$fsr[, "mmOtherOverMean"]
          ) %>%
            # All marker genes will have a **** significance level.
            filter(lfsr < 0.0001) %>%
            arrange(desc(`l2FC95-`)) %>%
            head(100) %>%
            list
        else NULL
      ),
      pattern = map(cluster.names)
    )
  ),
  tar_target(Upd_sc, filter_integrate_data(list(nos.1,nos.2,tj.1,tj.2))),
  tar_target(Upd_model_matrix, build_model_matrix(FetchData(Upd_sc, c("ident", "batch")))),
  tar_target(
    shuffle_feature_plot,
    sample(Cells(Upd_sc))
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
          fpkm_quarter_density(log(Upd_cpm) / log(10), ylim = c(-2.75, 4.5))
          +
          annotate(
            "segment",
            -Inf, log(7) / log(10),
            xend = Inf, yend = log(7) / log(10),
            color = "#bb2e0b",
            linewidth = 0.75
          )
          + scale_y_continuous(
            breaks = seq(-2, 4)
          ),
          4, 4,
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
          4, 4
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
          "eya"
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
  tar_target(
    Upd_sc_size_factors,
    Upd_sc %>% pooled_size_factors_seurat_ident_autosome
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
    fit_glm(Upd_sc, Upd_model_matrix, metadata),
  ),
  tar_target(
    Upd_regression_somatic,
    apeglm_coef_table(Upd_glm)
  ),
  tar_target(
    Upd_regression_tid,
    apeglm_coef_table(Upd_glm, coef = 3)
  ),
  tar_target(
    Upd_regression_sompre,
    apeglm_coef_table(Upd_glm, coef = 4)
  ),
  tar_target(
    Upd_regression_mscl,
    apeglm_coef_table(Upd_glm, coef = 5)
  ),
  sce_targets,
  tar_combine(
    Upd_pseudobulk,
    sce_targets[["pseudobulk.tx_germline"]],
    sce_targets[["pseudobulk.tx_somatic"]],
    sce_targets[["pseudobulk.tx_spermatocyte"]],
    sce_targets[["pseudobulk.tx_somaticprecursor"]],
    sce_targets[["pseudobulk.tx_muscle"]],
    command = cbind(sce.data, tx_file = c(!!!.x))
  ),
  tar_combine(
    Upd_pseudobulk_sf,
    sce_targets[["pseudobulk.library.size_germline"]],
    sce_targets[["pseudobulk.library.size_somatic"]],
    sce_targets[["pseudobulk.library.size_spermatocyte"]],
    sce_targets[["pseudobulk.library.size_somaticprecursor"]],
    sce_targets[["pseudobulk.library.size_muscle"]],
    command = cbind(sce.data, size_factor = c(!!!.x))
  ),
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
      + geom_point(color = "#06470c", size = 3)
      + scale_x_continuous(breaks = c(1, 5, 10))
      + scale_y_continuous(
        trans = "sqrt", breaks = c(4, 36, 100, 150), limits = c(3, 155)
      )
      + labs(x = "Principal Component (Validation PCA)", y = "Explained Variance")
      + theme_cowplot()
  ),
  # For cell cycle scoring.
  tar_download(
    cell_cycle_scoring_human_supplemental,
    "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4481139/bin/NIHMS687993-supplement-supp_data_2.xlsx",
    "scRNA-seq/NIHMS687993-supplement-supp_data_2.xlsx"
  ),
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
    supplemental_pca_figure,
    {
      filename <- "figure/Integrated-scRNAseq/Germline-Somatic-Pairs.pdf"
      dir.create(dirname(filename), recursive = TRUE, showW = FALSE)
      CairoPDF(filename, width = 16, height = 9)
      print(plot_Upd_pca_components(Upd_sc, load_cell_cycle_score_drosophila(cell_cycle_drosophila, assay.data.sc)))
      dev.off()
      filename
    },
    format = "file"
  ),
  tar_target(
    sct_gene_lists,
    determine_gene_list(sctransform_quantile, supplemental_bulk_cpm)
  ),
  tar_target(
    sct_gene_venn,
    plot_gene_lists(
      sct_gene_lists
    )
  ),
  tar_target(
    cpm_gene_lists,
    apply(
      Upd_cpm[, c("germline", "somatic")],
      2,
      \(v) (v >= 7) %>% which %>% names,
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
          "RNAseq-Validation-Elbow", supplemental_elbow_figure, 5, 2.5,
          "RNAseq-CPM-Gene-List-Venn-Area-Blank", cpm_gene_venn, 3.6, 2.7,
          "RNAseq-SCT-Gene-List-Venn-Area-Blank", sct_gene_venn, 3.6, 2.7,
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
      Upd_cpm_transcripts, Upd_cpm_transcript_to_use, assay.data.sc
    )
  ),
  # TPM values based on exon length correction, then scale to sum to 1MM.
  tar_target(
    Upd_tpm_do_not_use_for_quantification,
    join_cpm_data(
      assay.data.sc, Upd_cpm_transcripts,
      Upd_cpm_transcript_to_use, Upd_fpkm_regression) %>%
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
      Upd_cpm_regression,
      corr=F) %>%
      table_to_tpm
  ),
  tar_target(
    sc_excel,
    publish_excel_results(
      Upd_regression_somatic, Upd_regression_tid,
      Upd_regression_sompre, Upd_regression_mscl,
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

  # Section: Parameter selection for omics data, tuning MAPQ threshold and
  # estimated fragment size using cross-correlation.
  tar_target(
    chic.pileup_tj,
    read_pileup_df("chic_tj_sample_reads.markdup.bam")
  ),
  tar_target(
    chic.coverage_tj,
    list(watson = mapq_table_strand(chic.pileup_tj, "A"), crick = mapq_table_strand(chic.pileup_tj, "a"))
  ),
  tar_target(
    chic.coverage.corr.plot_tj,
    analyze_mapq(chic.coverage_tj, test_mapq=c(0, 5, 10, 15, 20, 30))
  ),
  tar_target(
    chic.coverage.corr.plot.png_tj,
    ggsave('QC-ChIC-tj-Direction-Corr.png', chic.coverage.corr.plot_tj, width=8, height=6, dpi=80),
    format = "file"
  ),
  tar_target(
    repli.pileup_tj_202310,
    read_pileup_df("repli_tj_sample_reads.markdup.bam")
  ),
  tar_target(
    repli.coverage_tj_202310,
    list(watson = mapq_table_strand(repli.pileup_tj_202310, "A"), crick = mapq_table_strand(repli.pileup_tj_202310, "a"))
  ),
  tar_target(
    repli.coverage.corr.plot_tj_202310,
    analyze_mapq(repli.coverage_tj_202310)
  ),
  tar_target(
    repli.coverage.corr.plot.png_tj_202310,
    ggsave('QC-Repli-202310-tj-Corr.png', repli.coverage.corr.plot_tj_202310, width=8, height=6, dpi=80),
    format = "file"
  ),

  # ChIC paired-end alignment targets.
  tar_file(align_chic_lightfiltering, "scripts/align_chic_lightfiltering.sh"),
  tar_file(
    align_chic_chromatin_specific_filtering,
    "scripts/align_chic_chromatin_specific_filtering.sh"
  ),
  tar_file(align_repli, "scripts/align_repli.sh"),
  chic.raw.tracks,

  # ChIC coverage tracks.
  tar_target(
    chic.macs.bandwidth, list(500, 1000, 2000)
  ),
  # Background value for uniform FPKM. We have a multiplier of 1MM but this
  # amount is spread across the reference sequence length.
  tar_target(
    chic.macs.background.value,
    1000 * 1000 * 1000 / sum(feature.lengths)
  ),
  tar_target(
    chic.macs.background.track,
    feature.lengths %>%
      sapply(\(len) Rle(values = chic.macs.background.value, len)) %>%
      RleList
  ),
  tar_map(
    chic.lookup,
    names = lookup_name,
    tar_map(
      data.frame(bw = c(25, 125, 250)),
      tar_target(
        chic.smooth,
        setNames(samples, sample_names) %>%
          sapply(
            \(v) v %>%
              smooth_sparse_vector_to_rle_list(
                feature.rle,
                bw = bw,
                sample_size = ifelse(bw < 50, 10, 50)
              )
          ) %>%
          purrr::reduce(`+`) %>%
          `/`(length(samples)) %>%
          rle_list_round_log_scale_pretty %>%
          pmax(
            . * 0 + 0.1
          ) %>%
          set_attr("n", length(samples))
      ),
      tar_target(
        chic.sd,
        {
          s <- setNames(samples, sample_names) %>%
            sapply(
              \(v) v %>%
                smooth_sparse_vector_to_rle_list(
                  feature.rle,
                  bw = bw,
                  sample_size = ifelse(bw < 50, 10, 50)
                )
            ) %>%
            rle_lists_sd
          attr(s, "n") <- length(samples)
          s
        }
      )
    ),
    tar_target(
      chic.wide,
      samples %>%
        sapply(
          \(x) x %>%
            smooth_sparse_vector_to_rle_list(
              feature.rle,
              bw = 1000,
              sample_size = 200
            ),
          simplify = FALSE
        ) %>%
        # Arithmetic mean of reps
        purrr::reduce(`+`) %>%
        `/`(length(samples)),
    ),
    tar_target(chic.macs.tagcounts, samples %>% sapply(sum)),
    tar_target(
      chic.macs.byparam,
      # MACS adaptive parameter using a large window (chic.macs.bandwidth), now
      # with multiple replicates (samples).
      samples %>%
        sapply(
          \(x) x %>%
            smooth_sparse_vector_macs(
              feature.rle,
              chic.macs.bandwidth[[1]],
              sample_size = 200,
              normalize_tag_count = TRUE
            ),
          simplify = FALSE
        ) %>%
        # In MACS Poisson distribution fitting, we can take the arithmetic mean
        # of the parameter normalized by tag count, for the biological
        # replicates.
        # TODO: Optimize the parameter using glm instead of using arithmetic
        # mean. We have applied an offset (tag count - FPKM calculation), so
        # unlike the case with no offset, the optimum might be different than
        # arithmetic mean.
        purrr::reduce(`+`) %>%
        `/`(length(samples)),
      pattern = map(chic.macs.bandwidth),
      iteration = "list"
    ),
    tar_target(
      # MACS Poisson parameter is the max of the input track created by 3
      # different large sliding window sizes. We removed the uniform coverage
      # track.
      chic.macs.fpkm,
      purrr::reduce(chic.macs.byparam, pmax)
    ),
    tar_target(
      chic.macs.peakdetect,
      samples %>%
        sapply(
          \(x) x %>%
            smooth_sparse_vector_macs(
              feature.rle,
              diameter = 150,
              sample_size = 50,
              normalize_tag_count = FALSE
            ),
          simplify = FALSE
        )
    ),
    # Alternative pipeline merging BAM files before proceeding.
    tar_target(
      chic.merge.bam,
      {
        filename <- paste0(
          dirname(files[1]),
          "/",
          lookup_name,
          ".bam"
        )
        processx::run(
          "samtools",
          c(
            "merge",
            "-f",
            filename,
            files
          )
        )
        processx::run("samtools", c("index", filename))
        filename
      },
      format = "file",
      cue = tar_cue("never")
    )
  ),

  tar_target(
    chic.squeezeVar_Germline,
    list(
      H3K4=chic.sd_250_input_H3K4_Germline,
      H3K27=chic.sd_250_input_H3K27_Germline,
      H3K9=chic.sd_250_input_H3K9_Germline
    ) %>%
      chic_squeeze_var,
    packages = "limma"
  ),
  tar_target(
    chic.squeezeVar_Somatic,
    list(
      H3K4=chic.sd_250_input_H3K4_Somatic,
      H3K27=chic.sd_250_input_H3K27_Somatic,
      H3K9=chic.sd_250_input_H3K9_Somatic
    ) %>%
      chic_squeeze_var,
    packages = "limma"
  ),
  tar_target(
    chic.squeezeVar,
    list(Germline = chic.squeezeVar_Germline, Somatic = chic.squeezeVar_Somatic)
  ),

  tar_target(
    chic.macs.peak_H3K4_Germline,
    chic_macs_test(
      chic.macs.fpkm_input_H3K4_Germline,
      chic.wide_input_H3K4_Germline,
      chic.macs.peakdetect_mod_H3K4_Germline
    )
  ),
  tar_target(
    chic.macs.peak_H3K27_Germline,
    chic_macs_test(
      chic.macs.fpkm_input_H3K27_Germline,
      chic.wide_input_H3K27_Germline,
      chic.macs.peakdetect_mod_H3K27_Germline
    )
  ),
  tar_target(
    chic.macs.peak_H3K9_Germline,
    chic_macs_test(
      chic.macs.fpkm_input_H3K9_Germline,
      chic.wide_input_H3K9_Germline,
      chic.macs.peakdetect_mod_H3K9_Germline
    )
  ),
  tar_target(
    chic.macs.peak_H3K4_Somatic,
    chic_macs_test(
      chic.macs.fpkm_input_H3K4_Somatic,
      chic.wide_input_H3K4_Somatic,
      chic.macs.peakdetect_mod_H3K4_Somatic
    )
  ),
  tar_target(
    chic.macs.peak_H3K27_Somatic,
    chic_macs_test(
      chic.macs.fpkm_input_H3K27_Somatic,
      chic.wide_input_H3K27_Somatic,
      chic.macs.peakdetect_mod_H3K27_Somatic
    )
  ),
  tar_target(
    chic.macs.peak_H3K9_Somatic,
    chic_macs_test(
      chic.macs.fpkm_input_H3K9_Somatic,
      chic.wide_input_H3K9_Somatic,
      chic.macs.peakdetect_mod_H3K9_Somatic
    )
  ),

  tar_map(
    chic.fpkm.data,
    names = name,
    tar_target(
      bed,
      reference_sort_by_fpkm_table(
        Upd_fpkm, tolower(name), assay.data.sc,
        paste0("scRNA-seq-Regression/", name, "-FPKM.bed")
      ),
      format = "file"
    ),
    tar_target(
      quartile.factor,
      quant_quartile_factor(Upd_cpm[, tolower(name)], q1_threshold=7),
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

  tar_map(
    chic.experiments,
    names = experiment_name,

    # Molecule/input FPKM ratio track, in bigwig format.
    tar_target(
      chic.bw,
      {
        filename <- paste0('chic/', driver, '_', mark, '.FE.bw')
        export(
          chic_mod_sym / chic_input_sym
          * RleList(
            sapply(
              chic_input_sym,
              \(v) ifelse(v >= 1, 1, 0),
              simplify = FALSE
            )
          ),
          filename,
          'bigwig'
        )
        filename
      },
      format = 'file'
    ),
    tar_target(
      tss_mark_matrix,
      flybase_big_matrix(
        load_flybase_bed(bed_sym),
        chic.bw
      )
    ),
    tar_target(
      tss_mark_heatmap,
      display_tss_tile_matrix(
        tss_mark_matrix,
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
        tss_mark_matrix,
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
        tss_mark_matrix,
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
        tss_mark_matrix,
        paste0("figure/", name, "/INV-FPKM-", mark, ".png"),
        scale_image_filter = rep(1/50, 50),
        fc_max = 4,
        fc_filter = Inf,
        direction = -1
      ),
      format = "file"
    ),

    tar_target(
      bam_replicates_input,
      as.list(chic_input_files)
    ),
    tar_target(
      bam_replicates_mod,
      as.list(chic_mod_files)
    ),
    tar_target(
      chic.enrich.grid,
      enrichR(
        bam_replicates_mod[[1]],
        bam_replicates_input[[1]],
        chic.genome,
        countConfig = countConfigPairedEnd(tlenFilter = c(70L, 400L), midpoint = FALSE)
      ) %>% enrichr_set_names(
        bam_replicates_mod[[1]],
        bam_replicates_input[[1]]
      ),
      pattern = cross(bam_replicates_input, bam_replicates_mod),
      iteration = "list"
    ),
    tar_target(
      chic.input.count,
      as.numeric(
        processx::run(
          "samtools",
          c("view", "-c", bam_replicates_input[[1]])
        )$stdout
      ),
      pattern = map(bam_replicates_input),
      iteration = "list"
    ),
    tar_target(
      chic.mod.count,
      as.numeric(
        processx::run(
          "samtools",
          c("view", "-c", bam_replicates_mod[[1]])
        )$stdout
      ),
      pattern = map(bam_replicates_mod),
      iteration = "list"
    ),
    tar_target(
      plot.chic.multiplier,
      enrichr_grid_ratio(
        chic.enrich.grid,
        setNames(unlist(chic.input.count), unlist(bam_replicates_input)),
        setNames(unlist(chic.mod.count), unlist(bam_replicates_mod))
      )
    ),
    tar_target(
      enrichr.size.factor,
      enrichr_ratio_per_reads_mapped_adjustment(plot.chic.multiplier)
    ),
    tar_target(
      plot.chic.multiplier_png,
      save_figures(
        paste0("figure/", name), ".png",
        tribble(
          ~name, ~figure, ~width, ~height,
          paste0("EnrichR-Norm-", name, "-", mark),
          plot.chic.multiplier,
          6,
          4
        )
      )
    ),

    tar_target(
      chic.var.limma,
      chic.squeezeVar[[name]]$var[[mark]]
    ),
    tar_target(
      chic.test.limma,
      chic_ttest(
        pmax(chic_smooth_125_mod_sym, chic_smooth_250_mod_sym) %>%
          set_attr("n", attr(chic_smooth_125_mod_sym, "n")),
        chic_smooth_250_input_sym,
        sd = sqrt(chic.var.limma),
        df = chic.squeezeVar[[name]]$df
      )
    ),
    tar_target(
      chic.broad.peaks.stat,
      chic_quantify_broad_peaks(
        pmax(chic_smooth_125_mod_sym, chic_smooth_250_mod_sym)
        / chic_smooth_250_input_sym,
        track_mask = chic_smooth_250_input_sym >= 1,
        features = read.csv(assay.data.sc, row.names = 1)
      )
    ),
    tar_target(
      chic.test.limma.bed,
      chic.test.limma %>%
        chic_track_generate_table_by_enrichment %>%
        subset(q < 0.05)
    ),
    tar_target(
      chic.peak.location.stat,
      classify_feature_overlaps(chic.test.limma.bed, genomic_feature_factor)
    )
  ),

  tar_target(
    plot.chic.peak.location_Germline,
    display_peak_location_stats(
      list(
        H3K4=chic.peak.location.stat_H3K4_Germline,
        H3K27=chic.peak.location.stat_H3K27_Germline,
        H3K9=chic.peak.location.stat_H3K9_Germline
      )
    )
  ),
  tar_target(
    plot.chic.peak.location_Somatic,
    display_peak_location_stats(
      list(
        H3K4=chic.peak.location.stat_H3K4_Somatic,
        H3K27=chic.peak.location.stat_H3K27_Somatic,
        H3K9=chic.peak.location.stat_H3K9_Somatic
      )
    )
  ),
  tar_target(
    limits_plot.chic.peak.location,
    c(
      0,
      sapply(
        list(plot.chic.peak.location_Germline, plot.chic.peak.location_Somatic),
        \(pl) pl$data$value
      ) %>% max
    )
  ),

  tar_target(
    plot.chic.input.anova_H3K4_Germline,
    plot_chic_anova(
      chic.smooth_125_mod_H3K4_Germline %>%
        set_attr("standard_deviation", chic.sd_125_mod_H3K4_Germline),
      list(
        H3K4=chic.smooth_125_input_H3K4_Germline %>%
          set_attr("standard_deviation", chic.sd_125_input_H3K4_Germline),
        H3K27=chic.smooth_125_input_H3K27_Germline %>%
          set_attr("standard_deviation", chic.sd_125_input_H3K27_Germline),
        H3K9=chic.smooth_125_input_H3K9_Germline %>%
          set_attr("standard_deviation", chic.sd_125_input_H3K9_Germline)
      ),
      bin_size = 50
    )
  ),
  tar_target(
    plot.chic.input.anova_H3K4_Somatic,
    plot_chic_anova(
      chic.smooth_125_mod_H3K4_Somatic %>%
        set_attr("standard_deviation", chic.sd_125_mod_H3K4_Somatic),
      list(
        H3K4=chic.smooth_125_input_H3K4_Somatic %>%
          set_attr("standard_deviation", chic.sd_125_input_H3K4_Somatic),
        H3K27=chic.smooth_125_input_H3K27_Somatic %>%
          set_attr("standard_deviation", chic.sd_125_input_H3K27_Somatic),
        H3K9=chic.smooth_125_input_H3K9_Somatic %>%
          set_attr("standard_deviation", chic.sd_125_input_H3K9_Somatic)
      ),
      bin_size = 50
    )
  ),

  tar_map(
    data.frame(macs_background = c(500, 1000, 5000, 10000)),
    tar_target(
      chic_macs_mean_variance,
      chic_illustrate_mean_variance(
        list(
          chic.raw_GC2894001_S1_L001, chic.raw_GC2894002_S2_L001,
          chic.raw_GC2894003_S3_L001
        ),
        macs_background
      ),
    )
  ),
  tar_target(
    chic_macs_mean_variance_s2_files,
    c("chip_seq_s2/SRR067914.bw", "chip_seq_s2/SRR067915.bw"),
    format = "file"
  ),
  tar_target(
    chic_macs_mean_variance_s2,
    chic_illustrate_s2_mean_variance(chic_macs_mean_variance_s2_files)
  ),
  tar_target(
    chic_macs_mean_variance_gaussian,
    chic_illustrate_mean_variance_gaussian(
      list(
        chic.raw_GC2894001_S1_L001, chic.raw_GC2894002_S2_L001,
        chic.raw_GC2894003_S3_L001
      ),
      bw = 250
    )
    + coord_cartesian(c(0,12), c(0,30), expand=F)
  ),

  tar_target(
    chic.h3.tj.track.plot,
    plot_chic_rect_window_replicates(
      c(
        chic.macs.peakdetect_input_H3K4_Germline,
        chic.macs.peakdetect_input_H3K27_Germline,
        chic.macs.peakdetect_input_H3K9_Germline
      ),
      "2L", 19462000:19468800
    ) %>% ggplot(aes(x, y))
      + geom_line(aes(group=sample, color=sample), alpha=0.5)
      + geom_smooth(method="loess", span=0.2)
      + scale_color_manual(values=c("purple", "forestgreen", "cyan", "goldenrod", "red", "limegreen", "violet", "steelblue", "brown"), guide=NULL)
      + coord_cartesian(NULL, c(-2, 60), expand=F)
      + scale_x_continuous(labels=\(x) paste0(x / 1000 / 1000, " Mb"))
      + labs(x = NULL, y = "count")
      + theme_bw()
      + theme(aspect.ratio = 0.33)
  ),
  tar_target(
    chic.h3.histone.track.plot,
    plot_chic_rect_window_replicates(
      c(
        chic.macs.peakdetect_input_H3K4_Germline,
        chic.macs.peakdetect_input_H3K27_Germline,
        chic.macs.peakdetect_input_H3K9_Germline
      ),
      "2L",
      21398000:21549000,
      by = 100
    ) %>% ggplot(aes(x, y))
      + geom_line(aes(group=sample, color=sample), alpha=0.5)
      + geom_smooth(method="loess", span=0.2)
      + scale_color_manual(values=c("purple", "forestgreen", "cyan", "goldenrod", "red", "limegreen", "violet", "steelblue", "brown"), guide=NULL)
      + coord_cartesian(NULL, c(-2, 10000), expand=F)
      + scale_x_continuous(labels=\(x) paste0(x / 1000 / 1000, " Mb"))
      + labs(x = NULL, y = "count")
      + theme_bw()
      + theme(aspect.ratio = 0.33)
  ), tar_target(
    chic.h3.histone.track.plot.fpkm,
    plot_chic_rect_window_replicates(
      c(
        chic.macs.peakdetect_input_H3K4_Germline,
        chic.macs.peakdetect_input_H3K27_Germline,
        chic.macs.peakdetect_input_H3K9_Germline
      ) %>%
        sapply(\(rl) rl * 1000 * 1000 * 1000 / 150 / attr(rl, "num_tags")),
      "2L",
      21398000:21549000,
      by = 100
    ) %>% ggplot(aes(x, y))
      + geom_line(aes(group=sample, color=sample), alpha=0.5)
      + geom_smooth(method="loess", span=0.2)
      + scale_color_manual(values=c("purple", "forestgreen", "cyan", "goldenrod", "red", "limegreen", "violet", "steelblue", "brown"), guide=NULL)
      + coord_cartesian(NULL, c(-2, 1000), expand=F)
      + scale_x_continuous(labels=\(x) paste0(x / 1000 / 1000, " Mb"))
      + labs(x = NULL, y = "count")
      + theme_bw()
      + theme(aspect.ratio = 0.33)
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

  tar_target(demo.f.distribution, demo_f_distribution()),
  tar_target(plot.scaled.f.distribution, plot_scaled_f()),

  tar_map(
    data.frame(extension = c(".pdf", ".png")),
    tar_target(
      chic_poisson_illustration,
      save_figures(
        "figure/Germline", extension,
        tribble(
          ~name, ~figure, ~width, ~height,
          "CHIC-H3K4-Sample-Simulation-1", 
          illustrate_coverage_poisson(
            # Center on tj gene
            chic.smooth_25_mod_H3K4_Germline$`2L`[seq(19463500,19466500,by=10)] %>%
              # Square the track (shrinks the variations where the value is smaller)
              `^`(2) %>%
              # Subtract off the minimum value
              `-`(min(.) * 0.75) %>%
              # rle to vector
              as.numeric
          ),
          4,
          4,
          "CHIC-H3K4-Sample-Simulation-2",
          illustrate_poisson_variable(
            # Center on tj gene
            chic.smooth_25_mod_H3K4_Germline$`2L`[seq(19463500,19466500,by=10)] %>%
              # Square the track (shrinks the variations where the value is smaller)
              `^`(2) %>%
              # Subtract off the minimum value
              `-`(min(.) * 0.75) %>%
              # rle to vector
              as.numeric
          ),
          4,
          4,
          "CHIC-H3K4-Mean-Variance-Window-500",
          chic_macs_mean_variance_500,
          4,
          8,
          "CHIC-H3K4-Mean-Variance-Window-1K",
          chic_macs_mean_variance_1000
          + coord_cartesian(c(0,350), c(0,25000), expand=FALSE),
          4,
          8,
          "CHIC-H3K4-Mean-Variance-Window-5K",
          chic_macs_mean_variance_5000,
          4,
          8,
          "CHIC-S2-H3-Mean-Variance-Window-1K",
          chic_macs_mean_variance_s2,
          4,
          8,
          "CHIC-H3K4-Mean-Variance-Gaussian",
          chic_macs_mean_variance_gaussian,
          4,
          8,
          "CHIC-H3-TJ",
          chic.h3.tj.track.plot,
          6,
          2,
          "CHIC-F-Test-Gaussian",
          plot.chic.input.anova_H3K4_Germline
          + theme(aspect.ratio = 1),
          6, 4,
          "F-Statistic",
          (demo.f.distribution
          + annotate(
            "text", -3, 0.5, label = bquote(s[X]^2*"  = SSB"), color = "red"
          )
          + annotate(
            "text", 3, -6.5, label = bquote(s[Y]^2*"  = SSB"), color = "red"
          )
          + annotate(
            "text", -3, -6.5, label = bquote("F =  "*s[Y]^2*"/"*s[X]^2)
          )
          )
          ,
          3, 3,
          "F-Scaled-Distribution",
          plot.scaled.f.distribution,
          8, 6
        )
      )
    ),
    tar_target(
      chic_poisson_illustration_somatic,
      save_figures(
        "figure/Somatic", extension,
        tribble(
          ~name, ~figure, ~width, ~height,
          "CHIC-F-Test-Gaussian",
          plot.chic.input.anova_H3K4_Somatic
          + theme(aspect.ratio = 1),
          6, 4
        )
      )
    ),
    tar_target(
      supplemental_chic_model,
      save_figures(
        "chic/Model-Selection", extension,
        tribble(
          ~name, ~figure, ~width, ~height,
          "Germline-F-Var",
          plot.chic.input.anova_H3K4_Germline
          + labs(
            title = "ChIC Germline Equality of Variance",
            subtitle = "H3K4 variance following null hypothesis is uncertain (dotted)"
          ) + theme(aspect.ratio = 1),
          6, 4,
          "Somatic-F-Var",
          plot.chic.input.anova_H3K4_Somatic
          + labs(
            title = "ChIC Somatic Equality of Variance",
            subtitle = "H3K4 variance following null hypothesis is uncertain (dotted)"
          ) + theme(aspect.ratio = 1),
          6, 4
        )
      ),
      format = "file"
    )
  ),

  tar_target(
    chic.peaks.bed_Germline,
    write_chic_peaks(
      list(
        H3K4=chic.test.limma_H3K4_Germline %>% chic_track_generate_table_by_enrichment %>% subset(q < 0.10),
        H3K9=chic.test.limma_H3K9_Germline %>% chic_track_generate_table_by_enrichment %>% subset(q < 0.10),
        H3K27=chic.test.limma_H3K27_Germline %>% chic_track_generate_table_by_enrichment %>% subset(q < 0.10)
      ),
      "chic/Germline-Peaks.bed"
    ),
    format = "file"
  ),
  tar_target(
    chic.peaks.bed_Somatic,
    write_chic_peaks(
      list(
        H3K4=chic.test.limma_H3K4_Somatic %>% chic_track_generate_table_by_enrichment %>% subset(q < 0.10),
        H3K9=chic.test.limma_H3K9_Somatic %>% chic_track_generate_table_by_enrichment %>% subset(q < 0.10),
        H3K27=chic.test.limma_H3K27_Somatic %>% chic_track_generate_table_by_enrichment %>% subset(q < 0.10)
      ),
      "chic/Somatic-Peaks.bed"
    ),
    format = "file"
  ),

  tar_target(
    chic.macs.bed_Germline,
    write_chic_peaks(
      list(
        H3K4=chic.macs.peak_H3K4_Germline %>% chic_macs_generate_table,
        H3K27=chic.macs.peak_H3K27_Germline %>% chic_macs_generate_table,
        H3K9=chic.macs.peak_H3K9_Germline %>% chic_macs_generate_table(1e-3)
      ),
      "chic/Germline-Peaks-MACS.bed"
    ),
    format = "file"
  ),
  tar_target(
    chic.macs.bed_Somatic,
    write_chic_peaks(
      list(
        H3K4=chic.macs.peak_H3K4_Somatic %>% chic_macs_generate_table,
        H3K27=chic.macs.peak_H3K27_Somatic %>% chic_macs_generate_table,
        H3K9=chic.macs.peak_H3K9_Somatic %>% chic_macs_generate_table(1e-3)
      ),
      "chic/Somatic-Peaks-MACS.bed"
    ),
    format = "file"
  ),

  tar_target(chic.1kb, window_chic('chic', 1000)),

  tar_map(
    values = repli.processing,
    tar_target(repli.raw, load_repli('repli/bedgraph.in', suffix=name)),
    tar_target(repli, pseudobulk_repli(repli.raw))
  ),
  tar_target(
    repli.output.track,
    write_repli_bedgraph_out(repli_.q20.near10.bedgraph_DeleteNearbyPositions, 'repli/bedgraph.out')
  ),
  tar_target(
    repli.pct.track,
    write_repli_pct(repli_.q20.near10.bedgraph_DeleteNearbyPositions, 'repli/bedgraph.out')
  ),
  tar_target(
    repli_qc,
    repli_qc_make_figures(
      list(
        DupePosition=repli_.q20.markdup.bedgraph_DupePosition,
        DownsampleAndDupePosition=repli_.q20.p10.markdup.bedgraph_DownsampleAndDupePosition,
        DeleteNearbyPositions=repli_.q20.near10.bedgraph_DeleteNearbyPositions,
        OrigPileup=repli_.q20.bedgraph_OrigPileup # ,
        # DeleteNearbyPositions25=repli_.q20.near25.bedgraph_DeleteNearbyPositions25
      ),
      output_path='repli/supp_figures'
    ),
    format='file'
  ),
  tar_target(name = repli.hmm, command = repli_fit_hmm(repli_.q20.near10.bedgraph_DeleteNearbyPositions)),
  tar_target(repli.pct.bp, make_tj_bp_early_score('repli/bedgraph.in')),

  tar_map(
    values = repli.coverage,
    tar_target(
      repli.hdf5,
      repli_coverage_contrast(
        'repli/bedgraph.in',
        paste0('repli/hdf5/', ident, '_', contrast, '.h5'),
        ident,
        repli.coverage.contrasts[[contrast]],
        num_sorted_bases = 10000
      ),
      format = "file",
      cue = tar_cue("never")
    ),
    tar_target(
      repli.quarters,
      analyze_repli_quarters(repli.hdf5, assay.data.sc, flybase.gtf)
    )
  ),

  tar_target(
    name = repli.chic.tj,
    analyze_repli_chic(
      repli.hdf5_Tj_Weighted,
      'chic',
      flybase.gtf,
      'tj',
      'repli/heatmap/tj',
      num_sorted_bases = 10000
    ),
    format = 'file'
  ),

  tar_target(
    name = fpkm.chic.quarter.tss.plot_Germline,
    chic_average_profiles(
      quartile.factor_Germline,
      dirname(
        c(
          chic.bw_H3K4_Germline,
          chic.bw_H3K27_Germline,
          chic.bw_H3K9_Germline
        )[1]
      ),
      assay.data.sc,
      'Nos',
      'FPKM Quartile',
      setNames(sc_quartile_annotations, NULL)
    )
  ),
  tar_target(
    name = fpkm.chic.quarter.tss.plot_Somatic,
    chic_average_profiles(
      quartile.factor_Somatic,
      dirname(
        c(
          chic.bw_H3K4_Somatic,
          chic.bw_H3K27_Somatic,
          chic.bw_H3K9_Somatic
        )[1]
      ),
      assay.data.sc,
      'tj',
      'FPKM Quartile',
      setNames(sc_quartile_annotations, NULL)
    )
  ),
  tar_target(
    name = fpkm.chic.quarter.plot_Germline,
    chic_quartile_gene_list_paneled_profiles(
      quartile.factor_Germline,
      dirname(
        c(
          chic.bw_H3K4_Germline,
          chic.bw_H3K27_Germline,
          chic.bw_H3K9_Germline
        )[1]
      ),
      assay.data.sc,
      'Nos',
      'FPKM Quartile',
      setNames(sc_quartile_annotations, NULL)
    )
  ),
  tar_target(
    name = fpkm.chic.quarter.plot_Somatic,
    chic_quartile_gene_list_paneled_profiles(
      quartile.factor_Somatic,
      dirname(
        c(
          chic.bw_H3K4_Somatic,
          chic.bw_H3K27_Somatic,
          chic.bw_H3K9_Somatic
        )[1]
      ),
      assay.data.sc,
      'tj',
      'FPKM Quartile',
      setNames(sc_quartile_annotations, NULL)
    )
  ),

  tar_target(
    sct_gene_lists_extended,
    sct_gene_lists %>%
      with(
        list(
          AllGermlineGenes = germline,
          AllSomaticGenes = somatic,
          ExclusiveGermlineGenes = setdiff(germline, somatic),
          ExclusiveSomaticGenes = setdiff(somatic, germline),
          AllGermlineOrSomaticGenes = union(germline, somatic),
          OffGenes = rownames(read.csv(assay.data.sc, row.names=1)) %>% setdiff(germline) %>% setdiff(somatic)
        )
    )
  ),
  tar_target(sct_OffGenes, sct_gene_lists_extended$OffGenes),
  tar_target(
    cpm_gene_lists_extended,
    cpm_gene_lists %>%
      with(
        list(
          AllGermlineGenes = germline,
          AllSomaticGenes = somatic,
          ExclusiveGermlineGenes = setdiff(germline, somatic),
          ExclusiveSomaticGenes = setdiff(somatic, germline),
          AllGermlineOrSomaticGenes = union(germline, somatic),
          OffGenes = rownames(read.csv(assay.data.sc, row.names=1)) %>% setdiff(germline) %>% setdiff(somatic)
        )
    )
  ),
  tar_target(cpm_OffGenes, cpm_gene_lists_extended$OffGenes),
  tar_map(
    cross_join(
      tribble(
        ~genelist, ~gene_obj, ~OffGenes,
        "SCT", rlang::sym("sct_gene_lists_extended"), rlang::sym("sct_OffGenes"),
        "CPM", rlang::sym("cpm_gene_lists_extended"), rlang::sym("cpm_OffGenes")
      ),
      chic.fpkm.data %>%
        mutate(driver = driver %>% replace(. == "Nos", "nos"))
    ),
    names = genelist | name,
    tar_target(
      name = fpkm.chic.plots,
      list(
        # Fix deps on each bw file in this folder...
        dirname = "chic"
      ) %>%
        with(
          tibble(
            experiment = name,
            gene_list = names(gene_obj),
            tss_plot = chic_average_gene_list_profiles(
              gene_obj[[1]],
              OffGenes,
              dirname,
              assay.data.sc,
              driver,
              track_color = if (grepl("Germline", names(gene_obj)))
                chic_line_track_colors$germline
              else if (grepl("Somatic", names(gene_obj)))
                chic_line_track_colors$somatic
            ) %>%
              list,
            paneled_plot = chic_custom_gene_list_paneled_profile(
              gene_obj[[1]],
              OffGenes,
              dirname,
              assay.data.sc,
              driver,
              track_color = if (grepl("Germline", names(gene_obj)))
                chic_line_track_colors$germline
              else if (grepl("Somatic", names(gene_obj)))
                chic_line_track_colors$somatic
            ) %>%
              list
          )
        ),
      pattern = map(gene_obj)
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
            filename = paste0(prefix, "AllMarks-RNAseq-", genelist, "-", gene_list),
            figure,
            width,
            height = 4
          )
      )
    )
  ),

  tar_map(
    data.frame(extension = c(".pdf", ".png")),
    tar_target(
      name = repli.chic.quarter_Somatic,
      chic_average_profiles(
        repli.quarters_Tj_Weighted %>%
          rownames_to_column %>%
          pull(quarter, rowname) %>%
          # Some genes were not quantified in 10X Genomics BAM TX. The reason why
          # those genomic features were unsuitable for quantifying at the isoform
          # level by 10X Cell Ranger is unclear. To make Repli and FPKM somewhat
          # more comparable (although there may still be genomic features with low
          # read count in the Repli window and we will remove those features
          # before any Repli computation), we will filter Repli features using 10X
          # FPKM feature criteria.
          subset(names(.) %in% names(quartile.factor_Somatic)),
        dirname(
          c(
            chic.bw_H3K4_Somatic,
            chic.bw_H3K27_Somatic,
            chic.bw_H3K9_Somatic
          )[1]
        ),
        assay.data.sc,
        'tj',
        'Repli Quartile',
        setNames(repli_quartile_fills, NULL)
      ) %>%
        list %>%
        tibble(
          name = "Somatic_marks",
          figure = .,
          width = 9,
          height = 4
        ) %>%
        save_figures("repli/profile", extension, ., dpi = 300),
      format = "file"
    ),
    tar_target(
      chic.results_Germline,
      save_figures(
        "figure/Germline",
        extension,
        tribble(
          ~name, ~figure, ~width, ~height,
          "CHIC-TSS-AllMarks-RNAseq-Quartile",
          fpkm.chic.quarter.tss.plot_Germline,
          9, 4,
          "CHIC-AllMarks-RNAseq-Quartile",
          fpkm.chic.quarter.plot_Germline,
          15, 4,
          "CHIC-AllMarks-Peak-Annotation",
          plot.chic.peak.location_Germline
          + scale_y_continuous(
            limits = limits_plot.chic.peak.location,
            expand = expansion(mult = c(0, 0.05))
          ),
          6, 3,
          "CHIC-H3-Periodicity",
          chic.plot.psd_Germline,
          4, 2
        ),
        dpi = 300
      ),
      format = 'file'
    ),
    tar_target(
      chic.results_Somatic,
      save_figures(
        "figure/Somatic",
        extension,
        tribble(
          ~name, ~figure, ~width, ~height,
          "CHIC-TSS-AllMarks-RNAseq-Quartile",
          fpkm.chic.quarter.tss.plot_Somatic,
          9, 4,
          "CHIC-AllMarks-RNAseq-Quartile",
          fpkm.chic.quarter.plot_Somatic,
          15, 4,
          "CHIC-AllMarks-Peak-Annotation",
          plot.chic.peak.location_Somatic
          + scale_y_continuous(
            limits = limits_plot.chic.peak.location,
            expand = expansion(mult = c(0, 0.05))
          ),
          6, 3,
          "CHIC-H3-Periodicity",
          chic.plot.psd_Somatic,
          4, 2
        ),
        dpi = 300
      ),
      format = 'file'
    )
  ),

  # CHIC Power Spectral Density
  apply(
    chic.fpkm.data %>% dplyr::rename(driver. = "driver"),
    1,
    \(v) with(
      as.list(v),
      list(
        tar_target_raw(
          paste0("chic.raw.list_input_", name),
          chic.samples %>%
            dplyr::filter(molecule == "H3", driver == driver.) %>%
            pull(sample) %>%
            paste0("chic.raw_", .) %>%
            rlang::syms() %>%
            append(list("list"), .) %>%
            do.call(call, ., quote=T) %>%
            list(
              "setNames",
              .,
              chic.samples %>%
                dplyr::filter(molecule == "H3", driver == driver.) %>%
                pull(sample) %>%
                append(list("c"), .) %>%
                do.call(call, ., quote=T)
            ) %>%
            do.call(call, ., quote=T)
        ),
        tar_target_raw(
          paste0("chic.psd.gene.list_", name),
          if (name == "Germline")
            quote(sct_gene_lists$germline)
          else quote(sct_gene_lists$somatic)
        ),
        tar_target_raw(
          paste0("chic.psd.obs_", name),
          list("cbind_features_at_tss") %>% append(
            rlang::syms(
              list(
                paste0("chic.raw.list_input_", name),
                paste0("chic.psd.gene.list_", name),
                "assay.data.sc"
              )
            )
          ) %>%
            do.call(call, ., quote=T)
        ),
        tar_target_raw(
          paste0("chic.psd_", name),
          call(
            "psd_centered_features",
            rlang::sym(paste0("chic.psd.obs_", name))
          ),
          packages = "gsignal"
        ),
        tar_target_raw(
          paste0("chic.plot.psd_", name),
          call(
            "heatmap_psd_centered_features",
            rlang::sym(paste0("chic.psd_", name))
          )
        )
      )
    )
  ),

  #,

  # tar_target(repli_dir, 'Repli-Seq/Validation', format='file'),
  # tar_target(name = repli.raw, command = load_repli(repli_dir), cue=tar_cue('never')),
  # tar_target(name = repli, command = normalize_repli(repli.raw), cue=tar_cue('never')),
  # tar_target(name = repli.hmm, command = repli_fit_hmm(repli), cue=tar_cue('never'))

  # Repli-Seq
  tar_file(run_fastqc_sh, "scripts/run_fastqc.sh"),
  repli_targets,

  targets.chic,
  targets.flybase,
  targets.quantification,
  targets.sce
)
 
