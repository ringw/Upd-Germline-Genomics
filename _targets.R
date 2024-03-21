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
library(targets)
library(tarchetypes)
library(tibble)
library(viridis)

pdfFonts(sans = pdfFonts()$Helvetica)
postscriptFonts(sans = postscriptFonts()$Helvetica)

# Set target options:
tar_option_set(
  packages = c(
    'AnnotationDbi',
    'apeglm',
    'BiocIO',
    'Cairo',
    'cowplot',
    'DESeq2',
    'dplyr',
    'drosophila2.db',
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

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

sce.data = data.frame(
  batch=c('nos.1','nos.2','tj.1','tj.2'),
  tenx_path=paste0(
    'scRNA-seq/',
    c('Nos','Nos','tj','tj'),
    '-Upd_H3-GFP_Rep',
    c(1,2,1,2),
    '/outs/filtered_feature_bc_matrix'
  )
)

sce.clusters = data.frame(cluster = c('germline','somatic','spermatocyte','muscle'))
sce.clusters <- tribble(
  ~ cluster, ~ contrast,
  'germline', c(1, CySCoverGSC=-0.5, 0,0,0,0,0),
  'somatic', c(1, CySCoverGSC=0.5, 0,0,0,0,0),
  'spermatocyte', c(1,0,1,0,0,0,0),
  'muscle', c(1,0,0,1,0,0,0)
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

chic.samples = read.csv('chic/chic_samples.csv') %>%
  subset(sample != "") %>%
  # Further subsetting due to a sample with small library count
  subset(sample != "GC3768013_S13_L001")

chic.fpkm.data <- tribble(
  ~name, ~contrast, ~driver,
  'Germline', c(1,0,0,0,0,0,0), 'Nos',
  'Somatic', c(1,1,0,0,0,0,0), 'tj'
)
# We are going to pull a "bed" target in the cross product below, which defines
# the genomic features shown in a heatmap for this particular experiment.
chic.fpkm.data$bed_sym <- rlang::syms(paste0("bed_", chic.fpkm.data$name))

chic.mark.data = tribble(~mark, 'H3K4', 'H3K27', 'H3K9')

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
# ChIC experiments: For every driver x mark.
chic.experiments <- chic.fpkm.data %>% cross_join(chic.mark.data)
# Pull the "chic" (track) target by name for the driver x mark targets.
chic.experiments <- chic.experiments %>% within({
  experiment_name <- paste(mark, name, sep="_")
  chic_input_sym <- rlang::syms(paste("chic.smooth_250", "input", mark, name, sep="_"))
  chic_mod_sym <- rlang::syms(paste("chic.smooth_25", "mod", mark, name, sep="_"))
  chic_smooth_input_sym <- rlang::syms(paste("chic.smooth", "input", mark, name, sep="_"))
  chic_smooth_250_input_sym <- rlang::syms(paste("chic.smooth_250", "input", mark, name, sep="_"))
  chic_smooth_mod_sym <- rlang::syms(paste("chic.smooth", "mod", mark, name, sep="_"))
  chic_smooth_125_mod_sym <- rlang::syms(paste("chic.smooth_125", "mod", mark, name, sep="_"))
  chic_smooth_250_mod_sym <- rlang::syms(paste("chic.smooth_250", "mod", mark, name, sep="_"))
  chic_input_bam <- rlang::syms(paste("chic.merge.bam", "input", mark, name, sep="_"))
  chic_mod_bam <- rlang::syms(paste("chic.merge.bam", "mod", mark, name, sep="_"))
})
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
  sce.data,
  names = batch,
  tar_target(
    tenx_file, tenx_path
    # Don't checksum the 10X outputs. Unfortunately we will want to avoid using
    # "file" with the tar cue "never", which is still reading in all of the
    # files (presumably to do a checksum) every time we call tar_make!
    # format='file'
  ),
  tar_target(
    seurat_qc,
    load_flybase(
      tenx_file, batch, sce.present.features, sce.mt.features, metafeatures
    )
  ),
  tar_target(
    sctransform_quantile,
    load_flybase(
      tenx_file, batch, sce.features, sce.mt.features, metafeatures,
      return.only.var.genes = FALSE, run_pca = FALSE
    ) %>%
      quantify_quantiles(metadata)
  ),
  tar_target(
    lognormalize_quantile,
    load_flybase(
      tenx_file, batch, sce.present.features, sce.mt.features, metafeatures,
      return.only.var.genes = FALSE, run_pca = FALSE
    ) %>%
      quantify_quantiles(
        metadata, assay = "RNA", slot = "data", per_million = FALSE
      )
  ),

  # Map over clusters (below, from metadata) and extract pseudobulk counts.
  tar_map(
    unlist = FALSE,
    sce.clusters,
    names = cluster,
    tar_target(
      pseudobulk.tx,
      download_bam_transcripts(
        paste0(scripts_dir, '/download_bam_tx.sh'),
        metadata,
        batch,
        cluster,
        paste0('scRNA-seq-Regression/pseudobulk_', batch, '_', cluster, '.txt')
      ),
      format = 'file'
    ),
    tar_target(
      pseudobulk.library.size,
        seurat_qc$nCount_RNA %>%
        subset(
          read.csv(metadata, row.names = 1)[names(.), 'ident'] == cluster
        ) %>% sum
    )
  )
)

list(
  tar_target(
    scRNAseq,
    'scRNA-seq',
    format='file',
    # Assume that 10X files do not change.
    cue = tar_cue("never")
  ),
  tar_target(
    sce.mt.features,
    read_mt_genome(paste0(sce.data$tenx_path[1], '/features.tsv.gz'))
  ),
  tar_target(
    sce.features,
    load_feature_names(str_replace(sce.data$tenx_path, 'scRNA-seq', scRNAseq))
  ),
  tar_target(
    sce.present.features,
    select_features(str_replace(sce.data$tenx_path, 'scRNA-seq', scRNAseq))
  ),
  tar_file(
    flybase.annotations,
    'fbgn_annotation_ID_fb_2022_04.tsv.gz'
  ),
  tar_file(
    flybase.gtf,
    # https://ftp.flybase.org/releases/FB2022_04/dmel_r6.47/gtf/dmel-all-r6.47.gtf.gz
    'dmel-all-r6.47.gtf.gz'
  ),
  tar_file(
    flybase.fa,
    'references/dmel-r6.47.fa.gz'
  ),
  tar_file(
    flybase.transposon,
    # https://ftp.flybase.org/releases/FB2022_04/precomputed_files/transposons/transposon_sequence_set.fa.gz
    'references/transposon_sequence_set.fa.gz'
  ),
  tar_target(
    transposon.lengths,
    flybase.transposon %>% read.table(quote="", sep="\x1B") %>% make_chr_lengths
  ),
  tar_file(
    flybase.genome,
    tibble(input_file=c(flybase.fa, flybase.transposon)) %>%
      rowwise() %>%
      mutate(contents = list(read.table(input_file, quote="", sep="\x1B"))) %>%
      ungroup() %>%
      summarize(
        contents = list(do.call(rbind, contents)),
        output_file = "references/chic_reference.fa",
        write_table = write.table(contents[[1]], output_file, row.names=F, col.names=F, quote=F)
      ) %>%
      pull(output_file),
    cue = tar_cue("never")
  ),
  tar_file(
    flybase.genome.index,
    tibble(
      input_file = flybase.genome,
      index_table = processx::run(
        "samtools", c("faidx", input_file)
      ) %>% list,
      faidx = paste0(input_file, ".fai")
    ) %>%
      pull(faidx)
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
  tar_file(
    flybase.bowtie,
    tibble(input_file=list(list(flybase.fa, flybase.transposon)), output_file="references/chic_bowtie2", output_base="chic_bowtie2") %>%
      mutate(
        rename_folder =
          if (file.exists(output_file))
            file.rename(output_file, paste0(output_file, "~")),
        make_folder = dir.create(output_file),
        run_bowtie2 = processx::run(
          "bowtie2-build",
          c(
            input_file[[1]] %>% append(list(sep=",")) %>% do.call(paste, .),
            paste(output_file, output_base, sep="/")
          )
        ) %>% list
      ) %>%
      pull(output_file)
  ),
  tar_target(
    metafeatures,
    create_meta_features(
      tenx_file_nos.1,
      flybase.annotations, flybase.gtf, sce.present.features, sce.mt.features,
      'scRNA-seq-Meta-Features.csv'), format='file'),
  tar_target(nos.1, call_nos.1(seurat_qc_nos.1)),
  tar_target(nos.2, call_nos.2(seurat_qc_nos.2)),
  tar_target(tj.1, call_tj.1(seurat_qc_tj.1)),
  tar_target(tj.2, call_tj.2(seurat_qc_tj.2)),
  tar_target(Upd_sc, filter_integrate_data(list(nos.1,nos.2,tj.1,tj.2))),
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
          "RNAseq-Quantification-Quarters",
          fpkm_quarter_density(Upd_fpkm),
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
          "AGO3","RpL22-like","RpL22","vas","nos","tj","lncRNA:roX1","lncRNA:roX2","Mst87F","soti","sunz","w-cup","Act57B","Dl","E(spl)m3-HLH"
        ) %>%
          mutate(Upd_sc = list(Upd_sc %>% NormalizeData)) %>%
          rowwise %>%
          mutate(
            gene_short_name = gene %>% str_replace("lncRNA:", ""),
            name = paste0("RNAseq-FeaturePlot-", gene_short_name),
            figure = list(Upd_sc %>% Upd_sc_feature_plot(gene) + labs(tag = gene_short_name)),
            width = 6,
            height = 4
          ) %>%
          subset(select=c(name, figure, width, height))
      )
    )
  ),
  tar_target(metadata, analyze_sce_to_csv(list(nos.1=nos.1, nos.2=nos.2, tj.1=tj.1, tj.2=tj.2), 'scRNA-seq-Metadata.csv'), format='file'),
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
    model_matrix,
    build_model_matrix(metadata)
  ),
  tar_target(
    Upd_glm,
    fit_glm(Upd_sc, model_matrix, metadata),
    # glm is quite slow and doesn't need to re-run as we are adding other colData to the seurat
    ), # cue=tar_cue('never')),
  tar_target(
    Upd_regression_somatic,
    apeglm_coef_table(Upd_glm)
  ),
  tar_target(
    Upd_regression_tid,
    apeglm_coef_table(Upd_glm, coef = 3)
  ),
  tar_target(
    Upd_regression_mscl,
    apeglm_coef_table(Upd_glm, coef = 4)
  ),
  # Split this tar_file, then remove tar_cue("never")
  tar_target(scripts_dir, "scripts", format = "file", cue = tar_cue("never")),
  sce_targets,
  tar_combine(
    Upd_pseudobulk,
    sce_targets[["pseudobulk.tx_germline"]],
    sce_targets[["pseudobulk.tx_somatic"]],
    sce_targets[["pseudobulk.tx_spermatocyte"]],
    sce_targets[["pseudobulk.tx_muscle"]],
    command = cbind(sce.data, tx_file = c(!!!.x))
  ),
  tar_combine(
    Upd_pseudobulk_sf,
    sce_targets[["pseudobulk.library.size_germline"]],
    sce_targets[["pseudobulk.library.size_somatic"]],
    sce_targets[["pseudobulk.library.size_spermatocyte"]],
    sce_targets[["pseudobulk.library.size_muscle"]],
    command = cbind(sce.data, size_factor = c(!!!.x))
  ),
  # Pseudobulk by the genotype (Nos-GAL4 or tj-GAL4).
  tar_target(
    supplemental_genes,
    "tj,vas" %>% strsplit(",") %>% unlist
  ),
  tar_map(
    tribble(
      ~name, ~batch_names, ~seurats,
      "nos", list("nos.1", "nos.2"), rlang::syms(c("nos.1", "nos.2")),
      "tj", list("tj.1", "tj.2"), rlang::syms(c("tj.1", "tj.2"))
    ),
    names = name,
    tar_target(
      supplemental_bulk_fpkm,
      gene_pseudobulk_fpkm(
        paste0(
          "scRNA-seq-Regression/pseudobulk_",
          name,
          c(".1", ".2"),
          "_filtered.txt"
        ),
        seurats,
        unlist(batch_names),
        flybase.gtf,
        metadata,
        metafeatures
      )
    )
  ),
  tar_target(
    supplemental_bulk_fpkm,
    gene_pseudobulk_fpkm(
      paste0(
        "scRNA-seq-Regression/pseudobulk_",
        c("nos", "nos", "tj", "tj"),
        c(".1", ".2", ".1", ".2"),
        "_filtered.txt"
      ),
      mapply(
        # Load all of the 10X matrices using sce.features (all features instead
        # of our min # cells features).
        \(n, f) load_flybase(
          f, n, sce.features, sce.mt.features, metafeatures, sctransform=FALSE
        ),
        c("nos.1", "nos.2", "tj.1", "tj.2"),
        c(tenx_file_nos.1, tenx_file_nos.2, tenx_file_tj.1, tenx_file_tj.2),
        SIMPLIFY = FALSE
      ),
      c("nos.1", "nos.2", "tj.1", "tj.2"),
      flybase.gtf,
      metadata,
      metafeatures
    )
  ),
  tar_target(
    supplemental_gene_list,
    c(
      "Act5C", "Act42A", "AGO3", "vas", "tj", "zfh1", "lncRNA:roX1",
      "lncRNA:roX2", "soti", "w-cup", "Act57B", "Mlp60A", "lncRNA:Hsromega"
    )
  ),
  tar_target(
    supplemental_bulk_figure,
    dot_plot_fpkm(
      list(Nos=supplemental_bulk_fpkm_nos, tj=supplemental_bulk_fpkm_tj),
      supplemental_gene_list
    )
  ),
  tar_map(
    sce.clusters,
    names = cluster,
    tar_target(
      supplemental_cluster_fpkm,
      gene_cluster_fpkm(
        Upd_fpkm,
        list(nos.1=nos.1, nos.2=nos.2, tj.1=tj.1, tj.2=tj.2),
        cluster,
        metadata
      )
    )
  ),
  tar_target(
    supplemental_cluster_dot_plot_figure,
    dot_plot_fpkm(
      list(
        GSC=supplemental_cluster_fpkm_germline,
        CySC=supplemental_cluster_fpkm_somatic,
        `other(germ)`=supplemental_cluster_fpkm_spermatocyte,
        mus=supplemental_cluster_fpkm_muscle
      ),
      supplemental_gene_list,
      oob_squish=TRUE
    )
  ),
  tar_target(
    supplemental_elbow_figure,
    data.frame(stdev = Upd_sc[['pca.subset']]@stdev, x = 1:50) %>%
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
    {
      filename <- "scRNA-seq-Regression/Cell-Cycle-Scoring.xlsx"
      write_excel_tables_list(
        write_Upd_sc_cell_cycle_phases(Upd_sc, cell_cycle_drosophila, metafeatures),
        filename
      )
      filename
    },
    format = "file"
  ),
  tar_target(
    supplemental_pca_figure,
    {
      filename <- "figure/Integrated-scRNAseq/Germline-Somatic-Pairs.pdf"
      dir.create(dirname(filename), recursive = TRUE, showW = FALSE)
      CairoPDF(filename, width = 16, height = 9)
      print(plot_Upd_pca_components(Upd_sc, load_cell_cycle_score_drosophila(cell_cycle_drosophila, metafeatures)))
      dev.off()
      filename
    },
    format = "file"
  ),
  tar_target(
    sct_gene_lists,
    determine_gene_list(sctransform_quantile, supplemental_bulk_fpkm)
  ),
  tar_target(
    sct_gene_venn,
    plot_gene_lists(
      sct_gene_lists
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
          "RNAseq-Genotype-Bulk-Dots", supplemental_bulk_figure, 4, 9,
          "RNAseq-Integrated-Cluster-Dots", supplemental_cluster_dot_plot_figure, 6, 9,
          "RNAseq-Validation-Elbow", supplemental_elbow_figure, 5, 2.5,
          "RNAseq-SCT-Gene-List-Venn-Area-Blank", sct_gene_venn, 4, 3
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
  # Estimate gene isoform means using summary of 10X BAM transcript id tags.
  tar_target(
    Upd_cpm_transcripts,
    pseudobulk_cpm(
      Upd_pseudobulk,
      Upd_pseudobulk_sf,
      flybase.gtf
    )
  ),
  # Call the most abundant gene isoform which maximizes FPKM value.
  tar_target(
    Upd_fpkm,
    tx_cpm_to_fpkm(Upd_cpm_transcripts, metafeatures)
  ),
  # Final table of gene abundances: CPM of the abundant isoform based on FPKM
  # calculation. Also join with H3-GFP gene mean coming from the GLM.
  tar_target(
    Upd_cpm,
    join_cpm_data(Upd_cpm_transcripts, Upd_fpkm, Upd_cpm_regression)
  ),
  tar_target(
    sc_excel,
    publish_excel_results(
      Upd_regression_somatic, Upd_regression_tid, Upd_regression_mscl,
      Upd_cpm_transcripts, Upd_cpm, Upd_fpkm, sctransform_quantile,
      supplemental_bulk_fpkm, metafeatures, flybase.gtf,
      'scRNA-seq-Regression/Enriched-Genes.xlsx'
    ),
    format='file'
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
  chic.raw.tracks,

  # ChIC coverage tracks.
  tar_target(
    chic.macs.bandwidth, list(500, 1000, 2000)
  ),
  # Background value for uniform FPKM. We have a multiplier of 1MM but this
  # amount is spread across the reference sequence length.
  tar_target(
    chic.macs.background.value,
    1000 * 1000 * 1000 / sum(c(chr.lengths, transposon.lengths))
  ),
  tar_target(
    chic.macs.background.track,
    c(chr.lengths, transposon.lengths) %>%
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
                Rle(
                  c(names(chr.lengths), names(transposon.lengths)),
                  c(chr.lengths, transposon.lengths)
                ),
                bw = bw,
                sample_size = ifelse(bw < 50, 10, 50)
              )
          ) %>%
          purrr::reduce(`+`) %>%
          round(digits = 1) %>%
          `/`(length(samples)) %>%
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
                  Rle(
                    c(names(chr.lengths), names(transposon.lengths)),
                    c(chr.lengths, transposon.lengths)
                  ),
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
              # TODO: put this expression into a target
              c(chr.lengths, transposon.lengths) %>%
                enframe("refseq", "len") %>%
                with(Rle(refseq, len)),
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
              # TODO: put this expression into a target
              c(chr.lengths, transposon.lengths) %>%
                enframe("refseq", "len") %>%
                with(Rle(refseq, len)),
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
              c(chr.lengths, transposon.lengths) %>%
                enframe("refseq", "len") %>%
                with(Rle(refseq, len)),
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
        Upd_fpkm, tolower(name), metafeatures,
        paste0("scRNA-seq-Regression/", name, "-FPKM.bed")
      ),
      format = "file"
    ),
    tar_target(
      quartile.factor,
      bed_flybase_quartile_factor(bed, metafeatures)
    )
  ),

  tar_target(
    chic.genome,
    c(chr.lengths, transposon.lengths) %>% data.frame(chr = names(.), len = .)
  ),
  tar_target(
    genomic_feature_factor,
    factor_genome(flybase.gtf, metafeatures)
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
        features = read.csv(metafeatures, row.names = 1)
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
      format = "file"
    ),
    tar_target(
      repli.quarters,
      analyze_repli_quarters(repli.hdf5, metafeatures, flybase.gtf)
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
      metafeatures,
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
      metafeatures,
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
      metafeatures,
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
      metafeatures,
      'tj',
      'FPKM Quartile',
      setNames(sc_quartile_annotations, NULL)
    )
  ),

  tar_target(
    name = fpkm.chic.sct.tss.plot_Germline,
    chic_average_gene_list_profiles(
      sct_gene_lists$germline,
      dirname(
        c(
          chic.bw_H3K4_Germline,
          chic.bw_H3K27_Germline,
          chic.bw_H3K9_Germline
        )[1]
      ),
      metafeatures,
      'Nos'
    )
  ),
  tar_target(
    name = fpkm.chic.sct.tss.plot_Somatic,
    chic_average_gene_list_profiles(
      sct_gene_lists$somatic,
      dirname(
        c(
          chic.bw_H3K4_Somatic,
          chic.bw_H3K27_Somatic,
          chic.bw_H3K9_Somatic
        )[1]
      ),
      metafeatures,
      'tj'
    )
  ), tar_target(
    name = fpkm.chic.sct.plot_Germline,
    chic_custom_gene_list_paneled_profile(
      sct_gene_lists$germline,
      dirname(
        c(
          chic.bw_H3K4_Germline,
          chic.bw_H3K27_Germline,
          chic.bw_H3K9_Germline
        )[1]
      ),
      metafeatures,
      'Nos'
    )
  ),
  tar_target(
    name = fpkm.chic.sct.plot_Somatic,
    chic_custom_gene_list_paneled_profile(
      sct_gene_lists$somatic,
      dirname(
        c(
          chic.bw_H3K4_Somatic,
          chic.bw_H3K27_Somatic,
          chic.bw_H3K9_Somatic
        )[1]
      ),
      metafeatures,
      'tj'
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
        metafeatures,
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
    ),
    tar_target(
      name = fpkm.chic.sct.tss_Germline,
      fpkm.chic.sct.tss.plot_Germline %>%
        list %>%
        tibble(
          name = "CHIC-TSS-AllMarks-RNAseq-SCT-AllGermlineGenes",
          figure = .,
          width = 8,
          height = 4
        ) %>%
        save_figures("figure/Germline", extension, ., dpi = 300),
      format = 'file'
    ),
    tar_target(
      name = fpkm.chic.sct_Germline,
      fpkm.chic.sct.plot_Germline %>%
        list %>%
        tibble(
          name = "CHIC-AllMarks-RNAseq-SCT-AllGermlineGenes",
          figure = .,
          width = 12,
          height = 4
        ) %>%
        save_figures("figure/Germline", extension, ., dpi = 300),
      format = 'file'
    ),
    tar_target(
      name = fpkm.chic.sct.tss_Somatic,
      fpkm.chic.sct.tss.plot_Somatic %>%
        list %>%
        tibble(
          name = "CHIC-TSS-AllMarks-RNAseq-SCT-AllSomaticGenes",
          figure = .,
          width = 8,
          height = 4
        ) %>%
        save_figures("figure/Somatic", extension, ., dpi = 300),
      format = 'file'
    ),
    tar_target(
      name = fpkm.chic.sct_Somatic,
      fpkm.chic.sct.plot_Somatic %>%
        list %>%
        tibble(
          name = "CHIC-AllMarks-RNAseq-SCT-AllSomaticGenes",
          figure = .,
          width = 12,
          height = 4
        ) %>%
        save_figures("figure/Somatic", extension, ., dpi = 300),
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
                "metafeatures"
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
  )

  #,

  # tar_target(repli_dir, 'Repli-Seq/Validation', format='file'),
  # tar_target(name = repli.raw, command = load_repli(repli_dir), cue=tar_cue('never')),
  # tar_target(name = repli, command = normalize_repli(repli.raw), cue=tar_cue('never')),
  # tar_target(name = repli.hmm, command = repli_fit_hmm(repli), cue=tar_cue('never'))
)
 
