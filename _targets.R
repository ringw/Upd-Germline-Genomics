# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# HDF5Array is used for glm and will write to custom location.
# HDF5Array::setHDF5DumpDir('D:/HDF5Array_dump')

# Load packages required to define the pipeline:
library(magrittr)
library(targets)
library(tarchetypes)
library(tibble)

# Set target options:
tar_option_set(
  packages = c(
    'AnnotationDbi',
    'apeglm',
    'BiocIO',
    'cowplot',
    'DESeq2',
    'dplyr',
    'drosophila2.db',
    'forcats',
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
    'purrr',
    'readxl',
    'reshape2',
    'reticulate',
    'rtracklayer',
    'scales',
    'scran',
    'scuttle',
    'Seurat',
    'SingleCellExperiment',
    'stringr',
    'tibble',
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

# chr lengths from dmel r6.47
chr.lengths = c(
  `2L`=23513712, `2R`=25286936, `3L`=28110227, `3R`=32079331, `4`=1348131, X=23542271, Y=3667352
)

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
  'germline', c(1,0,0,0,0,0,0),
  'somatic', c(1,1,0,0,0,0,0),
  'spermatocyte', c(1,0,1,0,0,0,0),
  'muscle', c(1,0,0,1,0,0,0)
)

repli.processing = data.frame(
  name=c('.q20.bedgraph', '.q20.markdup.bedgraph', '.q20.p10.markdup.bedgraph', '.q20.near10.bedgraph', '.q20.near25.bedgraph'),
  display=c('OrigPileup', 'DupePosition', 'DownsampleAndDupePosition', 'DeleteNearbyPositions', 'DeleteNearbyPositions25')
)

repli.coverage.contrasts = list(
  EarlyLate = c(1,0,0,-1),
  Weighted = c(7/8, 5/8, 3/8, 1/8)
)

repli.coverage <- list(
  ident = c('Tj','Tj','GSC','GSC'),
  contrast = c('EarlyLate', 'Weighted','EarlyLate','Weighted')
)

chic.fpkm.data = tribble(
  ~name, ~contrast, ~driver,
  'Germline', c(1,0,0,0,0,0,0), 'Nos',
  'Somatic', c(1,1,0,0,0,0,0), 'tj'
)

chic.mark.data = tribble(~mark, 'H3K4', 'H3K27', 'H3K9')

sce_targets <- tar_map(
  unlist = FALSE,
  sce.data,
  tar_target(tenx_file, tenx_path, format='file'),
  tar_target(sce, load_flybase(tenx_file, batch, sce.present.features, sce.mt.features, metafeatures)),

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
      sum(
        sce$nCount_RNA %>%
          subset(
            read.csv(
              metadata, row.names = 1
            )[
              paste(batch, names(.), sep='_'),
              'ident'
            ] == cluster
          )
      )
    )
  )
)

list(
  tar_target(
    scRNAseq,
    'scRNA-seq',
    format='file'
  ),
  tar_target(
    sce.mt.features,
    read_mt_genome(paste0(sce.data$tenx_path[1], '/features.tsv.gz'))
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
    'dmel-all-r6.47.gtf.gz'
  ),
  tar_target(
    metafeatures,
    create_meta_features(
      tenx_file_nos.1_scRNA.seq.Nos.Upd_H3.GFP_Rep1.outs.filtered_feature_bc_matrix,
      flybase.annotations, flybase.gtf, sce.present.features, sce.mt.features,
      'scRNA-seq-Meta-Features.csv'), format='file'),
  tar_target(nos.1, call_nos.1(sce_nos.1_scRNA.seq.Nos.Upd_H3.GFP_Rep1.outs.filtered_feature_bc_matrix)),
  tar_target(nos.2, call_nos.2(sce_nos.2_scRNA.seq.Nos.Upd_H3.GFP_Rep2.outs.filtered_feature_bc_matrix)),
  tar_target(tj.1, call_tj.1(sce_tj.1_scRNA.seq.tj.Upd_H3.GFP_Rep1.outs.filtered_feature_bc_matrix)),
  tar_target(tj.2, call_tj.2(sce_tj.2_scRNA.seq.tj.Upd_H3.GFP_Rep2.outs.filtered_feature_bc_matrix)),
  tar_target(Upd_sc, filter_integrate_data(list(nos.1,nos.2,tj.1,tj.2)), cue=tar_cue('never')),
  tar_target(scRNA_seq_figures, Upd_sc_figures('scRNA-seq-Figure', Upd_sc), format='file'),
  tar_target(metadata, analyze_sce_to_csv(list(nos.1=nos.1, nos.2=nos.2, tj.1=tj.1, tj.2=tj.2), 'scRNA-seq-Metadata.csv'), format='file'),
  tar_target(
    Upd_glm,
    fit_glm(
      c(nos.1=tenx_file_nos.1_scRNA.seq.Nos.Upd_H3.GFP_Rep1.outs.filtered_feature_bc_matrix,
        nos.2=tenx_file_nos.2_scRNA.seq.Nos.Upd_H3.GFP_Rep2.outs.filtered_feature_bc_matrix,
        tj.1=tenx_file_tj.1_scRNA.seq.tj.Upd_H3.GFP_Rep1.outs.filtered_feature_bc_matrix,
        tj.2=tenx_file_tj.2_scRNA.seq.tj.Upd_H3.GFP_Rep2.outs.filtered_feature_bc_matrix),
      metadata),
    # glm is quite slow and doesn't need to re-run as we are adding other colData to the seurat
    cue=tar_cue('never')),
  tar_target(
    Upd_regression_somatic,
    apeglm_coef_table(Upd_glm, sce.present.features)
  ),
  tar_target(
    Upd_regression_spmtc,
    apeglm_coef_table_contrast(Upd_glm, sce.present.features, apeglm_contrasts$spermatocyte)
  ),
  tar_target(
    Upd_regression_muscle,
    apeglm_coef_table_contrast(Upd_glm, sce.present.features, apeglm_contrasts$muscle)
  ),
  tar_target(scripts_dir, "scripts", format = "file"),
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
  tar_target(
    Upd_cpm,
    pseudobulk_cpm(
      Upd_pseudobulk,
      Upd_pseudobulk_sf,
      flybase.gtf
    )
  ),
  tar_target(
    Upd_fpkm_bam,
    tx_cpm_to_fpkm(Upd_cpm, metafeatures)
  ),
  tar_target(
    Upd_fpkm_longest_isoform,
    longest_isoform_fpkm(Upd_glm, metafeatures)
  ),
  # tar_target(
  #   Upd_fpkm,
  #   Upd_fpkm_bam %>%
  #     replace(!is.finite(.), Upd_fpkm_longest_isoform[!is.finite(.)])
  # ),
  tar_target(Upd_fpkm, Upd_fpkm_bam),
  tar_target(
    sc_excel,
    publish_excel_results(Upd_glm, Upd_regression_somatic, Upd_regression_spmtc, Upd_regression_muscle, metafeatures, 'scRNA-seq-Regression/Enriched-Genes.xlsx'),
    format='file'
  ),
  # tar_target(
  #   germline_fpkm,
  #   create_fpkm_reference_longest_isoform(
  #     metafeatures, Upd_glm, c(1,0,0,0,0,0,0), 'scRNA-seq-Regression/Germline-FPKM.bed'
  #   ),
  #   format='file'
  # ),
  # tar_target(
  #   somatic_fpkm,
  #   create_fpkm_reference_longest_isoform(
  #     metafeatures, Upd_glm, c(1,1,0,0,0,0,0), 'scRNA-seq-Regression/Somatic-FPKM.bed'
  #   ),
  #   format='file'
  # ),
  tar_map(
    chic.fpkm.data,
    names = name,
    tar_target(
      bed,
      reference_sort_by_fpkm_table(
        Upd_fpkm_bam, tolower(name), metafeatures,
        paste0("scRNA-seq-Regression/", name, "-FPKM.bed")
      ),
      format = "file"
    ),
    tar_map(
      chic.mark.data,
      names = mark,
      tar_target(
        tss_mark_matrix,
        flybase_big_matrix(
          load_flybase_bed(bed),
          paste0("chic/", driver, "_", mark, ".q5.bw")
        )
      ),
      tar_target(
        tss_mark_heatmap,
        display_tss_tile_matrix(
          tss_mark_matrix,
          paste0("figure/", name, "/FPKM-", mark, ".pdf")
        ),
        format = "file"
      ),
      tar_target(
        tss_mark_heatmap_png,
        display_tss_tile_matrix(
          tss_mark_matrix,
          paste0("figure/", name, "/FPKM-", mark, ".png")
        ),
        format = "file"
      )
    )
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
  tar_target(repli.concat.markdup, pseudobulk_repli(load_concat('repli/bedgraph.in'))),
  tar_target(
    repli_qc,
    repli_qc_make_figures(
      list(
        DupePositionReplicated=repli.concat.markdup,
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
  )

  #,

  # tar_target(repli_dir, 'Repli-Seq/Validation', format='file'),
  # tar_target(name = repli.raw, command = load_repli(repli_dir), cue=tar_cue('never')),
  # tar_target(name = repli, command = normalize_repli(repli.raw), cue=tar_cue('never')),
  # tar_target(name = repli.hmm, command = repli_fit_hmm(repli), cue=tar_cue('never'))
)
 