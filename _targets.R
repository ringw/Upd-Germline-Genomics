# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)

# Set target options:
tar_option_set(
  packages = c(
    'AnnotationDbi',
    'dplyr',
    'drosophila2.db',
    'forcats',
    'glmGamPoi',
    'irlba',
    'Matrix',
    'purrr',
    'reticulate',
    'scran',
    'scuttle',
    'Seurat',
    'SingleCellExperiment',
    'stringr',
    'tibble',
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

source('R/repli-loading.R')
source('R/repli-model.R')

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

list(
  tar_target(
    sce.mt.features,
    read_mt_genome(paste0(sce.data$tenx_path[1], '/features.tsv.gz'))
  ),
  tar_target(
    scRNAseq,
    'scRNA-seq',
    format='file'
  ),
  tar_target(
    sce.present.features,
    select_features(str_replace(sce.data$tenx_path, 'scRNA-seq', scRNAseq))
  ),
  tar_map(
    values = sce.data,
    tar_target(tenx_file, tenx_path, format='file'),
    tar_target(sce, load_flybase(tenx_file, batch, sce.present.features, sce.mt.features))
  ),
  tar_target(nos.1, call_nos.1(sce_nos.1_scRNA.seq.Nos.Upd_H3.GFP_Rep1.outs.filtered_feature_bc_matrix)),
  tar_target(nos.2, call_nos.2(sce_nos.2_scRNA.seq.Nos.Upd_H3.GFP_Rep2.outs.filtered_feature_bc_matrix)),
  tar_target(tj.1, call_tj.1(sce_tj.1_scRNA.seq.tj.Upd_H3.GFP_Rep1.outs.filtered_feature_bc_matrix)),
  tar_target(tj.2, call_tj.2(sce_tj.2_scRNA.seq.tj.Upd_H3.GFP_Rep2.outs.filtered_feature_bc_matrix)),
  tar_target(Upd_sc, filter_integrate_data(list(nos.1,nos.2,tj.1,tj.2))),
  tar_target(metadata, analyze_sce_to_csv(list(nos.1=nos.1, nos.2=nos.2, tj.1=tj.1, tj.2=tj.2), 'scRNA-seq-Metadata.csv'), format='file'),
  tar_target(Upd_glm, fit_glm(list(nos.1=nos.1, nos.2=nos.2, tj.1=tj.1, tj.2=tj.2), metadata)),

  tar_target(repli_dir, 'Repli-Seq/Validation', format='file'),
  tar_target(name = repli.raw, command = load_repli(repli_dir)),
  tar_target(name = repli, command = normalize_repli(repli.raw)),
  tar_target(name = repli.hmm, command = repli_fit_hmm(repli))
)
