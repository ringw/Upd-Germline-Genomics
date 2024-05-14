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

target_chic_load_raw_expr <- function(group, driver) {
  group. <- group
  driver. <- driver
  samples <- chic.samples %>%
    filter(group == group., driver == driver.) %>%
    reframe(molecule, rep, sample, pileup = rlang::syms(str_glue("chic.raw_{sample}")))
  rlang::call2("tibble", rowname = samples$sample, molecule = samples$molecule, rep = call("factor", samples$rep), pileup = samples$pileup)
}

targets.chic.aligned <- tar_map(
  bowtie.refs,
  names = name,
  tar_map(
    chic.samples,
    names = sample,
    unlist = FALSE,
    tar_file(
      chic.bam,
      with(
        list(name=name, group=group, sample=sample) %>%
          with(
            list(output_path = str_glue("chic/{name}/{group}/{sample}.bam"), batch=batch, sample=sample)
          ),
        {
          run(
            "bash",
            c(
              "-i",
              align_chic_lightfiltering,
              flybase.bowtie %>% paste("chic_bowtie2", sep="/"),
              paste0(batch, "/", sample, "_R1_001.fastq.gz"),
              paste0(batch, "/", sample, "_R2_001.fastq.gz"),
              output_path
            )
          )
          output_path
        }
      ),
      cue = tar_cue("never")
    )
  )
)

targets.chic <- list(
  targets.chic.aligned,
  apply(
    chic.experiments %>%
      mutate(driver = driver %>% replace(. == "Nos", "nos")),
    1,
    \(v) tar_target_raw(
      str_glue("chic.regression.result_{v['mark']}_{v['name']}"),
      call(
        "chic_perform_regression",
        target_chic_load_raw_expr(group = v["mark"], driver = v["driver"])
      )
    )
  ),
  tar_map(
    chic.experiments %>%
      rowwise %>%
      mutate(
        driver = driver %>% replace(. == "Nos", "nos"),
        chic.regression.result = rlang::syms(str_glue("chic.regression.result_{mark}_{name}")),
        chic_samples_df = target_chic_load_raw_expr(group = mark, driver = driver) %>%
          list
      ),
    names = mark | name,
    tar_file(
      chic.bw.2,
      tibble(
        data = list(
          chic_samples_df %>%
            rowwise %>%
            mutate(
              bw = if (molecule == "H3") 250 else 25,
              density = smooth_sparse_vector_to_density(pileup, feature.rle, bw=bw) %>% list
            )
        ),
        mymanova = chic_perform_regression(data[[1]]) %>% list,
        myfoldchange = chic_regression_fold_change(mymanova[[1]], data[[1]]$density[[1]]) %>% list,
        mywindow = feature_lengths_to_sliding_windows(feature.lengths, flybase.lengths, 10L, 10L) %>% list,
        mytrack = track_kde_to_sliding(myfoldchange[[1]], mywindow[[1]]) %>% list,
        myt2 = coverage_interp_obs_granges(
          mytrack[[1]],
          replace_na(mytrack[[1]]$coverage >= 1, FALSE)
          &
          is.finite(mytrack[[1]]$fold_change)
        ) %>%
          list,
        output_path = paste0("chic/", driver, "_", mark, ".new.FE.bw"),
        do_export = rtracklayer::export(
          {
            seqlengths(myt2[[1]]) <- flybase.lengths
            myt2[[1]]$fold_change[is.na(myt2[[1]]$fold_change) & !as.logical(myt2[[1]]@seqnames %in% names(chr.lengths))] <- 0
            myt2[[1]]$fold_change <- pmax(-0.1, myt2[[1]]$fold_change)
            GRanges(myt2[[1]], score = myt2[[1]]$fold_change)
          },
          output_path,
          "bigwig"
        ) %>%
          list
      ) %>%
        pull(output_path),
      packages = tar_option_get("packages") %>% c("tidyr")
    )
  )
)