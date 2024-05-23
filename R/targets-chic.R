library(dplyr)
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
              str_replace(bowtie[1], "\\..*", ""),
              paste0(batch, "/", sample, "_R1_001.fastq.gz"),
              paste0(batch, "/", sample, "_R2_001.fastq.gz"),
              output_path
            )
          )
          output_path
        }
      ),
      cue = tar_cue("never"),
      packages = c(
        "dplyr",
        "processx",
        "stringr"
      )
    ),
    tar_target(
      chic.readcount,
      as.integer(
        c(
          run(
            "bash",
            c(
              "-c",
              paste0(
                "rclone cat sharepoint:'Bio Data'/",
                batch,
                "/",
                sample,
                "_R1_001.fastq.gz",
                " | gunzip -c | wc -l"
              )
            )
          )$stdout,
          run(
            "bash",
            c(
              "-c",
              paste0(
                "rclone cat sharepoint:'Bio Data'/",
                batch,
                "/",
                sample,
                "_R2_001.fastq.gz",
                " | gunzip -c | wc -l"
              )
            )
          )$stdout
        )
      ),
      packages = "processx"
    )
  )
)

targets.chic <- list(
  targets.chic.aligned,

  # ChIC loaded BAM files to KDE and the bandwidth-equivalent rectangular windows.
  tar_map(
    cross_join(
      bowtie.refs,
      chic.samples
    ) %>%
      rowwise %>%
      mutate(
        suffix = str_glue("{sample}_{name}"),
        bulk_reads_misc = rlang::syms(str_glue("bulk_reads_misc_chic.bam_{sample}_{name}")),
        bulk_reads_idxstats = rlang::syms(str_glue("bulk_reads_idxstats_chic.bam_{sample}_{name}")),
        bulk_reads_split = mapply(
          \(ref_name, n) rlang::syms(
            str_glue(
              "bulk_reads_{names(if (name == \"masked\") masked.lengths else chr.lengths)}_{n}"
            )
          ) %>%
            setNames(names(if (name == "masked") masked.lengths else chr.lengths)) %>%
            append(list("list"), .) %>%
            do.call(call, ., quote=T),
          name,
          str_glue("chic.bam_{sample}_{name}"),
          SIMPLIFY=F
        )
      ),
    names = suffix,
    tar_target(
      chic.sparseVector,
      append(
        sapply(
          bulk_reads_split,
          \(df) bam_cover_paired_end_fragments_bp(
            df, min_mapq = 20, min_fl = 120, max_fl = 220, markdup = TRUE
          )[[1]],
          simplify=FALSE
        ),
        as.list(
          (
            bulk_reads_misc %>%
              paired_end_reads_to_fragment_lengths %>%
              # mutate(rname = rname %>% factor(setdiff(levels(.), names(if (name == "masked") masked.lengths else chr.lengths)))) %>%
              filter(between(mapq, 20, 254), between(length, 120, 220)) %>%
              split(.$rname) %>%
              sapply(\(df) nrow(df))
          )[setdiff(levels(bulk_reads_misc$rname), names(if (name == "masked") masked.lengths else chr.lengths))] %>%
            setNames(
              setdiff(levels(bulk_reads_misc$rname), names(if (name == "masked") masked.lengths else chr.lengths))
            )
        )
      )
    ),
    tar_target(
      chic.chr.fpkm,
      GRanges(
        seqnames = names(chic.sparseVector),
        IRanges(
          start = 1,
          width = pull(bulk_reads_idxstats, "rlength", "rname")[names(chic.sparseVector)]
        ),
        score = sapply(chic.sparseVector, sum) %>%
          `/`(sum(.)) %>%
          `*`(1000^3 / pull(bulk_reads_idxstats, "rlength", "rname")[names(chic.sparseVector)]),
        seqlengths = pull(bulk_reads_idxstats, "rlength", "rname")[names(chic.sparseVector)]
      ) %>%
        `metadata<-`(value = list(est_library_size = sum(sapply(chic.sparseVector, sum)))) %>%
        split(names(chic.sparseVector) %>% factor(., .))
    ),
    tar_target(
      chic.granges.fpkm,
      chic.chr.fpkm %>%
        replace(
          which(sapply(chic.sparseVector, length) > 1),
          chic.sparseVector %>%
            subset(sapply(chic.sparseVector, length) > 1) %>%
            mapply(
              \(rname, v) density_est_granges_approx(
                chic.chr.fpkm[[rname]],
                obs_vector = v,
                bw = 25,
                sample_rate = 10
              ),
              names(.),
              .
            ) %>%
            GRangesList
        )
    )
  ),

  tar_map(
    chic.experiments %>%
      cross_join(dplyr::rename(bowtie.refs, reference="name")) %>%
      mutate(driver = driver %>% replace(. == "Nos", "nos")) %>%
      rowwise %>%
      mutate(
        chic_sample_names = subset(
          chic.samples$sample,
          (chic.samples$driver %>% replace(. == "Nos", "nos")) == driver &
          (chic.samples$molecule == "H3" | chic.samples$group == mark)
        ) %>%
          list,
        chic.granges.fpkm = rlang::syms(str_glue("chic.granges.fpkm_{chic_sample_names}_{reference}")) %>%
          list,
        experimental_group = (
          chic.samples$group[match(chic_sample_names, chic.samples$sample)]
        ) %>%
          list,
        molecule = (
          chic.samples$molecule[match(chic_sample_names, chic.samples$sample)]
        ) %>%
          list,
        rep = (
          chic.samples$rep[match(chic_sample_names, chic.samples$sample)]
        ) %>%
          list
      ),
    names = mark | name | reference,
    tar_target(
      chic.experiment.sample.names, chic_sample_names
    ),
    tar_target(
      chic.experiment.granges,
      unlist(chic.granges.fpkm[[1]]) %>%
        attributes %>%
        with(
          GRanges(
            seqnames,
            ranges,
            seqlengths = seqlengths(chic.granges.fpkm[[1]]),
            score = sapply(
              chic.granges.fpkm,
              \(lst) unlist(lst)$score
            )
          ) %>%
            `metadata<-`(
              value = list(
                est_library_size = sapply(
                  chic.granges.fpkm,
                  \(lst) lst[[1]]@metadata$est_library_size
                )
              )
            )
        )
    ),
    tar_target(
      chic.experiment.model.data,
      tibble(
        molecule = factor(molecule),
        experimental_group = factor(experimental_group) %>%
          fct_relevel(mark),
        rep = factor(as.character(rep)) %>%
          `contrasts<-`(value = contr.helmert(length(levels(.)))),
        # Let's use one window in the exon of "tj" as the example locus for fitting an lm.
        score = chic.experiment.granges@elementMetadata[
          which.max(chic.experiment.granges@seqnames == "2L" & chic.experiment.granges@ranges@start == 19465001),
        ] %>%
          unlist,
        est_library_size = chic.experiment.granges@metadata$est_library_size
      ) %>%
        as.list %>%
        with(
          append(
            .,
            {
              experimental_group_rep <- interaction(experimental_group, rep)
              contrasts(experimental_group_rep) <- contr.helmert(length(levels(experimental_group_rep)))
              list(experimental_group_rep = experimental_group_rep)
            }
          )
        )
    ),
    tar_target(
      chic.experiment.lm.sample,
      with(
        chic.experiment.model.data,
        lm(
          "score ~ 0 + molecule + experimental_group_rep",
          list(score, molecule, experimental_group_rep),
          weights = 1 / sqrt(chic.experiment.model.data$est_library_size),
          subset = chic.experiment.model.data$experimental_group == marks
        )
      )
    )
  ),

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
  ),

  # Profiles of ChIC faceted by quantification
  tar_map(
    experiment.driver %>% mutate(quartile.factor = rlang::syms(str_glue("quartile.factor_{celltype}"))),
    names = celltype,
    tar_target(
      sc_chr_factor,
      tibble(
        rowname = names(quartile.factor),
        chr = read.csv(assay.data.sc)$chr %>%
          fct_recode(`2`="2L", `2`="2R", `3`="3L", `3`="3R") %>%
          factor(c("2", "3", "4", "X", "Y")),
        quant = quartile.factor %>%
          fct_recode(off="Q1", low="Q2", med="Q3", high="Q4"),
        fct = interaction(chr, quant)
      )
    ),
    tar_target(
      tss_sc_chr_quartile_data,
      chic_average_profiles(
        pull(sc_chr_factor, fct, rowname),
        dirname(
          c(
            chic.bw_H3K4_Germline,
            chic.bw_H3K27_Germline,
            chic.bw_H3K9_Germline
          )[1]
        ),
        assay.data.sc,
        driver,
        # "CPM Quartile",
        "",
        # fake colors
        sprintf("#%06d", 10 * seq_along(levels(sc_chr_factor$fct)))
      )$data,
      format = "parquet"
    ),
    tar_target(
      tss_sc_chr_data,
      chic_average_profiles(
        pull(sc_chr_factor, chr, rowname),
        dirname(
          c(
            chic.bw_H3K4_Germline,
            chic.bw_H3K27_Germline,
            chic.bw_H3K9_Germline
          )[1]
        ),
        assay.data.sc,
        driver,
        # "CPM Quartile",
        "",
        # fake colors
        sprintf("#%06d", 10 * seq_along(levels(sc_chr_factor$chr)))
      )$data,
      format = "parquet"
    )
  )
)