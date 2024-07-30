repli.samples <- read.csv("repli/repli_samples.csv")
repli.samples$replication_value <- repli.samples$replication_value %>% factor(unique(.))
repli.samples$full <- repli.samples$replication_value
levels(repli.samples$full) <- c("Early", "Early-Mid", "Mid-Late", "Late")
repli.samples$abbrev <- repli.samples$replication_value
levels(repli.samples$abbrev) <- c("E", "G", "J", "L")
repli.samples <- repli.samples %>%
  mutate(name = paste0(genotype, "_", abbrev, rep))

repli.contrasts <- repli.samples %>%
  subset(select = c(full, replication_value)) %>%
  subset(!duplicated(.)) %>%
  within(replication_value <- as.numeric(as.character(replication_value))) %>%
  as_tibble()

repli.exp.contrast.numerator <- c(0.75, 0.25, -0.25, -0.75)
repli.exp.contrast.denominator <- c(1, 1, 1, 1)

chic_fseq_l2fc <- function(chic.bw.tracks) {
  lst <- as.list(
    grep("FSeq_Mark_L2FC", chic.bw.tracks, val = T)
  )
  lst[[1]]
}

targets.repli <- list(
  # FASTQ files: Align to BAM and count lines in FASTQ.
  tar_map(
    cross_join(
      bowtie.refs,
      dplyr::rename(repli.samples, repli_target = "name")
    ) %>%
      rowwise() %>%
      mutate(
        command_line = if (isTRUE(is_paired_end)) {
          list(
            call(
              "c",
              "-i",
              quote(align_chic_lightfiltering),
              substitute(str_replace(ref_paths[1], "\\..*", ""), list(ref_paths = bowtie)),
              str_glue("{filename}_R1_001.fastq.gz"),
              str_glue("{filename}_R2_001.fastq.gz"),
              quote(output_path)
            )
          )
        } else {
          list(
            call(
              "c",
              "-i",
              quote(align_repli_lightfiltering),
              substitute(str_replace(ref_paths[1], "\\..*", ""), list(ref_paths = bowtie)),
              str_glue("Upd_Tumor/Repli/{filename}"),
              quote(output_path)
            )
          )
        }
      ),
    names = repli_target | name,
    tar_file(
      repli.bam,
      with(
        list(name = name, repli_target = repli_target) %>%
          with(
            list(output_path = str_glue("repli/{name}/{repli_target}.bam"))
          ),
        {
          run(
            "bash",
            command_line
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
    )
  ),

  # From BAM pileups: Create an initial GRanges for each sample. It uses bins of
  # read-centers (non-overlapping bins), no markdup, MAPQ >= 20.
  tar_map(
    cross_join(
      bowtie.refs,
      dplyr::rename(repli.samples, repli_target = "name")
    ) %>%
      mutate(
        suffix = str_glue("{repli_target}_{name}"),
        bulk_reads_misc = rlang::syms(str_glue("bulk_reads_misc_repli.bam_{repli_target}_{name}")),
        bulk_reads_idxstats = rlang::syms(str_glue("bulk_reads_idxstats_repli.bam_{repli_target}_{name}")),
        bulk_reads_split = mapply(
          \(ref_name, n) rlang::syms(
            str_glue(
              "bulk_reads_{names(chr.lengths)}_{n}"
            )
          ) %>%
            setNames(names(chr.lengths)) %>%
            append(list("list"), .) %>%
            do.call(call, ., quote = T),
          name,
          str_glue("repli.bam_{repli_target}_{name}"),
          SIMPLIFY = F
        ),
        chic.tile.diameter_1000 = str_glue("chic.tile.diameter_1000_{name}") %>%
          rlang::syms()
      ) %>%
      rowwise() %>%
      # Include 2L_Histone_Repeat_Unit - split refs and sum up - when splitting by rname and summing.
      mutate(
        bulk_reads_shortrefs = if (name == "masked") {
          list(call("rbind", bulk_reads_misc, rlang::sym(str_glue("bulk_reads_2L_Histone_Repeat_Unit_repli.bam_{repli_target}_masked"))))
        } else {
          list(bulk_reads_misc)
        },
        markdup = isTRUE(is_paired_end),
        df_pos_midpoint_callable = list(
          if (isTRUE(is_paired_end)) quote(paired_end_pos_to_midpoint()) else quote(single_end_pos_to_midpoint())
        ),
        count_overlaps_bases_callable = list(
          if (markdup) quote(count_overlaps_bases) else quote(count_overlaps_bases_no_markdup)
        )
      ),
    names = suffix,
    tar_target(
      repli.granges,
      tibble(
        reads = append(bulk_reads_split, list(bulk_reads_misc)) %>%
          sapply(
            \(df) df %>%
              df_pos_midpoint_callable,
            simplify=F
          ) %>%
          setNames(NULL),
        tile_granges = reads %>%
          sapply(
            \(df) chic.tile.diameter_1000[
              if (!nrow(df))
                integer(0)
              else if (df$rname[1] %in% names(masked.lengths))
                seqnames(chic.tile.diameter_1000) == df$rname[1]
              else
                !(seqnames(chic.tile.diameter_1000) %in% names(masked.lengths))
            ]
          ),
        partition_tiles = stopifnot(sum(sapply(tile_granges, length)) == length(chic.tile.diameter_1000)),
        granges = mapply(
          count_overlaps_bases_callable, tile_granges, reads
        )
      ) %>%
        with(
          granges %>%
            GRangesList %>%
            unlist %>%
            `metadata<-`(value = list(est_library_size = sum(sapply(reads, nrow))))
        )
    )
  ),

  # From each sample GRanges: Apply regression.
  tar_map(
    mutate(
      rowwise(
        cross_join(
          rbind(
            experiment.driver,
            tibble(driver = c("Kc167", "S2")) %>%
              mutate(celltype = driver)
          ),
          dplyr::rename(bowtie.refs, reference = "name")
        )
      ),
      samples = list(parse(text = deparse(subset(repli.samples, genotype == driver)))),
      sample_granges = rlang::syms(
        with(
          subset(repli.samples, genotype == driver),
          paste0(str_glue("repli.granges_{name}_{reference}"))
        )
      ) %>%
        list(),
      chic.track = rlang::syms(str_glue("chic.tile.diameter_40_score_{reference}")),
      chic.results = list(
        call(
          "c",
          H3K4 = call("chic_fseq_l2fc", rlang::sym(str_glue("chic.bw.tracks_H3K4_{celltype}_CN_{reference}"))),
          H3K27 = call("chic_fseq_l2fc", rlang::sym(str_glue("chic.bw.tracks_H3K27_{celltype}_CN_{reference}"))),
          H3K9 = call("chic_fseq_l2fc", rlang::sym(str_glue("chic.bw.tracks_H3K9_{celltype}_CN_{reference}")))
        )
      )
    ),
    names = celltype | reference,
    tar_target(
      repli.experiment,
      as_bulk_summarized_experiment(
        sample_granges,
        colData = mutate(
          as_tibble(eval(samples)),
          rep = as.factor(rep),
          rep = {
            contrasts(rep) <- contr.helmert(length(levels(rep)))
            rep
          }
        )
      ) %>%
        # As we will use Beta dist Bayesian Inference, create the metadata ahead
        # of time.
        init_beta_dm_experiment
    ),
    tar_target(
      repli.glm,
      glm_gp(
        repli.experiment,
        ~ 0 + full + as.factor(rep),
        offset = repli.experiment@metadata$offset,
        size_factors = 1,
        verbose = TRUE
      )
    ),
    tar_target(
      repli.tracks,
      tiles_to_fseq(repli.glm, "full", levels(repli.experiment$full), repli.experiment@metadata$granges, bw = 2000)
    ),
    tar_file(
      repli.bw,
      export(
        GRanges(
          unlist(repli.tracks),
          score = (
            as.matrix(unlist(repli.tracks)@elementMetadata) %*% repli.exp.contrast.numerator
              / as.matrix(unlist(repli.tracks)@elementMetadata) %*% repli.exp.contrast.denominator
          ) %>%
            as.numeric() %>%
            plogistanh() %>%
            # Correct the timing values so that half of positions are in Q1&Q2
            # (very close to achieving this outcome, when we consider that not
            # every range has the same width - as 2Cen and 3Cen ref seqs do not
            # need sliding windows at all). These are a necessary post-hoc
            # offset applied to the regression models so that our GSC and CySC
            # experiments are comparable.
            `-`(median(., na.rm = T)) %>%
            qlogistanh() %>%
            replace(
              which(
                as.matrix(unlist(repli.tracks)@elementMetadata) %*% repli.exp.contrast.denominator
                < 1
              ),
              0
            )
        ),
        with(list(celltype = celltype, reference = reference), BigWigFile(str_glue("repli/Replication_Value_{celltype}_{reference}.bw")))
      ) %>%
        as.character(),
      packages = c("dplyr", "rtracklayer", "stringr", "tibble", "tidyr")
    ),

    tar_target(
      repli.beta,
      future_sapply(
        seq(nrow(repli.experiment)),
        \(i) tryCatch(beta_dm_regression(repli.experiment, i), error=\(e) beta_dm_regression_calculate_prior()),
        simplify=FALSE
      ),
      packages = c(tar_option_get("packages"), "extraDistr", "future.apply", "mvtnorm", "pracma")
    ),
    tar_target(
      repli.beta.2,
      repli.beta %>%
        split(rowData(repli.experiment)$seqnames) %>%
        sapply(
          \(lst) if (length(lst) > 1) laplace_approx_sliding_transform(lst) else lst,
          simplify=FALSE
        ) %>%
        unlist(rec=FALSE, use.names=FALSE) %>%
        future_sapply(\(elem) tryCatch(laplace_approx_expectation_ratio(elem), error=\(e) NA)) %>%
        `*`(-2) %>%
        `+`(1) %>%
        plogistanh() %>%
        `-`(median(., na.rm = T)) %>%
        qlogistanh() %>%
        GRanges(
          unlist(repli.granges_nos_E3_chr, use.names=FALSE),
          score=.
        ),
      packages = c(tar_option_get("packages"), "future.apply", "extraDistr", "mvtnorm", "pracma")
    ),
    tar_file(
      repli.beta.bw,
      repli.beta.2 %>%
        export(BigWigFile(str_glue("repli/Replication_Bayes_", celltype, "_", reference, ".bw"))) %>%
        as.character
    ),

    # Repli ChIC heatmap
    tar_target(
      repli.value.rank,
      rtracklayer::import(rtracklayer::BigWigFile(repli.bw)) %>%
        approx_track(chic.track) %>%
        rank_track() %>%
        elementMetadata() %>%
        as.data.frame(),
      packages = tar_option_get("packages") %>% c("S4Vectors"),
      format = "parquet"
    ),
    tar_target(
      repli.value.bindata,
      append(
        list(chic.track),
        repli.value.rank %>% as.list()
      ) %>%
        do.call(GRanges, .) %>%
        bin_track_by_rank(250) %>%
        elementMetadata() %>%
        as.data.frame(),
      format = "parquet"
    ),
    tar_target(
      repli.value.projection,
      granges_bin_to_projection(GRanges(chic.track, bin = repli.value.bindata$bin, size_of_bin_bp = repli.value.bindata$size_of_bin_bp))
    ),
    tar_target(
      repli.chic.projection.profile,
      cbind(
        ranking = seq(0, 1, length.out = nrow(repli.value.projection[[1]])),
        repli = apply_granges_projection(repli.value.projection, GRanges(chic.track, score = repli.value.rank$score)),
        sapply(
          chic.results,
          \(f) f %>%
            BigWigFile() %>%
            rtracklayer::import() %>%
            attributes() %>%
            with(
              elementMetadata %>%
                subset(as.logical(seqnames %in% seqlevels(chic.track))) %>%
                as.data.frame() %>%
                pull(score) %>%
                `*`(log(2)) %>%
                exp()
            ) %>%
            GRanges(chic.track, score = .) %>%
            apply_granges_projection(repli.value.projection, .)
        ),
        sample_size_bp = repli.value.bindata$size_of_bin_bp[
          match(seq(nrow(repli.value.projection[[1]])), repli.value.bindata$bin)
        ]
      ) %>%
        matrix(
          nrow = nrow(.),
          dimnames = list(timing = NULL, series = colnames(.))
        )
    ),
    tar_file(
      fig.repli.chic,
      save_figures(
        paste0("figure/", celltype),
        ".pdf",
        tribble(
          ~name, ~figure, ~width, ~height,
          paste0("Repli-CHIC-Whole-Genome-", reference),
          plot_repli_track_raster(repli.chic.projection.profile),
          6,
          4.5
        )
      )
    )
  ),
  tar_target(
    repli.beta.prior.draws,
    beta_prior_draws()
  )
)
