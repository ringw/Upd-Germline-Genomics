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
              "bulk_reads_{if (ref_name == 'masked') names(masked.lengths) else names(chr.lengths)}_{n}"
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
        markdup = isTRUE(is_paired_end),
        df_pos_midpoint_callable = list(
          if (isTRUE(is_paired_end)) quote(paired_end_pos_to_midpoint) else quote(single_end_pos_to_midpoint)
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
              df_pos_midpoint_callable(),
            simplify = F
          ) %>%
          setNames(NULL),
        tile_granges = reads %>%
          sapply(
            \(df) chic.tile.diameter_1000[
              if (!nrow(df)) {
                integer(0)
              } else if (df$rname[1] %in% names(masked.lengths)) {
                seqnames(chic.tile.diameter_1000) == df$rname[1]
              } else {
                !(seqnames(chic.tile.diameter_1000) %in% names(masked.lengths))
              }
            ]
          ),
        partition_tiles = stopifnot(sum(sapply(tile_granges, length)) == length(chic.tile.diameter_1000)),
        granges = mapply(
          count_overlaps_bases_callable, tile_granges, reads
        )
      ) %>%
        with(
          granges %>%
            GRangesList() %>%
            unlist() %>%
            `metadata<-`(value = list(est_library_size = sum(sapply(reads, nrow))))
        )
    )
  ),
  # Histone-mapped reads analysis: These are put into a specific
  # SummarizedExperiment that gets analyzed with 5.1 kb window (otherwise we
  # used as little input as possible for Repliseq analysis and arrived at a
  # 3 kb window).
  tar_target(
    repli.tile.diameter_1000_histone_repeat_unit,
    GRanges(
      "2L_Histone_Repeat_Unit",
      IRanges(
        start = seq(1, 4001, by=1000),
        end = c(
          seq(1000, 4000, by=1000),
          setNames(masked.feature.lengths["2L_Histone_Repeat_Unit"], NULL)
        )
      ),
      seqlengths = masked.feature.lengths
    )
  ),
  tar_map(
    mutate(
      rowwise(dplyr::rename(repli.samples, repli_target = "name")),
      bulk_reads_2L_Histone_Repeat_Unit = rlang::syms(str_glue("bulk_reads_2L_Histone_Repeat_Unit_repli.bam_{repli_target}_masked")),
      markdup = isTRUE(is_paired_end),
      df_pos_midpoint_callable = list(
        if (isTRUE(is_paired_end)) quote(paired_end_pos_to_midpoint) else quote(single_end_pos_to_midpoint)
      ),
      count_overlaps_bases_callable = list(
        if (markdup) quote(count_overlaps_bases) else quote(count_overlaps_bases_no_markdup)
      )
    ),
    names = repli_target,
    tar_target(
      histone.repeat.unit.granges,
      count_overlaps_bases_callable(
        repli.tile.diameter_1000_histone_repeat_unit,
        df_pos_midpoint_callable(bulk_reads_2L_Histone_Repeat_Unit)
      )
    )
  ),
  # Sliding plate in regression. We need multiple independent observations of
  # replication coverage which are not too overdispersed.
  tar_target(repli.sliding.weights, rep(1, 3)),

  # Mesh grid to be used for Repliseq inference (Bayesian regression).
  tar_target(
    repli.polar.coordinates,
    meshgrid(
      x = seq(0, pi/2, length.out = 75),
      y = seq(0.1, 10, by=0.1)
    ),
    packages = "pracma"
  ),
  tar_target(
    repli.prior.distribution,
    repli.polar.coordinates %>% beta_dm_prior_filter()
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
      chic.tile.diameter_1000 = rlang::syms(
        str_glue("chic.tile.diameter_1000_", reference)
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
        init_beta_dm_experiment()
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
      repli.posterior.unscaled,
      analyze_repli_experiment(
        repli.experiment,
        repli.sliding.weights,
        repli.polar.coordinates,
        repli.prior.distribution,
        xform_scale = 1,
        xform_center = 0
      ),
      packages = tar_option_get("packages") %>% c("extraDistr", "future.apply", "pracma"),
      format = "parquet"
    ),
    tar_target(
      repli.posterior.unscaled.logit,
      GRanges(
        chic.tile.diameter_1000,
        score = repli.posterior.unscaled %>%
          summarise(
            value = quantify_repli_experiment(prob, prior = repli.prior.distribution$X),
            .by = "rowname"
          ) %>%
          pull(value) %>%
          qlogis()
      )
    ),
    # As an alternative to the Beta quantiles regression, we will fit the %
    # enrichment of each fraction using DM likelihood with no Beta function.
    # Simple maximum a posteriori estimator with df = 4.
    tar_target(
      repli.dm,
      sapply(
        split(
          seq(nrow(repli.experiment)),
          rowData(repli.experiment)$seqnames
        ),
        \(inds) future_lapply(
          seq_along(inds),
          dm_regression_gaussian_plate,
          exper = repli.experiment[inds, ],
          wts = repli.sliding.weights
        ),
        simplify = FALSE
      ) %>%
        unlist(rec = FALSE),
      packages = c(tar_option_get("packages"), "extraDistr", "future.apply", "mvtnorm", "pracma")
    ),
    tar_target(
      repli.posterior.xform.centering,
      c(
        scale=sd(repli.posterior.unscaled.logit$score),
        center=mean(repli.posterior.unscaled.logit$score)
      )
    ),
    tar_target(
      repli.posterior,
      analyze_repli_experiment(
        repli.experiment,
        repli.sliding.weights,
        repli.polar.coordinates,
        repli.prior.distribution,
        theta = repli.posterior.unscaled$theta[
          match(seq(nrow(repli.experiment)), repli.posterior.unscaled$rowname)
        ],
        xform_scale = 1,
        xform_center = repli.posterior.xform.centering["center"]
      ),
      packages = tar_option_get("packages") %>% c("extraDistr", "future.apply", "pracma"),
      format = "parquet"
    ),
    # Apply Bayesian estimator to the posterior and create GRanges track.
    tar_target(
      repli.timing,
      repli.posterior %>%
        summarise(
          value = quantify_repli_experiment(prob, prior = repli.prior.distribution$X),
          .by = "rowname"
        ) %>%
        with(value[match(seq_along(chic.tile.diameter_1000), rowname)]) %>%
        repli_logistic_beta_to_tanh() %>%
        GRanges(chic.tile.diameter_1000, score = .) # %>%
        # repli_logistic_link_apply_scale(
        #   num_fractions = repli.experiment@metadata$beta_regression_num_fractions
        # )
        ,
      packages = c(tar_option_get("packages"), "future.apply", "extraDistr", "mvtnorm", "pracma")
    ),
    tar_file(
      repli.bw,
      repli.timing$score %>%
        replace(is.na(.) | !is.finite(.), 0) %>%
        GRanges(chic.tile.diameter_1000, score = .) %>%
        export(BigWigFile(str_glue("repli/Replication_Bayes_", celltype, "_", reference, ".bw"))) %>%
        as.character()
    ),
    tar_file(
      repli.map.bw,
      GRanges(
        chic.tile.diameter_1000,
        score = repli.posterior %>%
          group_by(rowname) %>%
          summarise(
            value = mean(
              (1 - 2*value)[rank(prob, ties="max") == length(prob)]
            )
          ) %>%
          pull(value)
      ) %>%
        export(BigWigFile(str_glue("repli/Replication_MAP_", celltype, "_", reference, ".bw"))) %>%
        as.character()
    ),
    tar_target(
      repli.autocorrelation.beta,
      repli.timing %>% autocorrelate_centered_granges(1001)
    )
  ),
  # From each masked reference: Apply 2L_Histone_Repeat_Unit
  tar_map(
    mutate(
      rowwise(
        tibble(
          rbind(
            experiment.driver,
            tibble(driver = c("Kc167", "S2")) %>%
              mutate(celltype = driver)
          ),
          reference = "masked",
        )
      ),
      samples = list(parse(text = deparse(subset(repli.samples, genotype == driver)))),
      sample_granges = rlang::syms(
        with(
          subset(repli.samples, genotype == driver),
          paste0(str_glue("histone.repeat.unit.granges_{name}"))
        )
      ) %>%
        list(),
      repli.experiment = rlang::syms(str_glue("repli.experiment_{celltype}_{reference}")),
      repli.posterior.unscaled = rlang::syms(str_glue("repli.posterior.unscaled_{celltype}_{reference}")),
      repli.posterior.xform.centering = rlang::syms(str_glue("repli.posterior.xform.centering_{celltype}_{reference}")),
    ),
    names = celltype | reference,
    tar_target(
      repli.experiment_2L_Histone_Repeat_Unit,
      {
        exper <- as_bulk_summarized_experiment(
          sample_granges,
          colData = mutate(as_tibble(eval(samples)), rep = rep %>% as.factor())
        ) %>%
          # As we will use Beta dist Bayesian Inference, create the metadata ahead
          # of time.
          init_beta_dm_experiment()
        exper@metadata$beta_regression_size_factors <- repli.experiment@metadata$beta_regression_size_factors
        exper
      }
    ),
    tar_target(
      repli.posterior_2L_Histone_Repeat_Unit,
      beta_dm_regression_post_probability(
        repli.experiment_2L_Histone_Repeat_Unit,
        3,
        # 5 kb of observations!
        rep(1, 5),
        repli.polar.coordinates,
        repli.prior.distribution,
        xform_scale = 1,
        xform_center = repli.posterior.xform.centering["center"]
      ),
      packages = tar_option_get("packages") %>% c("extraDistr", "future.apply", "pracma")
    ),
    tar_target(
      repli.likelihood_2L_Histone_Repeat_Unit,
      beta_dm_regression_likelihood(
        repli.experiment_2L_Histone_Repeat_Unit,
        3,
        # 5 kb of observations!
        rep(1, 5),
        repli.polar.coordinates,
        xform_scale = 1,
        xform_center = repli.posterior.xform.centering["center"]
      ),
      packages = tar_option_get("packages") %>% c("extraDistr", "future.apply", "pracma")
    ),
    tar_target(
      repli.mmse_2L_Histone_Repeat_Unit,
      with(
        repli.polar.coordinates %>%
          append(
            list(d_timing_d_angle = 1/(1 + sin(2 * .$X)))
          ),
        c(
          weighted.mean(
            sin(X)/(sin(X)+cos(X)),
            repli.likelihood_2L_Histone_Repeat_Unit * d_timing_d_angle
          ),
          weighted.mean(
            Y,
            repli.likelihood_2L_Histone_Repeat_Unit * repli.polar.coordinates$Y * repli.prior.distribution$Y
          )
        )
      )
    ),
    tar_target(
      repli.glm_2L_Histone_Repeat_Unit,
      SummarizedExperiment(
        list(
          counts = repli.experiment_2L_Histone_Repeat_Unit %>%
            assay() %>%
            apply(1, identity) %>%
            matrix(nrow = 1)
        ),
        colData = repli.experiment_2L_Histone_Repeat_Unit %>%
          colData() %>%
          with(
            tibble(
              full = full %>% rep(length(repli.tile.diameter_1000_histone_repeat_unit)),
              size_factors = repli.experiment %>%
                assay() %>%
                colSums() %>%
                `*`(1 / 1000 / 1000) %>%
                rep(length(repli.tile.diameter_1000_histone_repeat_unit)),
            )
          )
      ) %>%
        glm_gp(
          ~ 0 + full,
          size_factors = .$size_factors,
          verbose = TRUE
        )
    )
  ),
  tar_file(
    fig.histone.timing,
    save_figures(
      "figure/Both-Cell-Types",
      ".pdf",
      tibble(
        rowname="Repli-Histone-Locus-Posterior-Dist",
        figure=bind_rows(
          list(
            Germline=repli.posterior_2L_Histone_Repeat_Unit_Germline_masked,
            Somatic=repli.posterior_2L_Histone_Repeat_Unit_Somatic_masked,
            Kc167=repli.posterior_2L_Histone_Repeat_Unit_Kc167_masked,
            S2=repli.posterior_2L_Histone_Repeat_Unit_S2_masked
          ),
          .id = "celltype"
        ) %>%
          mutate(celltype = celltype %>% factor(str_to_title(names(repli_posterior_bar_colors)))) %>%
          plot_posterior(repli.polar.coordinates) %>%
          `+`(
            theme(
              plot.margin = margin(0, 15.5, 0, 5.5)
            )
          ) %>%
          list(),
        width=4,
        height=3
      )
    )
  ),
  tar_file(
    fig.model.illustration,
    save_figures(
      "figure/Both-Cell-Types",
      ".pdf",
      tibble(
        rowname="Repli-Model-Illustration",
        figure=plot_fpkm_bayes_mmse(
          fpkm_bars = repli.glm_2L_Histone_Repeat_Unit_Somatic_masked$Beta %>% as.numeric() %>% exp(),
          bayes_mmse_param = repli.mmse_2L_Histone_Repeat_Unit_Somatic_masked,
          tibble(
            celltype = "Somatic",
            repli.posterior_2L_Histone_Repeat_Unit_Somatic_masked
          ),
          repli.polar.coordinates
        ) %>%
          list(),
        width=9.5,
        height=2.5
      )
    ),
    packages = tar_option_get("packages") %>% c("egg")
  ),
  tar_target(
    histone.timing.bayes.factor,
    sapply(
      list(
        Somatic=repli.posterior_2L_Histone_Repeat_Unit_Somatic_masked,
        Kc167=repli.posterior_2L_Histone_Repeat_Unit_Kc167_masked,
        S2=repli.posterior_2L_Histone_Repeat_Unit_S2_masked
      ),
      \(df) bayes_factor_repli_experiment(
        tibble(
          Germline = repli.posterior_2L_Histone_Repeat_Unit_Germline_masked$prob,
          Other = df$prob,
        ),
        repli.prior.distribution$X
      )
    )
  ),

  # Repli ChIC heatmap
  tar_map(
    mutate(
      rowwise(
        cross_join(
          experiment.driver, dplyr::rename(bowtie.refs, reference = "name")
        )
      ),
      repli.experiment = rlang::syms(str_glue("repli.experiment_{celltype}_{reference}")),
      repli.timing = rlang::syms(str_glue("repli.timing_{celltype}_{reference}")),
      repli.targets = repli.samples$name %>%
        subset(repli.samples$genotype == driver) %>%
        list(),
      repli.idxstats = rlang::syms(
        str_glue("bulk_reads_idxstats_repli.bam_{repli.targets}_{reference}")
      ) %>%
        list(),
      chic.track = rlang::syms(str_glue("chic.tile.diameter_500_score_{reference}")),
      chic.results = list(
        call(
          "setNames",
          rlang::syms(
            str_glue("chic.experiment.quantify_H3K{c(4,27,9)}_{celltype}_peakcalling.broad_{reference}")
          ),
          as.character(str_glue("H3K{c(4,27,9)}"))
        )
      ),
      chromosome_arms_diameter_500_score = rlang::syms(
        str_glue("chromosome_arms_diameter_500_score_{reference}")
      ),
      chromosome_arms_diameter_1000 = rlang::syms(
        str_glue("chromosome_arms_diameter_1000_{reference}")
      ),
      quartile.factor = rlang::syms(str_glue("quartile.factor_{celltype}")),
      celltype_lower = tolower(celltype),
    ),
    names = celltype | reference,
    tar_target(
      repli.value.bindata,
      repli.timing %>%
        approx_track(chic.track) %>%
        cut_track(seq(-1, 1, by = 0.002)) %>%
        elementMetadata() %>%
        as.data.frame(),
      format = "parquet"
    ),
    tar_target(
      repli.value.projection,
      granges_bin_to_projection(GRanges(chic.track, bin = repli.value.bindata$bin, size_of_bin_bp = repli.value.bindata$size_of_bin_bp))
    ),
    tar_target(
      repli.chromosome.arms.factor_diameter_500_score,
      repli_chromosome_arms_factor_scaffolds(chromosome_arms_diameter_500_score)
    ),
    tar_target(
      repli.chromosome.arms.factor_diameter_1000,
      repli_chromosome_arms_factor_scaffolds(chromosome_arms_diameter_1000)
    ),
    tar_target(
      repli.chromosome.arms.profile,
      repli.value.bindata$bin %>%
        factor(seq(max(repli.value.bindata$bin))) %>%
        split(
          repli.chromosome.arms.factor_diameter_500_score,
          .
        ) %>%
        sapply(
          \(chr) chr %>%
            table() %>%
            prop.table()
        ) %>%
        t() %>%
        `colnames<-`(value = str_glue("chr{colnames(.)}"))
    ),
    tar_target(
      repli.chic.projection.profile,
      cbind(
        repli = seq(-1, 1, by = 0.002),
        sapply(
          chic.results,
          \(f) f %>%
            attributes() %>%
            with(
              elementMetadata %>%
                subset(as.logical(seqnames %in% seqlevels(chic.track)), select = c(1, 2)) %>%
                as.data.frame() %>%
                `colnames<-`(value = c("H3", "mark")) %>%
                with((mark / H3) %>% replace(!is.finite(.), 1))
            ) %>%
            GRanges(chic.track, score = .) %>%
            apply_granges_projection(repli.value.projection, .)
        ),
        split(
          quartile.factor[names(chic.tile.diameter_500_score_genes)],
          factor(
            repli.value.bindata$bin[chic.tile.diameter_500_score_genes$lookup],
            seq(max(repli.value.bindata$bin))
          )
        ) %>%
          sapply(table) %>%
          `/`(rowMaxs(.)) %>%
          t %>%
          matrix(
            ncol = 4,
            dimnames = list(
              NULL,
              c("TSS_off", "TSS_low", "TSS_medium", "TSS_high")
            )
          ),
        # Expression pattern is not actually examined yet.
        split(
          expression.factor[names(chic.tile.diameter_500_score_genes)],
          factor(
            repli.value.bindata$bin[chic.tile.diameter_500_score_genes$lookup],
            seq(max(repli.value.bindata$bin))
          )
        ) %>%
          sapply(table) %>%
          `/`(rowMaxs(.)) %>%
          t %>%
          matrix(
            ncol = 4,
            dimnames = list(
              NULL,
              c("TSS_GSC", "TSS_CySC", "TSS_both", "TSS_neither")
            )
          ),
        repli.chromosome.arms.profile,
        sample_size_bp = repli.value.bindata$size_of_bin_bp[
          match(seq(nrow(repli.value.projection[[1]])), repli.value.bindata$bin)
        ]
      ) %>%
        matrix(
          nrow = nrow(.),
          dimnames = list(timing = NULL, series = colnames(.))
        )
    ),
    tar_target(
      repli.mode.chr.weights,
      weight_repli_mode_prop_table(repli.experiment, repli.chromosome.arms.factor_diameter_1000)
    ),
    tar_file(
      fig.repli.chic,
      save_figures(
        paste0("figure/", celltype),
        ".pdf",
        tribble(
          ~name, ~figure, ~width, ~height,
          paste0("Repli-CHIC-Whole-Genome-", reference),
          plot_repli_track_raster(
            repli.chic.projection.profile,
            log2_limits = c(-0.65, 0.8),
            repli.mode.chr.weights = repli.mode.chr.weights
          ) %>%
            grob_add_arm_colors_legend(w = unit(1, "null"), h = unit(4, "in")),
          6,
          4.5,
          paste0("Repli-CHIC-Whole-Genome-", reference, "-Naive-Mode"),
          plot_repli_track_raster(
            repli.chic.projection.profile,
            log2_limits = c(-0.65, 0.8),
            repli.mode.chr.weights = 1
          ) %>%
            grob_add_arm_colors_legend(w = unit(1, "null"), h = unit(4, "in")),
          6,
          4.5
        )
      ),
      packages = tar_option_get("packages") %>% c("colorspace", "cowplot", "grid", "gtable")
    ),
    tar_file(
      fig.repli.chic.raster,
      save_png(
        str_glue("figure/", celltype, "/Repli-CHIC-Whole-Genome-", reference, ".png"),
        plot_repli_track_raster(
          repli.chic.projection.profile,
          log2_limits = c(-0.65, 0.8),
          repli.mode.chr.weights = repli.mode.chr.weights
        ) %>%
          as_grob() %>%
          `$`(grobs) %>%
          `[[`(value = 6) %>%
          plot_grid(),
        w = 750,
        h = 5
      ),
      packages = tar_option_get("packages") %>% c("colorspace")
    )
  ),
  tar_target(
    repli.posterior_chr,
    with(
      repli.posterior_Germline_chr,
      tibble(
        rowname,
        sapply(
          list(
            Germline = repli.posterior_Germline_chr,
            Somatic = repli.posterior_Somatic_chr,
            Kc167 = repli.posterior_Kc167_chr,
            S2 = repli.posterior_S2_chr
          ),
          \(df) df$prob
        ) %>%
          as_tibble()
      )
    ),
    format = "parquet"
  ),
  tar_target(
    repli.bayes.factor_chr,
    chic.tile.diameter_1000_chr %>%
      `elementMetadata<-`(
        value = repli.posterior_chr %>%
          repli_posterior_bayes_factor_all_pairs(
            prior = repli.prior.distribution$X
          )
      )
  ),
  tar_target(
    repli.peaks_chr,
    mapply(
      repli_timing_nested_peak_calling,
      elementMetadata(repli.bayes.factor_chr),
      repli_timing_factor_all_pairs(
        list(
          Germline = repli.timing_Germline_chr,
          Somatic = repli.timing_Somatic_chr,
          Kc167 = repli.timing_Kc167_chr,
          S2 = repli.timing_S2_chr
        )
      ),
      c(
        100,
        100,
        100,
        100,
        100,
        1e4
      )
    )
  ),
  tar_target(
    diff.replication.progression.gene,
    with(
      read.csv(assay.data.sc),
      names(repli.peaks_chr$Germline_Somatic)[
        findOverlaps(
          GRanges(
            chr %>%
              replace(
                is.na(chr) | !(chr %in% levels(droplevels(seqnames(repli.peaks_chr$Germline_Somatic)))),
                "*"
              ) %>%
              factor(levels = c(seqlevels(repli.peaks_chr$Germline_Somatic), "*")),
            IRanges(
              start = ifelse(
                strand == "+", start, end
              ) %>%
                replace(is.na(chr), 1),
              width = 1
            )
          ),
          repli.peaks_chr$Germline_Somatic
        ) %>%
          sapply(\(vec) vec[1])
      ]
    )
  ),

  # Replication main-figure gene-wise analysis
  tar_target(
    genewise,
    detail_repli_analysis(
      assay.data.sc,
      flybase.sequence.ontology,
      repli.timing_Germline_chr,
      repli.timing_Somatic_chr,
      repli.bayes.factor_chr,
      Upd_cpm[, "germline"],
      Upd_cpm[, "somatic"],
      Upd_regression_somatic,
      chic.gene.enrichment,
      chromosome_pericetromere_label
    ),
    packages = tar_option_get("packages") %>% c("tidyr")
  ),

  # Repliseq dmel-all-chromosomes supplementary data
  tar_file(
    sd_repliseq,
    publish_repli_analysis(
      assay.data.sc,
      repli.timing_Germline_chr,
      repli.timing_Somatic_chr,
      repli.bayes.factor_chr,
      diff.replication.progression.gene,
      Upd_cpm[, "germline"],
      Upd_cpm[, "somatic"],
      chromosome_pericetromere_label,
      "Supplemental_Data/SD04_Repliseq.xlsx"
    ),
    packages = tar_option_get("packages") %>% c("tidyr")
  ),
  tar_file(
    repliseq_bed,
    tibble(
      output_file = "repli/Diff_Rep_Program_Assay.bed",
      do_write = with_options(
        list(scipen=100),
        write.table(
          reframe(
            rownames_to_column(as.data.frame(repli.peaks_chr$Germline_Somatic)),
            seqnames,
            start - 1,
            end,
            rowname,
            score = ((NegDiff + PosDiff) / 2) %>%
              pmax(-1) %>%
              pmin(1) %>%
              `+`(1) %>%
              `*`(1000 / 2) %>%
              round
          ),
          output_file,
          quote = F,
          sep = "\t",
          row.names = F,
          col.names = F
        )
      )
    )$output_file
  ),
  tar_file(
    repliseq_bed_earlier,
    tibble(
      output_file = "repli/Diff_Rep_Germline_Earlier.bed",
      do_write = with_options(
        list(scipen=100),
        write.table(
          read.table(repliseq_bed) %>%
            subset(grepl("GermlineEarlier", V4)),
          output_file,
          quote = F,
          sep = "\t",
          row.names = F,
          col.names = F
        )
      )
    )$output_file
  ),
  tar_file(
    repliseq_bed_later,
    tibble(
      output_file = "repli/Diff_Rep_Germline_Later.bed",
      do_write = with_options(
        list(scipen=100),
        write.table(
          read.table(repliseq_bed) %>%
            subset(grepl("GermlineLater", V4)),
          output_file,
          quote = F,
          sep = "\t",
          row.names = F,
          col.names = F
        )
      )
    )$output_file
  ),

  # For each repli track, look up gene TSS values.
  tar_target(
    chic.tile.diameter_1000_genes,
    GRanges(
      seqnames(tss_location),
      ranges(tss_location),
      seqinfo = seqinfo(chic.tile.diameter_1000_chr),
      lookup = tss_location %>%
        findOverlaps(chic.tile.diameter_1000_chr) %>%
        sapply(\(v) head(c(v, NA), 1))
    ) %>%
      setNames(names(tss_location))
  ),
  tar_target(
    chic.tile.diameter_500_score_genes,
    GRanges(
      seqnames(tss_location),
      ranges(tss_location),
      seqinfo = seqinfo(chic.tile.diameter_500_score_chr),
      lookup = tss_location %>%
        findOverlaps(chic.tile.diameter_500_score_chr) %>%
        sapply(\(v) head(c(v, NA), 1))
    ) %>%
      setNames(names(tss_location))
  ),
  tar_target(
    diff.timing.factor,
    rep("null", length(chic.tile.diameter_500_score_genes)) %>%
      replace(
        !is.na(chic.tile.diameter_500_score_genes$lookup),
        findOverlaps(
          chic.tile.diameter_500_score_chr[chic.tile.diameter_500_score_genes$lookup %>% subset(!is.na(.))],
          repli.peaks_chr$Germline_Somatic
        ) %>%
          as.list() %>%
          sapply(
            \(peak_inds) ifelse(
              grepl("GermlineEarlier", names(repli.peaks_chr$Germline_Somatic)[peak_inds]),
              1,
              -1
            ) %>%
              median()
          )
      ) %>%
      replace(is.na(.), "null") %>%
      factor() %>%
      fct_recode(GermlineEarlier="1", GermlineLater="-1") %>%
      fct_relevel("null", after=4) %>%
      setNames(names(chic.tile.diameter_500_score_genes))
  ),

  # Repli graphic for dmel-all-chromosomes, not the masked bowtie reference.
  tar_file(
    fig.repli.skew,
    save_figures(
      "figure/Both-Cell-Types",
      ".pdf",
      tribble(
        ~rowname, ~figure, ~width, ~height,
        "Repli-Skew-Diff-Timing",
        plot_genomic_ranges_score(
          GRanges(
            chic.tile.diameter_1000_chr,
            score = repli.timing_Germline_chr$score - repli.timing_Somatic_chr$score
          )[
            to(findOverlaps(repli.peaks_chr$Germline_Somatic, chic.tile.diameter_1000_chr))
          ],
          repli_early_late_background$E, repli_early_late_background$L,
          name = "Diff.",
          limits=c(-1.8, 1.8),
          breaks=c(-1.8, 0, 1.8)
        ),
        5.75,
        4
      )
    ),
    packages = tar_option_get("packages") %>% c("grid", "gtable")
  ),
  tar_target(
    histone.locus,
    with(
      read.csv(assay.data.sc),
      GRanges(
        "2L",
        IRanges(
          start = start[match("His1:CG33801", X)],
          end = end[match("His3:CG33866", X)]
        )
      )
    )
  ),
  tar_file(
    fig.repli.timing,
    save_figures(
      "figure/Both-Cell-Types",
      ".pdf",
      tribble(
        ~rowname, ~figure, ~width, ~height,
        "Repli-Timing",
        plot_track_2score(
          GRanges(repli.timing_Germline_chr, score_1=repli.timing_Germline_chr$score, score_2=repli.timing_Somatic_chr$score),
          "transparent",
          "transparent",
          score_1 = "#6D9965",
          score_2 = "#B377B0",
          background_color = "#EAEAEB",
          line_color = "#010101",
          roi = tibble(
            as_tibble(repli.peaks_chr$Germline_Somatic),
            sgn = grepl("Earlier", names(repli.peaks_chr$Germline_Somatic)),
            chr = seqnames,
            xmin = start,
            xmax = end,
            ymin = -Inf,
            ymax = Inf,
            fill = c(
              "FALSE" = "#461B56",
              "TRUE" = "#255529"
            )[
              as.character(sgn)
            ]
          ) %>%
            arrange(fill)
        ),
        5.75,
        4,
      ) %>%
        mutate(
          figure = {
            # Chr 2: after the X (second element of TableGrob)
            figure[[1]]$grobs[[2]] <- figure[[1]]$grobs[[2]] %>%
              gtable_add_grob(
                linesGrob(
                  rep(mid(histone.locus) / chr.lengths["2L"], 2),
                  c(0, 1),
                  gp=gpar(col = "yellow")
                ),
                t = 7,
                l = 5
              )
            figure
          }
        )
    ),
    packages = tar_option_get("packages") %>% c("grid", "gtable")
  ),
  tar_map(
    tibble(
      experiment.driver,
      repli.dm = rlang::syms(str_glue("repli.dm_{celltype}_chr")),
      repli.timing = rlang::syms(str_glue("repli.timing_{celltype}_chr")),
      figure_dir = str_glue("figure/{celltype}"),
      repli_quartile_data_TSS = rlang::syms(str_glue("repli_quartile_data_{celltype}_TSS")),
      repli_quartile_active_data_TSS = rlang::syms(str_glue("repli_quartile_active_data_{celltype}_TSS")),
      celltype_lower = tolower(celltype),
    ),
    names = celltype,
    tar_file(
      fig.repli.skew,
      save_figures(
        str_glue("figure/", celltype),
        ".pdf",
        tribble(
          ~rowname, ~figure, ~width, ~height,
          "Repli-Skew-Chrs",
          plot_track(repli.timing, repli_early_late_background$E, repli_early_late_background$L),
          5.75,
          4
        )
      ),
      packages = tar_option_get("packages") %>% c("grid", "gtable")
    ),

    # Factor for genes.
    tar_target(
      repli.gene,
      repli.timing$score[chic.tile.diameter_1000_genes$lookup] %>%
        cut(c(-Inf, -0.5, 0, 0.5, Inf)) %>%
        `levels<-`(
          value = c("L", "ML", "EM", "E")
        ) %>%
        setNames(names(chic.tile.diameter_1000_genes))
    ),
    # In cutting genes, exclude the mitochondrion_genome at this point. If the
    # criteria are timing values then it doesn't matter, but the later exclusion
    # of the mitochondrion_genome is going to cause "quartiles" of different
    # sizes.
    tar_target(
      repli.gene.active,
      repli.timing$score[chic.tile.diameter_1000_genes$lookup] %>%
        replace(
          Upd_cpm[names(chic.tile.diameter_1000_genes), celltype_lower] < 5 |
            !as.logical(
              seqnames(chic.tile.diameter_1000_genes) %in% names(chr.lengths)
            ),
          NA
        ) %>%
        rank(na="keep", ties="first") %>%
        cut(c(0, quantile(., c(1/4, 1/2, 3/4), na.rm=T), Inf)) %>%
        `levels<-`(
          value = c("LA", "MLA", "EMA", "EA")
        ) %>%
        setNames(names(chic.tile.diameter_1000_genes))
    ),
    tar_file(
      repli.chic.results,
      save_figures(
        figure_dir,
        ".pdf",
        tribble(
          ~name, ~figure, ~width, ~height,
          "CHIC-TSS-AllMarks-Repli-Timing-Quartile",
          chic_plot_average_profiles_facet_grid(
            mutate(
              repli_quartile_data_TSS,
              genes = repli %>%
                fct_relabel(
                  \(n) paste0(
                    n,
                    " (n = ",
                    sapply(
                      n,
                      \(n) repli_quartile_data_TSS$n[
                        head(
                          which(
                            repli_quartile_data_TSS$repli == n
                          ),
                          1
                        )
                      ]
                    ),
                    ")"
                  )
                ),
              .keep = "unused"
            ),
            "Timing",
            unlist(rev(repli_level_colors), use.names = FALSE),
            c(0.4, 0.62, 0.73, 0.85),
            facet_wrap(vars(mark))
          ) +
            coord_cartesian(NULL, chic_average_profile_limits, ex=F),
          10, 3.25,
          "CHIC-TSS-AllMarks-Repli-Timing-Quartile-RNAseq-Active",
          chic_plot_average_profiles_facet_grid(
            mutate(
              repli_quartile_active_data_TSS %>%
                subset(!is.na(repli.active)),
              genes = repli.active %>%
                fct_relabel(
                  \(n) paste0(
                    n,
                    " (n = ",
                    sapply(
                      n,
                      \(n) repli_quartile_active_data_TSS$n[
                        head(
                          which(
                            repli_quartile_active_data_TSS$repli.active == n
                          ),
                          1
                        )
                      ]
                    ),
                    ")"
                  )
                ),
              .keep = "unused"
            ),
            "Timing",
            unlist(rev(repli_level_colors), use.names = FALSE),
            c(0.4, 0.62, 0.73, 0.85),
            facet_wrap(vars(mark))
          ) +
            coord_cartesian(NULL, chic_average_profile_limits, ex=F),
          10, 3.25,
        )
      )
    ),
    tar_file(
      fig.repli.dm.barchart,
      save_figures(
        paste0("figure/", celltype),
        ".pdf",
        tribble(
          ~rowname, ~figure, ~width, ~height,
          "Repli-Fraction-Bars-Chrs",
          dm_barplot(
            repli.dm,
            chromosome_arms_diameter_1000
          ) %>%
            set_panel_size(
              w = unit(3.25, "in"),
              h = unit(1.35, "in")
            ),
          4,
          2,
          "Repli-Fraction-Bars-Chromosome-Arms",
          dm_barplot(
            repli.dm,
            chromosome_arms_diameter_1000,
            chromosome_arms_diameter_1000$group
          ) %>%
            set_panel_size(
              w = unit(3.25, "in"),
              h = unit(1.35 * 8/7, "in")
            ),
          4,
          2,
        )
      ),
      packages = tar_option_get("packages") %>% c("egg")
    )
  ),

  # Violin graphic of chromosome arm features.
  tar_target(
    repli.feature.factor,
    GRanges(
      chromosome_arms_diameter_1000,
      group = chromosome_arms_diameter_1000$group %>%
        factor(c(levels(.), "4", "X", "Y")) %>%
        replace(which(seqnames(chic.tile.diameter_1000_chr) == "4"), "4") %>%
        replace(which(seqnames(chic.tile.diameter_1000_chr) == "X"), "X") %>%
        replace(which(seqnames(chic.tile.diameter_1000_chr) == "Y"), "Y")
    )
  ),
  tar_target(
    repli.timing.byfeature,
    bind_rows(
      list(
        Germline=tibble(
          chr = factor(seqnames(chic.tile.diameter_1000_chr), names(chr.lengths)),
          feature = repli.feature.factor$group,
          timing = repli.timing_Germline_chr$score
        ),
        Somatic=tibble(
          chr = factor(seqnames(chic.tile.diameter_1000_chr), names(chr.lengths)),
          feature = repli.feature.factor$group,
          timing = repli.timing_Somatic_chr$score
        )
      ),
      .id = "celltype"
    ) %>%
      tibble(
        occur_tss = assay.data.sc %>%
          read.csv() %>%
          with(
            GRanges(
              chr %>% replace(is.na(chr), "Y"),
              IRanges(ifelse(strand == "+", start, end) %>% replace(is.na(chr), -1), width = 1)
            )
          ) %>%
          countOverlaps(chic.tile.diameter_1000_chr, .) %>%
          rep(2)
      ) %>%
      subset(!is.na(feature)) %>%
      tibble(
        row = feature %>%
          fct_relabel(
            \(n) ifelse(
              grepl("^2", n),
              "2",
              ifelse(
                grepl("^3", n),
                "3",
                ""
              )
            )
          )
      ),
    format = "parquet"
  ),
  tar_file(
    fig.repli.timing.byfeature,
    save_figures(
      "figure/Both-Cell-Types",
      ".pdf",
      tribble(
        ~rowname, ~figure, ~width, ~height,
        "Repli-Chromosome-Arm-Violin",
        plot_repli_timing_byfeature(repli.timing.byfeature),
        6,
        3,
        "Repli-Chromosome-Arm-TSS-Violin",
        plot_repli_timing_byfeature(repli.timing.byfeature, tss_only = TRUE, box_ci = TRUE),
        6,
        3,
        "Repli-Chromosome-Violin",
        ggplot(
          # Reverse (top to bottom) feature. Reverse again celltype!
          repli.timing.byfeature %>%
            mutate(celltype = factor(celltype, c("Somatic", "Germline"))),
          aes(chr, timing, fill=celltype)
        ) +
          geom_violin() +
          scale_x_discrete(limits = rev) +
          scale_y_reverse() +
          scale_fill_manual(values = cell_type_violin_colors) +
          labs(x = "dm6 Chromosome") +
          coord_flip() +
          theme(
            aspect.ratio = 3,
            panel.grid.major.x = element_blank(),
            legend.position = "none"
          ),
        2,
        6,
      )
    )
  ),

  # Repli-CHIC graphic for one CHIC track.
  tar_map(
    tibble(
      chic.experiments,
      celltype = name,
      repli.timing = rlang::syms(str_glue("repli.timing_{celltype}_chr")),
      chic.bw.tracks = rlang::syms(str_glue("chic.bw.tracks_{mark}_{celltype}_CN_chr")),
      chic.experiment.quantify_peakcalling.broad_chr = rlang::syms(str_glue("chic.experiment.quantify_{mark}_{celltype}_peakcalling.broad_chr")),
    ),
    names = mark | celltype,
    # CHIC broad track, currently used for Repli heatmap.
    tar_target(
      chic.enrich,
      GRanges(
        chic.tile.diameter_500_score_chr,
        score = with(
          elementMetadata(chic.experiment.quantify_peakcalling.broad_chr)[, 1:2] %>%
            `colnames<-`(value = c("H3", "chromatinMark")),
          (log(chromatinMark / H3) / log(2)) %>%
            replace(!is.finite(.), NA) %>%
            `-`(
              median(subset(., as.logical(seqnames(chic.tile.diameter_500_chr) %in% c("2L", "2R", "3L", "3R", "4"))), na.rm = T)
            )
        )
      )
    ),
    tar_file(
      fig.repli.chic.raster,
      save_png(
        str_glue("figure/", celltype, "/Repli-CHIC-", mark, ".png"),
        plot_grid(
          # The name of grobs[[6]]$children[[3]] is "panel-1". It is the plot
          # area of the ggplot data, which is exactly a 120x80 raster.
          as_grob(plot_repli_chic_bin2d(repli.timing,
            chic.enrich,
            chic_step_size = 100
          ) +
            theme(panel.ontop = FALSE))$grobs[[6]]$children[[3]]
        ),
        w = 121,
        h = 81
      )
    ),
    tar_file(
      fig.repli.chic,
      save_figures(
        str_glue("figure/", celltype),
        ".pdf",
        tribble(
          ~filename, ~figure, ~width, ~height,
          as.character(str_glue("Repli-CHIC-", mark)),
          plot_repli_chic_bin2d(repli.timing,
            chic.enrich,
            chic_step_size = 100
          ),
          6,
          4
        )
      )
    ),
    tar_file(
      fig.violin.repli.chic,
      save_figures(
        str_glue("figure/", celltype),
        ".pdf",
        tribble(
          ~rowname, ~figure, ~width, ~height,
          as.character(str_glue("Repli-CHIC-", mark, "-Violin")),
          tibble(
            repli = approx_track(repli.timing, chic.enrich)$score,
            timing = structure(
              4 - (as.numeric(cut(repli, c(-Inf, -0.375, 0, 0.375, Inf))) - 1),
              levels = c("E", "EM", "ML", "L"),
              class = "factor"
            ),
            enrichment = chic.enrich$score
          ) %>%
            violin_plot_repli_chic(),
          5, 4.25
        )
      )
    )
  ),
  tar_target(
    repli.beta.prior.draws,
    beta_prior_draws()
  )
)
