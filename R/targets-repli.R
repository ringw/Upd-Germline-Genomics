repli.samples = read.csv('repli/repli_samples.csv')
repli.samples$replication_value = repli.samples$replication_value %>% factor(unique(.))
repli.samples$full = repli.samples$replication_value
levels(repli.samples$full) <- c("Early", "Early-Mid", "Mid-Late", "Late")
repli.samples$abbrev = repli.samples$replication_value
levels(repli.samples$abbrev) <- c("E", "G", "J", "L")
repli.samples <- repli.samples %>%
  mutate(name = paste0(genotype, "_", abbrev, rep))

repli.contrasts <- repli.samples %>%
  subset(select=c(full, replication_value)) %>%
  subset(!duplicated(.)) %>%
  within(replication_value <- as.numeric(as.character(replication_value))) %>%
  as_tibble

repli.exp.contrast.numerator <- c(0.875, 0.625, 0.375, 0.125)
repli.exp.contrast.denominator <- c(1, 1, 1, 1)

targets.repli <- list(
  # FASTQ files: Align to BAM and count lines in FASTQ.
  tar_map(
    cross_join(
      bowtie.refs,
      dplyr::rename(repli.samples, repli_target="name")
    ) %>%
      rowwise %>%
      mutate(
        command_line = if (isTRUE(is_paired_end))
            list(
              call(
                "c",
                "-i",
                quote(align_chic_lightfiltering),
                substitute(str_replace(ref_paths[1], "\\..*", ""), list(ref_paths=bowtie)),
                str_glue("{filename}_R1_001.fastq.gz"),
                str_glue("{filename}_R2_001.fastq.gz"),
                quote(output_path)
              )
            )
          else
            list(
              call(
                "c",
                "-i",
                quote(align_repli_lightfiltering),
                substitute(str_replace(ref_paths[1], "\\..*", ""), list(ref_paths=bowtie)),
                str_glue("Upd_Tumor/Repli/{filename}"),
                quote(output_path)
              )
            )
        ),
    names = repli_target | name,
    tar_file(
      repli.bam,
      with(
        list(name=name, repli_target=repli_target) %>%
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
    ),
    tar_target(
      repli.readcount,
      as.integer(
        run(
          "bash",
          c(
            "-c",
            paste0(
              "rclone cat sharepoint:'Bio Data'/Upd_Tumor/Repli/",
              filename,
              " | gunzip -c | wc -l"
            )
          )
        )$stdout
      ),
      packages = "processx"
    )
  ),

  # From BAM pileups: Create an initial GRanges for each sample. It uses bins of
  # read-centers (non-overlapping bins), no markdup, MAPQ >= 20.
  tar_map(
    cross_join(
      bowtie.refs,
      dplyr::rename(repli.samples, repli_target="name")
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
            do.call(call, ., quote=T),
          name,
          str_glue("repli.bam_{repli_target}_{name}"),
          SIMPLIFY=F
        )
      ) %>%
      rowwise %>%
      # Include 2L_Histone_Repeat_Unit - split refs and sum up - when splitting by rname and summing.
      mutate(
        bulk_reads_shortrefs = if (name == "masked")
          list(call("rbind", bulk_reads_misc, rlang::sym(str_glue("bulk_reads_2L_Histone_Repeat_Unit_repli.bam_{repli_target}_masked"))))
        else
          list(bulk_reads_misc)
      ),
    names = suffix,
    tar_target(
      repli.granges,
      GRanges(
        seqnames = levels(bulk_reads_misc$rname),
        IRanges(
          start = 1,
          width = pull(bulk_reads_idxstats, "rlength", "rname")[levels(bulk_reads_misc$rname)]
        ),
        score = bulk_reads_shortrefs %>%
          subset(between(mapq, 20, 254)) %>%
          split(.$rname) %>%
          sapply(\(df) sum(df$dc)),
        seqlengths = pull(bulk_reads_idxstats, "rlength", "rname")[levels(bulk_reads_misc$rname)]
      ) %>%
        split(seqnames(.)) %>%
        replace(
          match(names(bulk_reads_split), names(.)),
          GRangesList(
            bulk_reads_split %>%
            mapply(
              \(n, df) df %>%
                mutate(rname = factor(rname, levels=n)) %>%
                bam_cover_read_bp(min_mapq = 20, markdup = FALSE) %>%
                count_overlaps_sparse_vectors(tile_width=1000L),
              names(.),
              .,
              SIMPLIFY=F
            ) %>%
              unlist
          )
        )
    )
  ),

  # From each sample GRanges: Apply regression.
  tar_map(
    mutate(
      rowwise(
        cross_join(
          experiment.driver,
          dplyr::rename(bowtie.refs, reference="name")
        )
      ),
      samples = list(parse(text = deparse(subset(repli.samples, genotype == driver)))),
      sample_granges = rlang::syms(
        with(
          subset(repli.samples, genotype == driver),
          paste0(str_glue("repli.granges_{name}_{reference}"))
        )
      ) %>%
        list,
      chic.track = rlang::syms(str_glue("chic.track_{reference}")),
      chic.results = rlang::syms(str_glue("chic.results_{driver}_{reference}")),
    ),
    names = driver | reference,
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
      )
    ),
    tar_target(
      repli.glm,
      glm_gp(
        repli.experiment,
        ~ 0 + full + as.factor(rep),
        offset = repli.experiment@metadata$offset,
        size_factors = 1,
        overdispersion = "global",
        verbose = TRUE
      )
    ),
    tar_target(
      repli.tracks,
      tiles_to_fseq(repli.glm, "full", levels(repli.experiment$full), repli.experiment@metadata$granges, bw=2000)
    ),
    # tar_target(
    #   repli.exp.contrast.numerator,
    #   left_join(
    #     tibble(coef = colnames(repli.glm$Beta)),
    #     reframe(repli.contrasts, coef = str_glue("full{full}"), replication_value),
    #     "coef"
    #   ) %>%
    #     deframe %>%
    #     replace_na(0.),
    #   packages = c("dplyr", "stringr", "tibble", "tidyr")
    # ),
    # tar_target(
    #   repli.exp.contrast.denominator,
    #   left_join(
    #     tibble(coef = colnames(repli.glm$Beta)),
    #     reframe(repli.contrasts, coef = str_glue("full{full}"), value = 1),
    #     "coef"
    #   ) %>%
    #     deframe %>%
    #     replace_na(0.),
    #   packages = c("dplyr", "stringr", "tibble", "tidyr")
    # ),
    tar_file(
      repli.bw,
      export(
        GRanges(
          unlist(repli.tracks),
          score = (
            as.matrix(unlist(repli.tracks)@elementMetadata) %*% repli.exp.contrast.numerator
            / as.matrix(unlist(repli.tracks)@elementMetadata) %*% repli.exp.contrast.denominator
          ) %>%
            as.numeric %>%
            replace(
              which(
                as.matrix(unlist(repli.tracks)@elementMetadata) %*% repli.exp.contrast.denominator
                < 1
              ),
              0
            )
        ),
        with(list(driver=driver, reference=reference), BigWigFile(str_glue("repli/Replication_Value_{driver}_{reference}.bw")))
      ) %>%
        as.character,
      packages = c("dplyr", "rtracklayer", "stringr", "tibble", "tidyr")
    ),

    # Repli ChIC heatmap
    tar_target(
      repli.value.rank,
      rtracklayer::import(rtracklayer::BigWigFile(repli.bw)) %>%
        approx_track(chic.track) %>%
        rank_track %>%
        elementMetadata %>%
        as.data.frame,
      packages = tar_option_get("packages") %>% c("S4Vectors"),
      format = "parquet"
    ),
    tar_target(
      repli.value.bindata,
      append(
        list(chic.track),
        repli.value.rank %>% as.list
      ) %>%
        do.call(GRanges, .) %>%
        bin_track_by_rank(250) %>%
        elementMetadata %>%
        as.data.frame,
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
          \(f) f %>% BigWigFile %>% rtracklayer::import() %>% attributes %>%
            with(
              elementMetadata %>%
                subset(as.logical(seqnames %in% seqlevels(chic.track))) %>%
                as.data.frame %>%
                pull(score)
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
  )
)
