repli.samples = read.csv('repli/repli_samples.csv')
repli.samples$replication_value = repli.samples$replication_value %>% factor(unique(.))
repli.samples$full = repli.samples$replication_value
levels(repli.samples$full) <- c("Early", "Early-Mid", "Mid-Late", "Late")
repli.samples$abbrev = repli.samples$replication_value
levels(repli.samples$abbrev) <- c("E", "G", "J", "L")
repli.samples <- repli.samples %>%
  mutate(name = paste0(genotype, "_", abbrev, rep))

targets.repli <- list(
  # FASTQ files: Align to BAM and count lines in FASTQ.
  tar_map(
    bowtie.refs,
    names = name,
    tar_map(
      dplyr::rename(repli.samples, repli_target="name"),
      names = repli_target,
      tar_file(
        repli.bam,
        with(
          list(name=name, repli_target=repli_target) %>%
            with(
              list(output_path = str_glue("repli/{name}/{repli_target}.bam"), filename=filename)
            ),
          {
            run(
              "bash",
              c(
                "-i",
                align_repli_lightfiltering,
                str_replace(bowtie[1], "\\..*", ""),
                str_glue("Upd_Tumor/Repli/{filename}"),
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
              "bulk_reads_{names(if (ref_name == 'masked') masked.lengths else chr.lengths)}_{n}"
            )
          ) %>%
            setNames(names(if (ref_name == 'masked') masked.lengths else chr.lengths)) %>%
            append(list("list"), .) %>%
            do.call(call, ., quote=T),
          name,
          str_glue("repli.bam_{repli_target}_{name}"),
          SIMPLIFY=F
        )
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
        score = bulk_reads_misc %>%
          subset(between(mapq, 20, 254)) %>%
          split(.$rname) %>%
          sapply(\(df) sum(df$dc))
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
  )
)