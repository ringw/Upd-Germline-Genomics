targets.flybase <- list(
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
  # old filename format:
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
  tar_file(
    dm6.bowtie,
    tibble(
      input_file=flybase.fa,
      output_dir="references/bowtie2",
      output_pattern="dm6",
      bowtie2_arg=paste0(output_dir, "/", output_pattern),
      glob_arg=paste0(output_dir, "/", output_pattern, "/*.bt2")
    ) %>%
      mutate(
        create_output_dir = dir.create(output_dir, showW=FALSE),
        delete_bowtie_files = file.remove(Sys.glob(glob_arg)) %>% list,
        run_bowtie2 = processx::run(
          "bowtie2-build",
          c(
            input_file %>% append(list(sep=",")) %>% do.call(paste, .),
            bowtie2_arg
          )
        ) %>% list,
        glob_result = list(Sys.glob(glob_arg))
      ) %>%
      pull(glob_result) %>%
      unlist
  ),
  tar_map(
    tribble(
      ~name, ~filename,
      "flybase", "dmel-all-chromosome-r6.47.fasta.gz",
      "masked", "RepeatMasker/dmel-r6.47-2L-histone-locus-all-te-masked.fasta.gz",
      "histone_unit", "histone_unit.fa",
      "transposon_sequence_set", "transposon_sequence_set.fa.gz"
    ),
    names = name,
    tar_file(fasta, paste("references", filename, sep="/"))
  ),
  tar_map(
    tribble(
      ~name, ~refs,
      "chr", rlang::syms("fasta_flybase"),
      "masked",
      rlang::syms(c("fasta_masked", "fasta_histone_unit", "fasta_transposon_sequence_set"))
    ),
    names = name,
    tar_file(
      bowtie,
      tibble(
        output_dir="references/bowtie2",
        output_pattern=name,
        bowtie2_arg=str_glue("{output_dir}/{output_pattern}/{output_pattern}"),
        glob_arg=str_glue("{output_dir}/{output_pattern}/{output_pattern}*.bt2"),
      ) %>%
        mutate(
          create_output_dir = dir.create(str_glue("{output_dir}/{output_pattern}"), rec=T, showW=FALSE),
          delete_bowtie_files = file.remove(Sys.glob(glob_arg)) %>% list,
          run_bowtie2 = processx::run(
            "bowtie2-build",
            c(
              refs %>% append(list(sep=",")) %>% do.call(paste, .),
              bowtie2_arg
            )
          ) %>% list,
          glob_result = list(Sys.glob(glob_arg))
        ) %>%
        pull(glob_result) %>%
        unlist
    )
  )
)