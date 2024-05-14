repli.samples = read.csv('repli/repli_samples.csv')
repli.samples$replication_value = repli.samples$replication_value %>% factor(unique(.))
repli.samples$full = repli.samples$replication_value
levels(repli.samples$full) <- c("Early", "Early-Mid", "Mid-Late", "Late")
repli.samples$abbrev = repli.samples$replication_value
levels(repli.samples$abbrev) <- c("E", "G", "J", "L")
repli.samples <- repli.samples %>%
  mutate(name = paste0(genotype, "_", abbrev, rep))

targets.repli <- tar_map(
  bowtie.refs,
  names = name,
  tar_map(
    rename(repli.samples, name="repli_target"),
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
      packages = c(
        "dplyr",
        "processx",
        "stringr"
      )
    )
  )
)