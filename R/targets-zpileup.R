bulk.samples <- rbind(
  reframe(
    chic.samples,
    bam = "chic.bam",
    symbol = sample
  ),
  reframe(
    repli.samples,
    bam = "repli.bam",
    symbol = name
  )
)

bulk.bam <- cross_join(
  bulk.samples,
  tibble(reference = c("chr", "masked"))
) %>%
  mutate(
    bam.target = rlang::syms(str_glue("{bam}_{symbol}_{reference}"))
  )

targets.bulk.samples <- tar_map(
  bulk.bam,
  names = any_of(c("bam", "symbol", "reference")),
  tar_map(
    tibble(refseq = names(masked.lengths)),
    tar_target(
      bulk_reads,
      bam_to_df(bam.target, refseq),
      format = "parquet",
      cue = tar_cue("never"),
      packages = c("GenomicRanges", "Rsamtools", "stringr", "tibble")
    )
  ),
  tar_target(
    bulk_reads_idxstats,
    run(
      "samtools",
      c(
        "idxstats",
        bam.target
      )
    )$stdout %>%
      textConnection %>%
      read.table(sep="\t", col.names=c("rname", "rlength", "mapped_unique_reads", "unmapped_unique_reads")) %>%
      as_tibble,
    format = "parquet",
    cue = tar_cue("never"),
    packages = c("magrittr", "processx", "tibble")
  ),
  tar_target(
    bulk_reads_misc,
    sapply(
      setdiff(bulk_reads_idxstats$rname, names(masked.lengths)),
      \(n) bam_to_df(bam.target, n),
      simplify=FALSE
    ) %>%
      bind_rows,
    format = "parquet",
    cue = tar_cue("never"),
    packages = c("dplyr", "GenomicRanges", "Rsamtools", "stringr", "tibble")
  )
)