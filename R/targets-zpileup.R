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
      format = "parquet"
    )
  )
)