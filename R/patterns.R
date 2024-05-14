bowtie.refs <- tribble(
  ~name, ~bowtie, ~refs,
  "chr",
  rlang::sym("bowtie_chr"),
  rlang::syms("fasta_flybase"),
  "masked",
  rlang::sym("bowtie_masked"),
  rlang::syms(c("fasta_masked", "fasta_histone_unit", "fasta_transposon_sequence_set"))
)
