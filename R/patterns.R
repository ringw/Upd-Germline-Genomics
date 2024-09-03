bowtie.refs <- tribble(
  ~name, ~bowtie, ~refs,
  "chr",
  rlang::sym("bowtie_chr"),
  rlang::syms("fasta_flybase"),
  # Currently no use for the RepeatMasker & histone repeat unit & transposon
  # sequence set reference!
  # "masked",
  # rlang::sym("bowtie_masked"),
  # rlang::syms(c("fasta_masked", "fasta_histone_unit", "fasta_transposon_sequence_set"))
)
