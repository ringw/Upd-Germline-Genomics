# Data frame to help with Targets build targets ----
# Bowtie refs: Pattern of 2 bowtie2 references with or without masking sequences of interest.
# The pipeline is to be run with each bowtie2 reference
bowtie.refs <- tribble(
  ~name, ~bowtie, ~refs,
  "chr",
  rlang::sym("bowtie_chr"),
  rlang::syms("fasta_flybase"),
  "masked",
  rlang::sym("bowtie_masked"),
  rlang::syms(c("fasta_masked", "fasta_histone_unit", "fasta_transposon_sequence_set"))
)
