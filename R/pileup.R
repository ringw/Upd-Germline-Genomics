# Rsamtools to GRanges of fragment midpoints ----
pileup_scan_bam_param <- function(...) ScanBamParam(
  what=c(
    "qname",
    "flag",
    "rname",
    "strand",
    "pos",
    "qwidth",
    "mapq",
    "cigar"
  ),
  tag="dc",
  ...
)

# Reads a bam file from disk using Rsamtools to data frame.
# All bam files will be indexed, read in using one rname at a time, and written
# to Apache Arrow for streamlining further analysis in pure R.
bam_to_df <- function(bam_file, rname, ...) {
  granges <- GRanges(str_glue("{rname}:1-{max(chr.lengths)}"))
  scan_bam_param <- pileup_scan_bam_param(which=granges, ...)
  bam_obj <- scanBam(bam_file, param=scan_bam_param)[[1]]
  with(
    bam_obj,
    with(
      bam_obj$tag,
      tibble(
        qname = factor(qname, unique(qname)),
        flag,
        rname,
        strand,
        pos,
        qwidth,
        mapq,
        cigar,
        dc
      ) %>%
        arrange(qname)
    )
  )
}

# Applies our indexed bam to df for a list of reference sequences to read all of them.
bam_to_df_multi_ref <- function(bam_file, rname, ...) {
  granges <- GRanges(rname, IRanges(1, width = rep(max(chr.lengths), length(rname))))
  scan_bam_param <- pileup_scan_bam_param(which=granges, ...)
  bam_results <- scanBam(bam_file, param=scan_bam_param)
  bam_results_df <- bam_results %>% as_tibble %>% as.data.frame
  rownames(bam_results_df) <- names(bam_results[[1]])
  bam_obj <- apply(bam_results_df, 1, \(lst) do.call(c, setNames(lst, NULL)))
  bam_obj$tag <- list(
    dc = do.call(
      c,
      bam_obj$tag[names(bam_obj$tag) == "dc"] %>%
        setNames(NULL)
    )
  )
  with(
    bam_obj,
    with(
      bam_obj$tag,
      tibble(
        qname = factor(qname, unique(qname)),
        flag,
        rname,
        strand,
        pos,
        qwidth,
        mapq,
        cigar,
        dc
      ) %>%
        arrange(qname)
    )
  )
}

# Watson and Crick end of Sequence Alignment (SAM) information. We need to parse
# the CIGAR alignment string to determine the Crick end of the alignment,
# because just the Watson coordinate is given to us.
cigar_ref_length <- function(cigars) {
  cigars %>%
    sapply(
      \(cigar) cigar %>%
        gregexpr("[0-9]+[MNDI]", .) %>%
        `[[`(1) %>%
        mapply(
          # Insertion is the only case where we are not counting reference base pairs.
          \(pos, length) if(substring(cigar, pos + length - 1, pos + length - 1) != "I")
            as.numeric(
              substring(cigar, pos, pos + length - 2)
            )
            else 0,
          .,
          attr(., "match.length")
        ) %>%
        sum
    )
}

# Given our BAM data frame, determine the alignment width on the genome using the cigar string.
# Written to a "width" column.
bam_reference_coords <- function(df) {
  df %>%
    mutate(
      cigar = factor(cigar),
      width = cigar %>%
        fct_relabel(\(v) as.character(cigar_ref_length(v))) %>%
        as.character %>%
        as.numeric,
      pos_crick = pos + width - 1
    )
}

# Index into a sorted paired end reads df to take the first read, and produce
# the fragment length of the proper pair for further fragment binning.
# Reduces the length of the data frame by half from proper pair reads to fragments.
# All bash scripts for bowtie2 paired-end alignment filter by proper pairs before
# writing to disk, so this is safe for this genomics project.
paired_end_reads_to_fragment_lengths <- function(df) {
  stopifnot(all(df$strand[seq(2, nrow(df), by=2)] == "-"))
  df$fragment_end_crick <- rep(
    bam_reference_coords(df[seq(2, nrow(df), by=2),, drop=F])$pos_crick,
    each=2
  )
  df <- df[seq(1, nrow(df), by=2),, drop=F]
  df$length <- df$fragment_end_crick - df$pos + 1
  df
}

# Adds the alignment length column from the cigar and then produces a GRanges
# with the fragment lengths in this alignment.
paired_end_reads_to_granges <- function(df, ...) {
  if (!("length" %in% colnames(df)))
    df <- df %>% paired_end_reads_to_fragment_lengths
  GRanges(df$rname, IRanges(start = df$pos, width = df$length), ...)
}

# For a fragment ends version of our analysis, replace the pos with the 5 prime
# end of the alignment if the read is reverse-aligned.
paired_end_pos_to_5_prime <- function(df) {
  stopifnot(all(df$strand[seq(2, nrow(df), by=2)] == "-"))
  df <- df %>% bam_reference_coords
  df$fragment_end_crick <- rep(
    df$pos_crick[seq(2, nrow(df), by=2)],
    each=2
  )
  df$length <- df$fragment_end_crick - rep(
    df$pos[seq(1, nrow(df), by=2)],
    each=2
  ) + 1
  df %>%
    mutate(
      pos = ifelse(strand == "+", pos, pos_crick)
    )
}

# Cuts the length of the proper pairs set of alignment in half and replaces all
# "pos" with the point halfway between the proper pairs, for each fragment.
paired_end_pos_to_midpoint <- function(df) {
  df %>%
    paired_end_reads_to_fragment_lengths %>%
    mutate(pos = pos + round((fragment_end_crick - pos) / 2))
}

# Given single-end reads, always take the full length of a single read to be good
# enough and calculate the point 50% between the start and end of the read. This is
# expected to be reasonably close to the fragment midpoint at the much larger scale
# of sliding windows applied for Repliseq. For sharp quantification (ChIC-seq), we
# expect that the technology will be paired-end not single-end.
single_end_pos_to_midpoint <- function(df) {
  df %>%
    bam_reference_coords() %>%
    mutate(pos = mid(IRanges(pos, width=width)))
}

# Runs find overlaps. Ignores the markdup "duplicate count" field so the markdup that
# we have run before writing to disk is effective at taking reads down to a single
# fragment being counted in case of duplicates.
count_overlaps_bases <- function(windows, bases_df) {
  bases_df <- bases_df %>% group_by(rname, pos) %>% summarise(count = length(pos), .groups="drop")
  gcts <- bases_df %>% with(GRanges(rname, IRanges(pos, width = 1), score = count))
  hits <- findOverlaps(windows, gcts)
  GRanges(windows, score = sum(extractList(gcts$score, hits)))
}

# Reads out the markdup "duplicate count" when it comes to finding overlaps. We ignore
# the fact that we have run samtools markdup, and in the fragment counts, we bring those
# duplicate counts back into the quantification to count all of the reads/fragments that
# have the duplicated alignment.
count_overlaps_bases_no_markdup <- function(windows, bases_df) {
  bases_df <- bases_df %>% group_by(rname, pos) %>% summarise(count = sum(dc), .groups="drop")
  gcts <- bases_df %>% with(GRanges(rname, IRanges(pos, width = 1), score = count))
  hits <- findOverlaps(windows, gcts)
  GRanges(windows, score = sum(extractList(gcts$score, hits)))
}