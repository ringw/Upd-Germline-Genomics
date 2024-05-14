pileup_scan_bam_param <- function(...) ScanBamParam(
  mapqFilter=20,
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

bam_to_df <- function(bam_file, rname, ...) {
  granges <- GRanges(str_glue("{rname}:1-{max(chr.lengths)}"))
  scan_bam_param <- pileup_scan_bam_param(which=granges, ...)
  bam_obj <- scanBam(bam_file, param=scan_bam_param)[[1]]
  with(
    bam_obj,
    with(
      bam_obj$tag,
      tibble(
        qname,
        flag,
        rname,
        strand,
        pos,
        qwidth,
        mapq,
        cigar,
        dc
      )
    )
  )
}