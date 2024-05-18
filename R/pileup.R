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

bam_reference_coords <- function(df) {
  df %>%
    mutate(
      cigar = factor(cigar),
      width = cigar %>%
        fct_relabel(\(v) as.character(cigar_ref_length(v))) %>%
        as.numeric,
      pos_crick = pos + width - 1
    )
}

bam_cover_read_bp <- function(df, min_mapq, markdup=FALSE) {
  stopifnot(all(df$rname %in% names(masked.lengths)))
  df <- df %>%
    bam_reference_coords %>%
    filter(between(mapq, min_mapq, 254)) %>%
    reframe(
      rname,
      pos = pos + round((pos_crick - pos) / 2),
      # TODO: Fix bug where bam_to_df with no reads might have NULL for the
      # dc column.
      dc = if (markdup || !("dc" %in% colnames(df))) 1 else dc
    ) %>%
    group_by(rname, pos) %>%
    summarise(nreads = sum(dc), .groups="drop") %>%
    arrange(rname, pos)
  if (!is.factor(df$rname)) df$rname <- factor(df$rname)
  mapply(
    \(rname, i, x) sparseVector(x, i, length = masked.lengths[rname]),
    levels(df$rname),
    split(df$pos, df$rname),
    split(df$nreads, df$rname),
    SIMPLIFY=FALSE
  )
}

count_overlaps_sparse_vectors <- function(sparse_vectors, tile_width) {
  tiles <- mapply(
    \(n, v) str_glue("{n}:1-{length(v)}"),
    names(sparse_vectors),
    sparse_vectors
  ) %>%
    GRanges %>%
    slidingWindows(width = as.integer(tile_width), step = as.integer(tile_width)) %>%
    unlist
  counts <- mapply(
    \(n, v) GRanges(seqnames = if (length(v@i)) n, ranges = if (length(v@i)) IRanges(v@i, width = 1) else IRanges(), score = v@x),
    names(sparse_vectors),
    sparse_vectors,
    SIMPLIFY=FALSE
  ) %>%
    GRangesList %>%
    unlist
  hits = findOverlaps(tiles, counts) %>% as("List")
  GRanges(tiles, score = sum(extractList(counts$score, hits)), seqlengths=sapply(sparse_vectors, length))
}