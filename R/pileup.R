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
    \(rname, i, x) sparseVector(x, i, length = masked.lengths[as.character(rname)]),
    levels(df$rname),
    split(df$pos, df$rname),
    split(df$nreads, df$rname),
    SIMPLIFY=FALSE
  )
}

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

paired_end_pos_to_midpoint <- function(df) {
  df %>%
    paired_end_reads_to_fragment_lengths %>%
    mutate(pos = pos + round((fragment_end_crick - pos) / 2))
}

bam_cover_paired_end_fragments_bp <- function(df, min_mapq, min_fl, max_fl, markdup=FALSE) {
  stopifnot(all(df$rname %in% names(masked.lengths)))
  df$rname <- df$rname %>% droplevels()

  df <- df %>% paired_end_reads_to_fragment_lengths
  df <- df %>%
    bam_reference_coords %>%
    filter(between(length, min_fl, max_fl), between(mapq, min_mapq, 254)) %>%
    reframe(
      rname,
      pos = pos + round((fragment_end_crick - pos) / 2),
      # TODO: Fix bug where bam_to_df with no reads might have NULL for the
      # dc column.
      dc = if (markdup || !("dc" %in% colnames(df))) 1 else dc
    ) %>%
    group_by(rname, pos) %>%
    summarise(nreads = sum(dc), .groups="drop") %>%
    arrange(rname, pos)
  if (!is.factor(df$rname)) df$rname <- factor(df$rname)
  mapply(
    \(rname, i, x) sparseVector(x, i, length = masked.lengths[as.character(rname)]),
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

count_overlaps_bases <- function(windows, bases_df) {
  bases_df <- bases_df %>% group_by(rname, pos) %>% summarise(count = length(pos), .groups="drop")
  gcts <- bases_df %>% with(GRanges(rname, IRanges(pos, width = 1), score = count))
  hits <- findOverlaps(windows, gcts)
  GRanges(windows, score = sum(extractList(gcts$score, hits)))
}