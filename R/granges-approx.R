granges_approx <- function(
    query_windows,
    subject_windows) {
  seqlevels(query_windows) <- seqlevelsInUse(query_windows)
  subject_windows <- subject_windows[seqnames(subject_windows) %in% seqlevels(query_windows)]
  seqlevels(subject_windows) <- seqlevels(query_windows)
  query_order <- order(seqnames(query_windows))
  GRanges(
    query_windows,
    score =
      do.call(
        c,
        mapply(
          \(query, subject) approx(
            mid(subject),
            subject$score,
            xout = mid(query)
          )$y,
          split(query_windows, seqnames(query_windows)),
          split(subject_windows, seqnames(subject_windows)),
          SIMPLIFY = FALSE,
          USE.NAMES = FALSE
        )
      )[query_order]
  )
}
