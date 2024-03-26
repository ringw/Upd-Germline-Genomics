# Log-scale rounding. With digits=1, this can round to pretty values 0.1 and 1
# and the values 0.2, 0.3, ..., 0.9 each round to unique values in the interval.
# The same is true for the values 0.01, 0.02, ..., 0.09, with digits=2.
# Eventually for very large values, we will not have decimal-level precision.
# This helps with creating an Rle track that is slightly compressed. When we
# will take the quotient of two Rle tracks later, then this is helpful so that
round_log_scale_pretty <- function(v, digits=2) {
  v %>%
    log %>%
    `*`(1/log(10)) %>%
    round(digits=digits) %>%
    `*`(log(10)) %>%
    exp
}

# For each rle in the rle list, compress it (using log-scale under the hood) if
# the length is greater than one.
rle_list_round_log_scale_pretty <- function(rle_list, digits=2) {
  contents <- rle_list %>%
    sapply(\(v) if (length(v) > 1) round_log_scale_pretty(v, digits=digits) else v)
  new_rle_list <- RleList(contents)
  for (n in names(attributes(rle_list)))
    if (is(attr(rle_list, n), "RleList") || is.numeric(attr(rle_list, n)))
      attr(new_rle_list, n) <- attr(rle_list, n)
  new_rle_list
}