rbind_rle_lists <- function(chic.replicates) {
  chic.table <- chic.replicates %>%
    sapply(
      \(rle_list) rle_list[names(chr.lengths)] %>%
        sapply(list, simplify=F) %>%
        as_tibble,
      simplify=F
    ) %>%
    do.call(rbind, .)
}

rle_lists_sd <- function(rles) {
  apply(
    rles %>% rbind_rle_lists,
    2,
    rles_summary
  ) %>%
    RleList
}

rles_summary <- function(rles, summary="sd") {
  summary_fns <- list(sd = rowSds)
  sapply(rles, as.numeric) %>% (summary_fns[[summary]]) %>% as("Rle")
}