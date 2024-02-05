make_chr_lengths <- function(fa_table) {
  all_features = NULL
  feature_name = NULL
  feature_length = 0
  save_feature <- function(all_features) all_features %>%
    c(setNames(feature_length, feature_name))
  for (line in fa_table[, 1]) {
    if (str_starts(line, ">")) {
      if (!is.null(feature_name)) all_features <- all_features %>% save_feature
      feature_name <- line %>% str_extract(">(\\S+)", group=1)
      feature_length <- 0
    }
    else feature_length <- feature_length + line %>%
      str_extract_all("[ACGTN]") %>%
      sapply(length)
  }
  all_features %>% save_feature
}