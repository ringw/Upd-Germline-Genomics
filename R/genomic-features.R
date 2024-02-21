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

# Create a factor-valued track for the genome. The levels are categories of loci
# so that we can show a bar chart of peak locations.
factor_genome <- function(gtf_path, metafeatures_path) {
  metafeatures <- metafeatures_path %>% read.csv(row.names=1)
  gtf <- read.table(
    gtf_path,
    sep = '\t',
    col.names = c('chr', 'source', 'type', 'start', 'end', 'sc', 'strand', 'fr', 'annotation'),
    header = F,
    quote = ''
  ) %>%
    subset(chr %in% c("2L", "2R", "3L", "3R", "4", "X", "Y"))
  gtf$gene_id <- gtf$annotation %>% str_extract(
    'gene_id "([^"]+)"',
    group = 1
  )
  gtf <- gtf %>% left_join(
    metafeatures %>%
      subset(!duplicated(flybase)) %>%
      mutate(gene_id = flybase, symbol = rownames(.), .keep = "none"),
    by = "gene_id"
  )
  chr_data <- gtf %>% split(gtf$chr) %>%
    sapply(\(df) df %>% split(df$strand), simplify=FALSE)
  rles <- chr_data %>%
    mapply(
      \(n, chr) {
        # Initialize chr to "intergene"
        call_chr <- Rle("intergene", chr.lengths[n])
        # Add "gene_body" to paper over later with other features.
        call_chr[
          chr$`+` %>%
            subset(type == "gene") %>%
            apply(1, \(v) seq(v['start'], v['end']), simplify=FALSE) %>%
            do.call(c, .)
        ] <- "gene_body"
        call_chr[
          chr$`-` %>%
            subset(type == "gene") %>%
            apply(1, \(v) seq(v['end'], v['start']), simplify=FALSE) %>%
            do.call(c, .)
        ] <- "gene_body"
        call_chr
      },
      names(.),
      .,
      SIMPLIFY = FALSE
    )
}