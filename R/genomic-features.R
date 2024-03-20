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
  chr_data <- chr_data[names(chr.lengths)]
  rles <- chr_data %>%
    mapply(
      \(n, chr) {
        # Initialize chr to "intergene"
        call_chr <- Rle("intergene", chr.lengths[n])

        # Lay down 5'UTR and 3'UTR, but 5'UTR may be overwritten as part of the
        # "promoter".
        call_chr[
          rbind(chr$`+`, chr$`-`) %>%
            subset(type == "5UTR") %>%
            apply(
              1,
              \(v) seq(v['start'], v['end']),
              simplify=FALSE
            ) %>%
            do.call(c, .)
        ] <- "5'UTR"
        call_chr[
          rbind(chr$`+`, chr$`-`) %>%
            subset(type == "3UTR") %>%
            apply(
              1,
              \(v) seq(v['start'], v['end']),
              simplify=FALSE
            ) %>%
            do.call(c, .)
        ] <- "3'UTR"

        # Per IRanges::promoters method, the genomic feature extracted for
        # potential promoters will start 2000 bp upstream and end 200 bp
        # downstream of TSS.
        call_chr[
          chr$`+` %>%
            subset(type == "gene") %>%
            apply(
              1,
              \(v) seq(
                max(as.numeric(v['start']) - 2000, 1),
                min(as.numeric(v['start']) + 200, chr.lengths[n])
              ),
              simplify=FALSE) %>%
            do.call(c, .)
        ] <- "promoter"
        call_chr[
          chr$`-` %>%
            subset(type == "gene") %>%
            apply(
              1,
              \(v) seq(
                max(as.numeric(v['end']) - 200, 1),
                min(as.numeric(v['end']) + 2000, chr.lengths[n])
              ),
              simplify=FALSE) %>%
            do.call(c, .)
        ] <- "promoter"

        # Lay down intron first using the hull of all of the CDS sequences. By
        # using CDS, we will not clobber the UTR regions.
        call_chr[
          rbind(chr$`+`, chr$`-`) %>%
            subset(type == "CDS") %>%
            group_by(gene_id) %>%
            summarise(start = min(c(start, end)), end = max(c(start, end))) %>%
            apply(
              1,
              \(v) seq(v['start'], v['end']),
              simplify=FALSE
            ) %>%
            do.call(c, .)
        ] <- "intron"

        # Name CDS as exon. It cuts out 5'UTR and 3'UTR which have already been
        # labeled above.
        call_chr[
          rbind(chr$`+`, chr$`-`) %>%
            subset(type == "CDS") %>%
            apply(
              1,
              \(v) seq(v['start'], v['end']),
              simplify=FALSE
            ) %>%
            do.call(c, .)
        ] <- "exon"
        call_chr %>% factor(
          c("5'UTR", "exon", "intron", "3'UTR", "promoter", "intergene")
        )
      },
      names(.),
      .,
      SIMPLIFY = FALSE
    )
}

classify_feature_overlaps <- function(bed_track, factor_track) {
  mapply(
    \(peaks, rle) {
      peaks <- peaks %>% rowwise %>% mutate(lookup = list(seq(start, end, by = 10)))
      peak_cat <- peaks$lookup %>% do.call(c, .)
      peak_splitter <- peaks$lookup %>%
        mapply(
          \(index, v) rep(index, length(v)),
          seq_along(.),
          .,
          SIMPLIFY = FALSE
        ) %>%
        do.call(c, .)
      categories <- rle[peak_cat] %>% split(peak_splitter)
      peaks %>%
        subset(select=-lookup) %>%
        cbind(
          region = categories %>%
            sapply(\(v) table(v) %>% which.max %>% names) %>%
            factor(levels(factor_track[[1]]))
        )
    },
    (bed_track %>% split(.$chr))[names(chr.lengths)],
    factor_track[names(chr.lengths)],
    SIMPLIFY=FALSE
  ) %>%
    do.call(rbind, .)
}
