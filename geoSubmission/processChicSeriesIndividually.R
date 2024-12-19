library(Matrix)
library(GenomicRanges)
library(withr)

source("_targets.R")

tar_load(chic.tile.diameter_40_score_chr)

save_bed <- function(name, tracks, filename) {
  data <- elementMetadata(tracks) %>% as.data.frame()
  colnames(data) <- colnames(data) %>%
    str_replace("^score.", paste0(name, "_"))
  data <- tibble(
    chr = as.character(seqnames(tracks)),
    start = -1 + start(chic.tile.diameter_40_score_chr),
    end = end(chic.tile.diameter_40_score_chr),
    data,
  )
  with_options(
    list(scipen=1000),
    write.table(data, filename, row.names=F, sep="\t", na="", quote=F)
  )
  filename <- filename
}

original_rep <- chic.samples %>% subset(driver == "nos" & group == "H3K4") %>% pull(rep)
new_rep <- original_rep %>% factor() %>% `levels<-`(value = seq_along(levels(.)))
H3K4_Germline <- tar_read(chic.experiment.granges_H3K4_Germline_CN_chr) %>%
  split(seqnames(.)) %>%
  sapply(\(gr) ksmooth_sliding_windows(gr, bw = 40)) %>%
  GRangesList() %>%
  unlist() %>%
  `elementMetadata<-`(
    value = elementMetadata(.) %>%
      apply(2, \(v) v %>% round(digits = 3)) %>%
      as.data.frame() %>%
      `colnames<-`(
        value = colnames(.) %>%
          mapply(
            \(n, original, new) n %>% str_replace(original, new),
            .,
            paste0("Rep", as.character(original_rep)),
            paste0("Rep", as.character(new_rep))
          )
      )
  )
save_bed("H3K4_Germline", H3K4_Germline, "geoSubmission/H3K4_Germline.tsv")
original_rep <- chic.samples %>% subset(driver == "nos" & group == "H3K27") %>% pull(rep)
new_rep <- original_rep %>% factor() %>% `levels<-`(value = seq_along(levels(.)))
H3K27_Germline <- tar_read(chic.experiment.granges_H3K27_Germline_CN_chr) %>%
  split(seqnames(.)) %>%
  sapply(\(gr) ksmooth_sliding_windows(gr, bw = 40)) %>%
  GRangesList() %>%
  unlist() %>%
  `elementMetadata<-`(
    value = elementMetadata(.) %>%
      apply(2, \(v) v %>% round(digits = 3)) %>%
      as.data.frame() %>%
      `colnames<-`(
        value = colnames(.) %>%
          mapply(
            \(n, original, new) n %>% str_replace(original, new),
            .,
            paste0("Rep", as.character(original_rep)),
            paste0("Rep", as.character(new_rep))
          )
      )
  )
save_bed("H3K27_Germline", H3K27_Germline, "geoSubmission/H3K27_Germline.tsv")
original_rep <- chic.samples %>% subset(driver == "nos" & group == "H3K9") %>% pull(rep)
new_rep <- original_rep %>% factor() %>% `levels<-`(value = seq_along(levels(.)))
H3K9_Germline <- tar_read(chic.experiment.granges_H3K9_Germline_CN_chr) %>%
  split(seqnames(.)) %>%
  sapply(\(gr) ksmooth_sliding_windows(gr, bw = 40)) %>%
  GRangesList() %>%
  unlist() %>%
  `elementMetadata<-`(
    value = elementMetadata(.) %>%
      apply(2, \(v) v %>% round(digits = 3)) %>%
      as.data.frame() %>%
      `colnames<-`(
        value = colnames(.) %>%
          mapply(
            \(n, original, new) n %>% str_replace(original, new),
            .,
            paste0("Rep", as.character(original_rep)),
            paste0("Rep", as.character(new_rep))
          )
      )
  )
save_bed("H3K9_Germline", H3K9_Germline, "geoSubmission/H3K9_Germline.tsv")
original_rep <- chic.samples %>% subset(driver == "tj" & group == "H3K4") %>% pull(rep)
new_rep <- original_rep %>% factor() %>% `levels<-`(value = seq_along(levels(.)))
H3K4_Somatic <- tar_read(chic.experiment.granges_H3K4_Somatic_CN_chr) %>%
  split(seqnames(.)) %>%
  sapply(\(gr) ksmooth_sliding_windows(gr, bw = 40)) %>%
  GRangesList() %>%
  unlist() %>%
  `elementMetadata<-`(
    value = elementMetadata(.) %>%
      apply(2, \(v) v %>% round(digits = 3)) %>%
      as.data.frame() %>%
      `colnames<-`(
        value = colnames(.) %>%
          mapply(
            \(n, original, new) n %>% str_replace(original, new),
            .,
            paste0("Rep", as.character(original_rep)),
            paste0("Rep", as.character(new_rep))
          )
      )
  )
save_bed("H3K4_Somatic", H3K4_Somatic, "geoSubmission/H3K4_Somatic.tsv")
original_rep <- chic.samples %>% subset(driver == "tj" & group == "H3K27") %>% pull(rep)
new_rep <- original_rep %>% factor() %>% `levels<-`(value = seq_along(levels(.)))
H3K27_Somatic <- tar_read(chic.experiment.granges_H3K27_Somatic_CN_chr) %>%
  split(seqnames(.)) %>%
  sapply(\(gr) ksmooth_sliding_windows(gr, bw = 40)) %>%
  GRangesList() %>%
  unlist() %>%
  `elementMetadata<-`(
    value = elementMetadata(.) %>%
      apply(2, \(v) v %>% round(digits = 3)) %>%
      as.data.frame() %>%
      `colnames<-`(
        value = colnames(.) %>%
          mapply(
            \(n, original, new) n %>% str_replace(original, new),
            .,
            paste0("Rep", as.character(original_rep)),
            paste0("Rep", as.character(new_rep))
          )
      )
  )
save_bed("H3K27_Somatic", H3K27_Somatic, "geoSubmission/H3K27_Somatic.tsv")
original_rep <- chic.samples %>% subset(driver == "tj" & group == "H3K9") %>% pull(rep)
new_rep <- original_rep %>% factor() %>% `levels<-`(value = seq_along(levels(.)))
H3K9_Somatic <- tar_read(chic.experiment.granges_H3K9_Somatic_CN_chr) %>%
  split(seqnames(.)) %>%
  sapply(\(gr) ksmooth_sliding_windows(gr, bw = 40)) %>%
  GRangesList() %>%
  unlist() %>%
  `elementMetadata<-`(
    value = elementMetadata(.) %>%
      apply(2, \(v) v %>% round(digits = 3)) %>%
      as.data.frame() %>%
      `colnames<-`(
        value = colnames(.) %>%
          mapply(
            \(n, original, new) n %>% str_replace(original, new),
            .,
            paste0("Rep", as.character(original_rep)),
            paste0("Rep", as.character(new_rep))
          )
      )
  )
save_bed("H3K9_Somatic", H3K9_Somatic, "geoSubmission/H3K9_Somatic.tsv")

X <- chic.samples %>%
  subset(str_starts(molecule, "H3")) %>%
  mutate(celltype = c(nos="Germline", tj="Somatic")[driver]) %>%
  mutate(group = group %>% factor(c("H3K4", "H3K27", "H3K9"))) %>%
  group_by(celltype, group) %>%
  summarise(
    sample_list = list(sample_glob[order(molecule, as.factor(rep))]),
    rep_list = list(levels(as.factor(rep))),
    molecule = list(sort(molecule)),
    rep_print = list(rep(seq_along(rep_list[[1]]), 2))
  )

lst <- list(
  H3K4_Germline,
  H3K27_Germline,
  H3K9_Germline,
  H3K4_Somatic,
  H3K27_Somatic,
  H3K9_Somatic
)
for (i in seq(nrow(X))) {
  cat("## ")
  cat(paste(X$group[i], X$celltype[i], sep="_"))
  cat("\nSample names (R1 and R2):")
  for (j in seq(length(X$sample_list[[i]]))) {
    s <- X$sample_list[[i]][[j]]
    cat("\n")
    cat("\n")
    cat(X$molecule[[i]][[j]])
    cat("_Rep")
    cat(X$rep_print[[i]][[j]])
    cat(": ")
    if (grepl("\\?", s)) {
      cat(s %>% str_replace("\\?", "2"))
      cat(", ")
      cat(s %>% str_replace("\\?", "3"))
    } else {
      cat(s)
    }
  }
  cat("\n")
  cat("\n")
}
