library(dplyr)
library(forcats)
library(ggnewscale)
library(ggplot2)
library(magrittr)
library(reshape2)
library(rtracklayer)
library(scales)
library(stringr)
library(targets)
library(viridisLite)

lift <- import.chain("illustrations/dm3ToDm6.over.chain")
repli <- tar_read(repli.timing_Kc167_chr)
source('R/1-constants.R')
source('R/chic-heatmap.R')
source('R/repli-chic-heatmap.R')
source('R/repli-logistic.R')

lift_dm3_to_flybase_dm6 <- \(gr)
  gr %>%
    liftOver(lift) %>%
    unlist() %>%
    attributes() %$%
    GRanges(
      as.factor(seqnames) %>%
        recode(
          chr2L="2L",
          chr2R="2R",
          chr3L="3L",
          chr3R="3R",
          chr4="4",
          chrX="X",
          chrY="Y",
          chrM="mitochondrion_genome"
        ),
      ranges,
      score = elementMetadata$score
    )
time_dm6 <- \(gr) {
  hits <- findOverlaps(
    GRanges(seqnames(gr), IRanges(mid(ranges(gr)), width=1)),
    repli
  )
  gr$timing <- repli$score[
    sapply(hits, \(v) v[1])
  ]
  gr
}
H3K4 <- read.table("illustrations/GSE45088_repset.17402831.smoothedM.bedgraph.gz", sep=" ", skip=1) %$%
  GRanges(V1, IRanges(V2 + 1, V3), score = V4) %>%
  lift_dm3_to_flybase_dm6() %>%
  time_dm6()
dm6_reference <- H3K4
H3K27 <- read.table("illustrations/GSE45083_H3K27.bedgraph.gz", sep=" ", skip=1) %$%
  GRanges(V1, IRanges(V2 + 1, V3), score = V4) %>%
  lift_dm3_to_flybase_dm6() %>%
  time_dm6()
H3K9 <- read.table("illustrations/GSE27796_H3K9.bedgraph.gz", sep=" ", skip=1) %$%
  GRanges(V1, IRanges(V2 + 1, V3), score = V4) %>%
  lift_dm3_to_flybase_dm6() %>%
  time_dm6()

K4bar <- tibble(
  timing = cut(H3K4$timing, seq(-0.75, 0.75, by=0.002)),
  L2FC = H3K4$score
) %>%
  subset(!is.na(timing) & as.numeric(timing) != 10) %>%
  group_by(timing) %>%
  summarise(L2FC = mean(L2FC))
 proj <- tibble(
  H3K4 = H3K4$score,
  H3K27 = H3K27$score,
  H3K9 = H3K9$score,
  bp = width(dm6_reference)
) %>%
  split(cut(H3K4$timing, seq(-1, 1, by=0.002))) %>%
  sapply(\(df) tibble(
    H3K4 = mean(df$H3K4), H3K27 = mean(df$H3K27), H3K9 = mean(df$H3K9),
    sample_size_bp = sum(df$bp)
  ), simplify=FALSE) %>%
  bind_rows(.id = "timing") %$%
  tibble(
    repli = seq(-1, 0.998, by=0.002),
    H3K4 = exp(H3K4 * log(2)),
    H3K27 = exp(H3K27 * log(2)),
    H3K9 = exp(H3K9 * log(2)),
    sample_size_bp
  ) %>%
  as.matrix()
proj <- proj %>%
  cbind(chr2L = 1, chr2LC = 1, chr2RC = 1, chr2R = 1, chr3L = 1, chr3LC = 1, chr3RC = 1, chr3R = 1, chr4 = 1, chrX = 1, chrY = 1, chrscaffolds = 1)
names(dimnames(proj)) <- c("timing", "series")
cairo_pdf("figure/Repli-CHIP-chip-Kc167.pdf", w = 6, h = 4.5)
print(plot_repli_track_raster(proj, log2_limits = NULL))
dev.off()
