# Wrapper for regression for our ChIC experiment ----
# Wraps glmGamPoi for our bulk HTS regression. We generally set overdispersion
# to "global" as this is fast and we don't see evidence of different dispersion
# in different genomic windows. Our two-tailed test option that we wrap is
# test_de (likelihood ratio test). Our one-tailed test option will be a Wald
# test using predict.glmGamPoi.
chic_quantify <- function(
    chic_experiment,
    molecule,
    replicate_nums,
    offset,
    global_overdispersion = TRUE,
    test_de = FALSE,
    test_wald = FALSE) {
  if (length(replicate_nums) <= 2) {
    return(
      GRanges(
        seqnames(chic_experiment),
        ranges(chic_experiment)
      ) %>%
        `elementMetadata<-`(
          value = (
            as.matrix(elementMetadata(chic_experiment))
            / exp(offset)
          ) %>%
            replace(. < 1e-8, 0) %>%
            as.data.frame() %>%
            cbind(
              p_peak = 1, L2FC = 0
            )
        )
    )
  }
  # We put "molecule" and "rep" in our tar_map, so create new names. As
  # "molecule" coefs come first, our p_peak contrast will be the log fold
  # change between the two levels of "molecule".
  if (length(unique(replicate_nums)) > 1) {
    R <- replicate_nums %>%
      factor() %>%
      `contrasts<-`(
        value = contr.helmert(length(levels(.)))
      )
  } else {
    R <- NA
  }
  if (length(unique(replicate_nums)) > 1) {
    fmla <- ~ 0 + mol + R
  } else {
    fmla <- ~ 0 + mol
  }
  fit <- glm_gp(
    as.matrix(chic_experiment@elementMetadata),
    fmla,
    tibble(
      mol = molecule,
      R
    ),
    size_factors = 1,
    offset = offset,
    overdispersion = list(TRUE, "global")[
      c(isFALSE(global_overdispersion), isTRUE(global_overdispersion))
    ][[1]],
    overdispersion_shrinkage = list(TRUE, FALSE)[
      c(isFALSE(global_overdispersion), isTRUE(global_overdispersion))
    ][[1]],
    verbose = TRUE
  )
  hypothesis_testing <- tibble(
    if (test_de) {
      test_de(
        fit,
        c(-1, 1) %>% c(rep(0, ncol(fit$Beta) - length(.))),
        verbose = T
      )
    } else {
      tibble(pval = 1, lfc = 0)
    },
    if (test_wald) {
      with(
        predict(
          fit,
          cbind(
            diag(c(1, 1)),
            matrix(0, nrow = 2, ncol = ncol(fit$Beta) - 2)
          ),
          offset = 0,
          se.fit = T,
          verbose = TRUE
        ),
        as_tibble(
          cbind(
            matrix(
              fit,
              nrow = nrow(fit),
              ncol = ncol(fit),
              dimnames = list(
                NULL, paste0("logMu_", levels(as.factor(molecule)))
              )
            ),
            matrix(
              se.fit,
              nrow = nrow(se.fit),
              ncol = ncol(se.fit),
              dimnames = list(
                NULL, paste0("se_", levels(as.factor(molecule)))
              )
            )
          )
        )
      )
    } else {
      as_tibble(
        cbind(
          matrix(
            0,
            nrow = nrow(fit$Beta),
            ncol = 2,
            dimnames = list(NULL, paste0("logMu_", levels(as.factor(molecule))))
          ),
          matrix(
            0,
            nrow = nrow(fit$Beta),
            ncol = 2,
            dimnames = list(NULL, paste0("se_", levels(as.factor(molecule))))
          )
        )
      )
    }
  )
  GRanges(
    seqnames(chic_experiment),
    ranges(chic_experiment),
    seqlengths = seqlengths(chic_experiment)
  ) %>%
    `elementMetadata<-`(
      value = with(
        hypothesis_testing,
        tibble(
          as.data.frame(exp(fit$Beta)) %>%
            replace(. < 1e-8, 0) %>%
            rename_with(~ str_glue("score.{.}")),
          p_peak = pval,
          L2FC = lfc,
          subset(hypothesis_testing, select = grep("logMu_|se_", colnames(hypothesis_testing)))
        )
      )
    ) %>%
    `metadata<-`(
      value = list(
        overdispersions = fit$overdispersions,
        overdispersion_shrinkage_list = fit$overdispersion_shrinkage_list,
        deviances = fit$deviances,
        ridge_penalty = fit$ridge_penalty,
        model_matrix = fit$model_matrix,
        design_formula = fit$design_formula
      )
    )
}

# One-tailed p-value calling. We ended up with tracks with two levels of
# granularity: window width 500bp or 40bp. For each window being tested, we can
# paint the entire width covered by the window as enriched, then reduce these
# windows at the end.
reduce_peaks_2_tracks <- function(
    broad_track,
    rough_track,
    broad_track_window_size = 500,
    rough_track_window_size = 40,
    p = 0.001) {
  chrs <- GRanges(
    seqlevels(broad_track),
    IRanges(1, width = seqlengths(broad_track))
  )
  bpeak <- broad_track[
    broad_track$p_peak < p & broad_track$L2FC > 0
  ] %>%
    resize(broad_track_window_size, fix = "center") %>%
    GenomicRanges::reduce()
  rpeak <- rough_track[
    rough_track$p_peak < p & rough_track$L2FC > 0
  ] %>%
    resize(rough_track_window_size, fix = "center") %>%
    GenomicRanges::reduce()
  rpeak <- rpeak[width(rpeak) >= 1.5 * rough_track_window_size]
  all_peaks <- GRangesList(list(bpeak, rpeak)) %>%
    unlist() %>%
    GenomicRanges::reduce()
  all_peaks %>%
    restrict(
      start = 1L, end = seqlengths(broad_track)[as.character(seqnames(all_peaks))]
    )
}

# Converts object GRanges to coverage GRanges.
# The objects must already be reduced. We will create a new "score" which is 1
# for these objects, and 0 for runs of 0 coverage across the genome. The new
# track is suitable for export to BigWig.
peaks_to_genome_coverage <- function(object_granges) {
  object_granges %>%
    split(seqnames(object_granges)) %>%
    as.list() %>%
    enframe() %>%
    rowwise() %>%
    reframe(
      value = coverage(value)[[name]] %>%
        attributes() %>%
        with(
          if (length(lengths) == 0) {
            GRanges()
          } else {
            GRanges(
              name,
              IRanges(
                cumsum(c(1, lengths[-length(lengths)])),
                width = lengths
              ),
              score = as.numeric(values)
            )
          }
        ) %>%
        list()
    ) %>%
    unlist(use.names = F) %>%
    GRangesList() %>%
    unlist(use.names = F) %>%
    GRanges(seqinfo = seqinfo(object_granges))
}

# Cleaning up visualization of stat. signif nucleosomes. ----
nucleosomes_cleanup_tracks <- function(granges, p_nuc, p_enriched, nuc_size = 147, step_size = 20) {
  p_enriched <- p_enriched %>%
    replace(which(p_nuc >= 0.05), 1) %>%
    replace(is.na(p_nuc), 1)
  peaks <- p_nuc %>%
    split(seqnames(granges)) %>%
    mapply(
      \(n, p) as(
        (p < 0.05) %>% replace_na(FALSE), "Rle"
      ) %>%
        attributes() %>%
        with(
          GRanges(
            n,
            IRanges(
              cumsum(c(1, lengths[-length(lengths)])),
              width = lengths
            )
          )[values]
        ),
      seqlevels(granges),
      .
    ) %>%
    GRangesList()
  peaks <- peaks %>%
    sapply(
      \(gr) gr[width(gr) <= floor(nuc_size / step_size)] %>%
        resize(floor(nuc_size / step_size), fix = "center") %>%
        GenomicRanges::reduce()
    ) %>%
    GRangesList()
  diff_enriched_logical <- mapply(
    \(pk, pdiff) mapply(
      \(start, end) min(
        pdiff[pmax(1, start):pmin(length(pdiff), end)],
        na.rm = T
      ) < 0.05,
      start(pk),
      end(pk)
    ) %>%
      as.logical(),
    peaks,
    p_enriched %>% split(seqnames(granges)),
    SIMPLIFY = FALSE
  ) %>%
    do.call(c, .)
  nucs <- mapply(
    \(gr, pk) GRanges(
      seqnames(pk),
      ranges = IRanges(
        start = start(gr)[pmax(1, as.numeric(start(pk)))],
        end = end(gr)[pmin(length(gr), as.numeric(end(pk)))]
      ),
      seqlengths = seqlengths(pk)
    ),
    granges %>% split(seqnames(granges)),
    peaks
  ) %>%
    GRangesList() %>%
    unlist(use.names = F) %>%
    resize(pmax(nuc_size, width(.)), fix = "center")
  diffs <- nucs[diff_enriched_logical]
  list(
    Nucleosomes = nucs[width(nucs) <= 2 * nuc_size] %>%
      GenomicRanges::reduce() %>%
      resize(nuc_size, fix = "center"),
    Diff_Enriched = diffs[width(diffs) <= 2 * nuc_size] %>%
      GenomicRanges::reduce() %>%
      resize(nuc_size, fix = "center")
  )
}

# Genomic ranges with width 500 at the head of the CDS (TSS). ----
gtf_granges_extended_from_tss <- function(granges, genes) {
  seqlevels(granges) <- seqlevels(granges) %>% c("*")
  overlaps <- GRanges(
    genes$chr %>% replace(is.na(.), "*"),
    IRanges(
      ifelse(genes$strand == "+", genes$start, genes$end - 499) %>%
        replace(is.na(.), -1),
      width = 500
    )
  ) %>%
    findOverlaps(granges) %>%
    as.list()
  mapply(
    \(overlaps, start, end, strand) if (length(overlaps) == 0) {
      0
    } else if (strand == "+") {
      min(end(granges)[overlaps], end) - start + 1
    } else {
      end - max(start(granges)[overlaps], start) + 1
    },
    overlaps,
    genes$start,
    genes$end,
    genes$strand
  ) %>%
    setNames(pull(genes, 1))
}

# Violin plots of histone code L2FC. ----
# Plot the combined chr feature violins once. Cannot use the same
# position = dodge on a boxplot as we used on the violin, so this is going to be
# a plainer plot.
plot_enriched_chromosomes_combined <- function(data, label_pattern, title = label_pattern) {
  data <- data %>% subset(grepl(label_pattern, label))
  ggplot(
    data,
    aes(label, L2FC, fill = celltype, group = interaction(label, celltype, mark))
  ) +
    geom_violin() +
    # geom_boxplot(outlier.shape = NA, fill = "transparent") +
    annotate(
      "text",
      (1 + seq(3 * length(unique(data$label))))/3, 1.1,
      label = rep(c("H3K4me3", "H3K27me3", "H3K9me3"), length(unique(data$label))),
      size = 2
    ) +
    scale_fill_manual(values = cell_type_violin_colors) +
    coord_cartesian(NULL, c(-1.2, 1.2)) +
    labs(title = title) +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank()
    )
}

# Use a built plot_enriched_chromosomes_combined and add additional information
# to the plot.
annotate_enriched_chromosomes_combined <- function(gg) {
  data <- ggplot_build(gg)$data[[1]]
  # 75% CI information. We persist the density results in ggplot_build data, as
  # well as the final coordinate space x and y, so we need to match the quantile
  # to plot to the cumsum of % density and then get the y value for this
  # particular variable value.
  box_data <- data %>%
    group_by(group) %>%
    summarise(
      xmin = xmin[1] + (xmax - xmin)[1] * (1 - width[1]) / 2,
      xmax = xmax[1] - (xmax - xmin)[1] * (1 - width[1]) / 2,
      ymin = y[
        findInterval(
          0.25,
          cumsum(density) / sum(density)
        )
      ],
      ymax = y[
        findInterval(
          0.75,
          cumsum(density) / sum(density)
        )
      ]
    )
  median_data <- data %>%
    group_by(group) %>%
    summarise(
      x = xmin[1] + (xmax - xmin)[1] * (1 - width[1]) / 2,
      x = xmin[1],
      y = y[
        findInterval(
          0.5,
          cumsum(density) / sum(density)
        )
      ],
      xend = xmax[1] - (xmax - xmin)[1] * (1 - width[1]) / 2,
      yend = y
    )
  # We will set the breaks and labels of the continuous x axis to match the
  # discrete x axis which we built.
  labels <- ggplot_build(gg)$layout$panel_params[[1]]$x$limits
  ggplot(data, aes()) +
    geom_violin(
      aes(x, y, fill = fill, xmin = xmin, xmax = xmax, violinwidth = violinwidth, group = group),
      stat="identity"
    ) +
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      box_data,
      color = "black",
      fill = "transparent"
    ) +
    geom_segment(
      aes(x, y, xend = xend, yend = yend),
      median_data,
      color = "black",
      linewidth = 1.25
    ) +
    annotate(
      "text",
      rep(c(-0.3, 0, 0.3), length(labels)) + rep(seq_along(labels), each = 3),
      1.5,
      label = rep(c("H3K4me3", "H3K27me3", "H3K9me3"), length(labels)),
      size = 2
    ) +
    scale_x_continuous(breaks = seq_along(labels), labels = labels) +
    scale_fill_identity() +
    coord_cartesian(NULL, c(-1.9, 1.9)) +
    labs(x = "label", y = "L2FC")
}

# Counts: Grouped by chromosome_arms_diameter_40 for each sample. The
# diameter-40 windows have exactly 50% overlap between successive windows; the
# fragment should be counted in exactly 2 of the windows that we are summing up,
# but the total count is not even because of the window at the start or end of
# the chromosome_arms group.
half <- \(x) x/2
enriched_chromosomes_build_summarized_experiment <- function(
  counts,
  chromosome_arms_diameter_40,
  model_frame,
  offset_germline,
  offset_somatic
) {
  counts = counts %>%
    split(chromosome_arms_diameter_40$group) %>%
    sapply(\(df) df %>% as.matrix %>% colSums) %>%
    t() %>%
    matrix(
      nrow = nrow(.),
      dimnames = list(
        rownames(.),
        model_frame$rowname
      )
    ) %>%
    half() %>%
    round()
  celltype. <- model_frame$celltype %>% factor(c("GSC", "CySC"))
  rep. <- interaction(model_frame$rep, model_frame$celltype) %>%
    droplevels()
  num_gsc_batches <- length(unique(rep.[celltype. == "GSC"]))
  num_cysc_batches <- length(unique(rep.[celltype. == "CySC"]))
  contrasts(rep.) <- cbind(
    Celltype = rep(
      c(1, -1),
      c(num_gsc_batches, num_cysc_batches)
    ),
    bdiag(
      contr.sum(num_gsc_batches), contr.sum(num_cysc_batches)
    ) %>%
      as.matrix()
  )
  sf <- exp(
    as.numeric(
      sapply(
        list(
          elementMetadata(offset_germline)[1,],
          elementMetadata(offset_somatic)[1,]
        ),
        unlist
      )
    )
  )
  SummarizedExperiment(
    counts,
    colData = mutate(
      model_frame,
      celltype = celltype.,
      rep = rep.,
      sf = sf,
    )
  )
}

# Aggregate significance test for the chromosomes compare model. This produces
# significance levels of q-value for 33 tests (chromosome arm genomic ranges &
# 3 marks being tested).
enriched_chromosomes_compare_signif <- function(data) {
  pvalues <- rbind(
    H3K4 = pull(data$H3K4, pval, name),
    H3K27 = pull(data$H3K27, pval, name),
    H3K9 = pull(data$H3K9, pval, name)
  )
  pvalues <- pvalues %>%
    replace(
      seq_along(.),
      p.adjust(as.numeric(.), "BH")
    )
  signif <- pvalues %>%
    as.numeric %>%
    cut(c(0, 1e-4, 1e-3, 1e-2, 0.05)) %>%
    `levels<-`(value = c("****", "***", "**", "*"))
  matrix(
    signif,
    nrow = nrow(pvalues),
    ncol = ncol(pvalues),
    dimnames = dimnames(pvalues)
  )
}

# Masked Transposable Elements Analysis ----
factor_transposable_elements <- function(
  te_ids,
  transposon_sequence_set_metadata
) {
  stopifnot(all(te_ids %in% transposon_sequence_set_metadata$V1))
  transposon_sequence_set_metadata <- transposon_sequence_set_metadata[
    match(te_ids, transposon_sequence_set_metadata$V1),
  ]
  specific_type <- (
    transposon_sequence_set_metadata$V15 %>%
      sapply(
        \(n) n %>%
          strsplit(" <newline> ") %>%
          unlist(use.names = F) %>%
          setdiff("transposable_element ; SO:0000101") %>%
          `[`(value = 1)
      )
  )
  ifelse(
    grepl("retrotransposon", specific_type),
    "retrotransposon",
    ifelse(
      grepl("foldback_element|terminal_inverted_repeat_element", specific_type),
      "DNA transposon",
      NA
    )
  ) %>%
    factor() %>%
    setNames(te_ids)
}

plot_transposable_element_enrich_chromatin <- function(
  H3K4_Germline, H3K4_Somatic, H3K27_Germline, H3K27_Somatic, H3K9_Germline, H3K9_Somatic
) {
  plot_features <- which(
    H3K4_Germline$score.molH3 >= 0.5 &
      H3K27_Germline$score.molH3 >= 0.5 &
      H3K9_Germline$score.molH3 >= 0.5 &
      H3K4_Somatic$score.molH3 >= 0.5 &
      H3K27_Somatic$score.molH3 >= 0.5 &
      H3K9_Somatic$score.molH3 >= 0.5
  )
  data <- tribble(
    ~mark, ~celltype, ~track,
    "H3K4", "Germline", H3K4_Germline,
    "H3K27", "Germline", H3K27_Germline,
    "H3K9", "Germline", H3K9_Germline,
    "H3K4", "Somatic", H3K4_Somatic,
    "H3K27", "Somatic", H3K27_Somatic,
    "H3K9", "Somatic", H3K9_Somatic,
  ) %>%
    rowwise() %>%
    reframe(
      mark = mark %>% factor(chic.mark.data$mark),
      celltype,
      value = track$L2FC[plot_features] %>% subset(between(., -5, 5))
    )
  data_dodge <- data %>%
    rowwise() %>%
    mutate(
      x = list(
        H3K4 = c(Germline=0.775, Somatic=1.225),
        H3K27 = c(Germline=1.775, Somatic=2.225),
        H3K9 = c(Germline=2.775, Somatic=3.225)
      )[[mark]][celltype] +
        runif(1, min = -0.175, max = 0.175)
    )
  (
    data %>%
      ggplot(aes(mark, value, fill=celltype)) +
      geom_violin(
        scale = "width",
        adjust = 0.5,
        width = 0.8,
        position = position_dodge(width = 0.9)
      ) +
      rasterise(
        geom_point(
          aes(x=x, fill=NULL),
          data = data_dodge,
          stroke = NA,
          size = 0.6
        ),
        dpi = 300
      ) +
      annotate(
        "rect",
        xmin = -0.02 + c(0.575, 1.025) + rep(0:2, each=2),
        xmax = 0.02 + c(0.975, 1.425) + rep(0:2, each=2),
        ymin = summarise(group_by(data, mark, celltype), y = quantile(value, 0.25))$y,
        ymax = summarise(group_by(data, mark, celltype), y = quantile(value, 0.75))$y,
        color = "black",
        fill = "transparent"
      ) +
      annotate(
        "segment",
        x = -0.02 + c(0.575, 1.025) + rep(0:2, each=2),
        xend = 0.02 + c(0.975, 1.425) + rep(0:2, each=2),
        y = summarise(group_by(data, mark, celltype), y = median(value))$y,
        yend = summarise(group_by(data, mark, celltype), y = median(value))$y,
        linewidth = 1.5
      ) +
      scale_fill_manual(values = cell_type_violin_colors) +
      coord_cartesian(c(0.5, 3.5), c(-1.55, 1.55), expand = FALSE) +
      labs(
        title = "Transposable Elements",
        subtitle = str_glue("Presence Filter: H3 FPKM >= 0.5 * H3 Autosome Median. n = {length(plot_features)}"),
        y = "L2FC"
      ) +
      theme(
        aspect.ratio = 3/2,
        legend.position = "none"
      )
  )
}
