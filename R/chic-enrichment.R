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

write_chic_peaks <- function(peak_table_list, output_path) {
  dir.create(dirname(output_path), recursive = TRUE, showW = FALSE)
  peaks_bed <- peak_table_list %>%
    sapply(
      \(tab) tab %>% subset(q < 0.1),
      simplify = F
    ) %>%
    bind_rows(.id = "mark") %>%
    mutate(
      mark = paste0(
        mark,
        " (",
        cut(
          q,
          c(0, 1e-4, 1e-3, 1e-2, 5e-2, 1e-1, Inf)
        ) %>%
          fct_recode(
            `****` = "(0,0.0001]",
            `***` = "(0.0001,0.001]",
            `**` = "(0.001,0.01]",
            `*` = "(0.01,0.05]",
            `~` = "(0.05,0.1]"
          ),
        ")"
      )
    ) %>%
    subset(select = c(chr, start, end, mark)) %>%
    arrange(chr, start, mark)
  with_options(
    list(scipen = 100),
    peaks_bed %>%
      write.table(output_path, sep = "\t", quote = F, row.names = F, col.names = F)
  )
  output_path
}

chic_track_generate_table_by_enrichment <- function(chic_df_list, enrichment_threshold = 1.5) {
  enrichment_track <- chic_df_list %>%
    sapply(\(df) df %>% with(Rle(enrichment, length))) %>%
    RleList()
  grouped_df <- (enrichment_track >= enrichment_threshold) %>%
    mapply(
      \(name, rle, p_values) tibble(
        start = 1 + c(0, cumsum(rle@lengths[-length(rle@lengths)])),
        length = rle@lengths,
        end = start + length - 1,
        keep_run = rle@values,
        p = p_values %>%
          split(
            Rle(
              seq_along(rle@lengths),
              rle@lengths
            )
          ) %>%
          list(p = .) %>%
          with(
            ifelse(
              keep_run,
              sapply(p, min),
              sapply(p, mean)
            )
          )
        # sapply(min)
      ),
      names(chic_df_list),
      .,
      chic_df_list %>% chic_df_list_to_pvalue_list(),
      SIMPLIFY = FALSE
    ) %>%
    bind_rows(.id = "chr") %>%
    ungroup() %>%
    mutate(
      q = p.adjust(p, "BH")
    ) %>%
    subset(keep_run) %>%
    `[`(c("chr", "start", "end", "p", "q"))
}

chic_df_list_to_pvalue_list <- function(df_list) {
  df_list %>%
    sapply(\(df) df %>% with(Rle(p, length))) %>%
    RleList()
}

display_peak_location_stats <- function(lst) {
  lst %>%
    sapply(
      \(df) df$region %>%
        table() %>%
        enframe(),
      simplify = FALSE
    ) %>%
    bind_rows(.id = "mark") %>%
    mutate(
      mark = mark %>% factor(chic.mark.data$mark),
      name = name %>% factor(lst[[1]]$region %>% levels())
    ) %>%
    ggplot(
      aes(x = name, y = value, fill = mark)
    ) +
    geom_bar(
      stat = "identity", position = "dodge"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_fill_viridis_d(
      option = "turbo", begin = 0.3, end = 0.9
    ) +
    labs(
      x = "Genomic Annotation",
      y = "Number of Peaks (q < 0.05)"
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

# Plot the combined chr feature violins once. Cannot use the same
# position = dodge on a boxplot as we used on the violin, so this is going to be
# a plainer plot.
plot_enriched_chromosomes_combined <- function(data, label_pattern, title = label_pattern) {
  fill_colors <- c(Germline="#96C586", Somatic="#A97AAC")
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
    scale_fill_manual(values = fill_colors) +
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