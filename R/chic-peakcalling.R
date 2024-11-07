# Broad Peaks analysis ----
calculate_chic_gene_enrichment_broad_logit_trend_fit <- function(
  assay.data.sc,
  Upd_cpm
) {
  logCDS <- with(assay.data.sc, log(abs(start - end) + 1) / log(10))
  reference_level <- mean(
    c(
      mean(subset(logCDS, Upd_cpm[, "germline"] >= 5 & !is.na(logCDS))),
      mean(subset(logCDS, Upd_cpm[, "somatic"] >= 5 & !is.na(logCDS)))
    )
  )
  lower_bound_fit <- c(
    yintercept=2 * log(10),
    slope=-log(10)
  )
  # The fitted line that we will subtract off. To use a double-negative, when
  # we subtract off the effect of logCDS, we don't want to standardize the
  # values as to a corrected value where logCDS == 0, but instead to the level
  # of logCDS == reference_level. So find the y-intercept of that line, which
  # will be greater than the starting yintercept because slope is negative.
  lower_bound_fit %>%
    replace(
      1,
      -reference_level * lower_bound_fit["slope"]
    )
}

plot_chic_gene_enrichment_broad_logit_trend_fit <- function(
  chic.gene.enrichment.broad.rawvalue,
  Upd_cpm,
  assay.data.sc
) {
  abline_par <- c(
    intercept=2 * log(10),
    slope=-log(10)
  )
  cds <- abs(assay.data.sc$start - assay.data.sc$end) + 1
  points <- bind_rows(
    list(
      Germline=data.frame(
        pct=c(
          chic.gene.enrichment.broad.rawvalue$H3K4_Germline,
          chic.gene.enrichment.broad.rawvalue$H3K27_Germline
        ),
        cpm=Upd_cpm[, "germline"] %>% setNames(NULL),
        cds
      ),
      Somatic=data.frame(
        pct=c(
          chic.gene.enrichment.broad.rawvalue$H3K4_Somatic,
          chic.gene.enrichment.broad.rawvalue$H3K27_Somatic
        ),
        cpm=Upd_cpm[, "somatic"] %>% setNames(NULL),
        cds
      )
    ),
    .id = "celltype"
  ) %>%
    subset(
      between(pct, 0.00001, 0.99999) &
        cpm >= 5
    ) %>%
    reframe(
      celltype,
      logit = qlogis(pct),
      logCPM = log(cpm) / log(10),
      logCDS = log(cds) / log(10)
    ) %>%
    as_tibble()
  points[sample(nrow(points)), ] %>%
    ggplot(aes(logCDS, logit)) +
    geom_abline(
      intercept = abline_par["intercept"], slope = abline_par["slope"],
      color = "#990000"
    ) +
    rasterise(geom_point(aes(color=celltype), stroke=NA, size=0.5), dpi=300) +
    scale_color_manual(values = unlist(chic_line_track_colors, use.names = FALSE)) +
    labs(
      x = bquote(log[10]*"(Gene Length)"),
      y = "logit(% Gene Body Enriched)"
    ) +
    theme(
      aspect.ratio = 3/4,
      legend.position = "none"
    )
}

calculate_chic_gene_enrichment_broad_correction <- function(
  chic.gene.enrichment.broad,
  chic.gene.enrichment.broad.logit.trend.fit
) {
  logCDS <- with(
    chic.gene.enrichment.broad,
    log(abs(start - end) + 1) / log(10)
  )
  correction <- -chic.gene.enrichment.broad.logit.trend.fit[1] -
    chic.gene.enrichment.broad.logit.trend.fit[2] * logCDS
  chic.gene.enrichment.broad[-(1:6)] <- chic.gene.enrichment.broad[-(1:6)] %>%
    ungroup() %>%
    summarise_all(
      \(v) v %>%
        replace(
          which(v != 0 & v != 1),
          v %>%
            qlogis() %>%
            `+`(correction) %>%
            plogis() %>%
            subset(v != 0 & v != 1)
        )
    )
  chic.gene.enrichment.broad
}
