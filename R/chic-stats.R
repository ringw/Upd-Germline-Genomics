chic_stats_boxplot <- function(track) {
  jitter_sd_width <- 0.12
  jitter_height <- 0.001
  tibble(
    x = droplevels(as.factor(seqnames(track))) %>%
      fct_recode(`2`="2L", `2`="2R", `3`="3L", `3`="3R") %>%
      relevel("X"),
    y = track$score
  ) %>%
    ggplot(aes(x, y)) +
    geom_boxplot(color="#3d3c3b", fill="#f0f1a1", outlier.shape = NA) +
    rasterise(
      geom_point(
        data=\(data) data %>%
          group_by(x) %>%
          slice_sample(n=5000),
        shape=".",
        position=ggforce::position_jitternormal(jitter_sd_width, jitter_height)
      ),
      dpi=200
    ) +
    scale_y_continuous(limits=c(0, 5)) +
    coord_cartesian(
      c(0.5 - 0.05, 5.5 + 0.05),
      c(0, 3),
      expand=FALSE
    ) +
    theme_bw() +
    labs(x = NULL, y = "Enrichment (vs Auto Monosome Median)")
}

histogram_paired_end_fragment_size <- function(bulk_reads) {
  fragment_sizes <- do.call(
    c,
    sapply(
      bulk_reads,
      \(df) df %>% paired_end_reads_to_fragment_lengths %>% pull("length"),
      simplify=FALSE
    )
  )
  tbl <- as.data.frame(table(fragment_sizes))
  tbl <- deframe(tbl)[
    as.character(seq(50, 500))
  ] %>%
    enframe("fragment_sizes", "Freq") %>%
    mutate(fragment_sizes = as.integer(fragment_sizes))
  tbl %>%
    ggplot(aes(fragment_sizes, y=0, height=2*Freq)) +
    geom_tile(linewidth=0.01, color="black", fill="black") +
    scale_y_continuous(labels = scales::unit_format(unit = "K", scale = 1e-3)) +
    coord_cartesian(NULL, c(0, max(tbl$Freq) * 1.05), expand=FALSE) +
    labs(x = "Fragment Size (bp)", y = "Num H3 Fragments") +
    theme_bw() +
    theme(aspect.ratio = 1/2)
}

histogram_paired_end_fragment_size <- function(bulk_reads, faceted = FALSE) {
  lvl_lookup <- c(
    X="X", `2L`="2", `2R`="2", `3L`="3", `3R`="3", `4`="4", Y="Y"
  )
  lvl_dfs <- sapply(
    bulk_reads,
    \(df) lvl_lookup[as.character(df$rname[1])]
  )
  fragment_sizes <- sapply(
    bulk_reads,
    \(df) df %>% paired_end_reads_to_fragment_lengths %>% pull("length"),
    simplify=FALSE
  )
  fragment_sizes <- do.call(
    rbind,
    mapply(
      tibble,
      chr = lvl_dfs,
      size = fragment_sizes,
      SIMPLIFY=FALSE,
      USE.NAMES=FALSE
    )
  )
  fragment_sizes$chr <- fragment_sizes$chr %>% factor(c("X", "2", "3", "4", "Y"))
  if (faceted)
    tbl <- as.data.frame(with(fragment_sizes, table(chr, size)))
  else
    tbl <- as.data.frame(table(fragment_sizes))
  tbl$size <- tbl$size %>% as.character %>% as.integer
  tbl <- tbl %>% subset(between(size, 50, 500))
  # tbl <- deframe(tbl)[
  #   as.character(seq(50, 500))
  # ] %>%
  #   enframe("fragment_sizes", "Freq") %>%
  #   mutate(fragment_sizes = as.integer(fragment_sizes))
  tbl %>%
    ggplot(aes(size, y=0, height=2*Freq)) +
    geom_tile(linewidth=0.01, color="black", fill="black") +
    facet_wrap(if (faceted) vars(chr) else vars(), scales="free", ncol=1) +
    scale_y_continuous(
      labels = scales::unit_format(unit = "K", scale = 1e-3),
      # As we want expand = FALSE, we will apply the y-expansion in a callable.
      # We only want to expand the + y margin, not the - y margin.
      limits = \(lm) lm + c(0, 1) * 0.05 * diff(lm)
    ) +
    coord_cartesian(NULL, c(0, NA), expand=F) +
    labs(x = "Fragment Size (bp)", y = "Num H3 Fragments") +
    theme_bw() +
    theme(aspect.ratio = 1/2, plot.margin = margin(5.5, 15, 5.5, 5.5))
}