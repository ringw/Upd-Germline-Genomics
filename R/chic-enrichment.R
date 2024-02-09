subsample_time_series <- function(series, bin_size=100) {
  sampling <- seq(floor(bin_size / 2), length(series), by=bin_size)
  # If the final window is less than half-full, then we need to sample an
  # extra value at the final window.
  if (between(length(series) %% bin_size, 1, floor(bin_size/2) - 1))
    sampling <- c(sampling, length(series) - (length(series) %% bin_size) + 1)
  setNames(as.numeric(series[sampling]), sampling)
}

subsample_unlist_rle_list <- function(rles, bin_size) do.call(
  c,
  sapply(
    rles,
    \(v) subsample_time_series(v, bin_size),
    simplify = FALSE
  )
)

chic_rep_t_test <- function(
  mod_track,
  mod_num_rep,
  input_track,
  input_num_rep,
  bin_size=100
) {
  unbin_track <- function(track) do.call(
    c,
    sapply(
      track,
      \(v) subsample_time_series(v),
      simplify = FALSE
    )
  )
  data <- tibble(
    mu_1 = unbin_track(mod_track),
    sd_1 = unbin_track(mod_track %>% attr("standard_deviation")),
    count_1 = mod_num_rep,
    mu_2 = unbin_track(input_track),
    sd_2 = unbin_track(input_track %>% attr("standard_deviation")),
    count_2 = input_num_rep
  ) %>%
    rowwise() %>%
    mutate(
      df = WelchSatter(c(sd_1, sd_2), df = c(count_1 - 1, count_2 - 1))$ws.df
    ) %>%
    ungroup() %>%
    mutate(
      sd_test = sqrt(sd_1^2 + sd_2^2),
      t = (mu_1 - mu_2) / sd_test,
      p = pt(t, df, lower.tail = FALSE) %>% replace(t < 0, 1),
      q = p.adjust(p)
    )
}

chic_rle_make_rowData <- function(rles, bin_size) {
  bind_rows(
    mapply(
      \(n, vec) rbind(
        if (length(vec) > bin_size)
          tibble(
            start = seq(0, length(vec) - bin_size, by = bin_size),
            end = start + bin_size
          )
        else tibble(start = numeric(0), end = numeric(0)),
        if ((length(vec) %% bin_size) != 0)
          tibble(
            start = length(vec) - (length(vec) %% bin_size),
            end = length(vec)
          )
        else tibble(start = numeric(0), end = numeric(0))
      ),
      names(rles),
      rles,
      SIMPLIFY = FALSE
    ),
    .id = "chr"
  )
}

# F-test of equality of variances of the input tracks (biological replicates).
# Are there many genomic windows where, when pairing up input tracks at the
# window, there is not homogeneity of variance? In this case,
# pct_flagged_genomic_windows will far exceed the 5% of genomic windows where
# this is expected by random chance. If pct_flagged_genomic_windows passes this
# QC test, then we can estimate the standard deviation for our t-test pooling
# the input tracks together.
mean_coverage_track_var_test <- function(track1, track2, bin_size = 100) {
  mu1 <- track1 %>% subsample_unlist_rle_list(bin_size = bin_size)
  mu2 <- track2 %>% subsample_unlist_rle_list(bin_size = bin_size)
  series1 <- track1 %>% attr("standard_deviation") %>% subsample_unlist_rle_list(bin_size = bin_size)
  series2 <- track2 %>% attr("standard_deviation") %>% subsample_unlist_rle_list(bin_size = bin_size)
  square <- \(x) x^2
  pval <- pf(
    square(series1 / series2),
    attr(track1, "n") - 1,
    attr(track2, "n") - 1
  )
  pval <- ifelse(
    pval < 0.5,
    2 * pval,
    2 - 2 * pval
  ) %>%
    replace(pmin(mu1, mu2) < 0.1, NA)
  list(pvals = pval, pct_flagged_genomic_windows = mean(pval < 0.05))
}

estimate_variance_input_tracks <- function(input_tracks, bin_size=50) {
  square <- \(x) x^2
  sd <- input_tracks %>%
    sapply(
      \(tr) tr %>%
        attr("standard_deviation") %>%
        subsample_unlist_rle_list(bin_size = bin_size)
    )
  sqrt(rowMeans(square(sd)))
}

plot_chic_anova <- function(mod_track, input_tracks, bin_size=50) {
  pair_tracks <- combn(length(input_tracks), 2, simplify=FALSE)
  pvals <- pair_tracks %>%
    sapply(
      \(inds) ecdf(mean_coverage_track_var_test(input_tracks[[inds[1]]], input_tracks[[inds[2]]], bin_size=bin_size)$pvals),
      simplify=FALSE
    )
  pair_names <- sapply(
    pair_tracks,
    \(v) v %>% sapply(
        \(idx) names(input_tracks)[idx] %>%
          str_extract("H3K[0-9]+")
      ) %>%
      replace(1, paste0("Input: ", .[1])) %>%
      as.list %>%
      append(list(sep="/")) %>%
      do.call(paste, .)
  )
  plot_data <- pvals %>% mapply(
    \(n, fn) tibble(
      comparison = n,
      F = seq(0, 1, length.out = 101),
      p = fn(F)
    ),
    pair_names, .,
    SIMPLIFY = FALSE
  ) %>% do.call(rbind, .)

  # Mod/Input f-test
  sd_est <- estimate_variance_input_tracks(input_tracks, bin_size)
  sd_mod <- mod_track %>% attr("standard_deviation") %>% subsample_unlist_rle_list(bin_size = bin_size)
  square <- \(x) x^2
  pval_est <- pf(
    square(sd_mod / sd_est),
    attr(mod_track, "n") - 1,
    sum(
      sapply(input_tracks, \(tr) attr(tr, "n") - 1)
    )
  )
  plot_data <- plot_data %>% rbind(
    tibble(
      comparison = "Mod/Inputs",
      F = seq(0, 1, length.out=101),
      p = ecdf(pval_est)(F)
    )
  )

  ggplot(
    plot_data, aes(F, p)
  ) + geom_line(aes(group=comparison, color=comparison), \(df) df %>% subset(comparison != "Mod/Inputs")) + geom_line(aes(group=comparison, color=comparison), \(df) df %>% subset(comparison == "Mod/Inputs"), linetype="dotted") + scale_color_manual(
    values = c("purple", "cyan", "orange", "magenta")
  ) + geom_line(
    data=data.frame(F=c(0,1), p=c(0,1)),
    linetype="dashed"
  ) + geom_line(
    data=data.frame(
      F=c(0.025,0.025,0, 0.975,0.975,1),
      p=c(0.025,0,0.025, 1,0.975,0.975),
      group=rep(c("b", "t"), each=3)
    ),
    aes(group=group),
    color="blue"
  ) + theme_bw() + theme(
  ) + coord_cartesian(
    expand=FALSE
  )
}

chic_track_welch <- function(mod_track, input_tracks, mark_name, bin_size=100) {
  track_list <- append(list(mod=mod_track), input_tracks)
  stopifnot(mark_name %in% names(track_list))
  bin_data <- as_tibble(sapply(track_list, list, simplify=FALSE)) %>% reframe(
    across(
      everything(),
      .fns = list(
        mu = \(rles) rles[[1]] %>%
          subsample_unlist_rle_list(bin_size = bin_size),
        sd = \(rles) rles[[1]] %>% attr("standard_deviation") %>%
          subsample_unlist_rle_list(bin_size = bin_size)
      )
    )
  )
  mu_values_matrix <- bin_data %>%
    dplyr::select(all_of(paste0(c("mod", names(input_tracks)), "_mu"))) %>%
    as.matrix
  sd_values_matrix <- bin_data %>%
    dplyr::select(all_of(paste0(c("mod", names(input_tracks)), "_sd"))) %>%
    as.matrix
  sd_values_list <- bin_data %>%
    dplyr::select(all_of(paste0(c("mod", names(input_tracks)), "_sd"))) %>%
    # Place into list of row entries
    apply(1, identity, simplify=FALSE)
  df_values <- sapply(
    track_list,
    \(rles) attr(rles, "n") - 1
  )
  # Linear combination of random variables: The mod track has a coefficient of
  # 1, while we are averaging the processed input tracks together.
  ci_values = c(1, -1 / length(input_tracks) %>% rep(length(input_tracks)))
  epsilon <- 1e-10
  square <- \(x) x^2
  stat_data <- tibble(
    df_welch = mapply(
      \(sd, ci, df) WelchSatter(ui = sd, ci = ci, df = df)$ws.df,
      sd_values_list,
      list(abs(ci_values)),
      list(df_values)
    ),
    df_welch_simple = sapply(
      sd_values_list,
      \(v) WelchSatter(
        ui = v[c("mod_sd", paste0(mark_name, "_sd"))],
        df = df_values[c("mod", mark_name)]
      )$ws.df
    ),
    sd = sqrt(as.numeric(sd_values_matrix^2 %*% ci_values^2)),
    sd_simple = sqrt(
      bin_data %>% pull("mod_sd") %>% square
      + bin_data %>% pull(paste0(mark_name, "_sd")) %>% square
    ),
    t = (
      bin_data %>% pull("mod_mu")
      - bin_data %>% pull(paste0(mark_name, "_mu"))
    ) / sd,
    enrichment = bin_data %>% pull("mod_mu")
    / bin_data %>% pull(paste0(mark_name, "_mu")),
    p = pt(t, df_welch, lower.tail = FALSE) %>% replace(t <= epsilon, NA)
  )
  row_data <- chic_rle_make_rowData(mod_track, bin_size)
  peak_enrichment_threshold <- 1.5
  peak_caller <- as(
    interaction(
      row_data$chr,
      ifelse(stat_data$enrichment >= peak_enrichment_threshold, "enriched", "bg")
    ),
    "Rle"
  )
  peak_caller <- Rle(
    factor(
      ifelse(
        str_ends(peak_caller@values, "enriched"),
        paste0(as.character(peak_caller@values), ".", seq_along(peak_caller@values)),
        as.character(peak_caller@values)
      )
    ),
    peak_caller@lengths
  )
  peak_table <- cbind(
    chic_rle_make_rowData(mod_track, bin_size),
    stat_data,
    peak_caller = as.factor(peak_caller)
  ) %>%
    group_by(peak_caller) %>%
    summarise(
      chr = min(chr),
      start = min(start),
      end = max(end),
      enrichment = max(enrichment),
      p = min(c(p, 1), na.rm=T)
    ) %>%
    filter(
      grepl("enriched", as.character(peak_caller))
    ) %>%
    dplyr::select(!all_of("peak_caller")) %>%
    mutate(q = p.adjust(p, "BH"))
}

write_chic_peaks <- function(peak_table_list, output_path) {
  dir.create(dirname(output_path), recursive = TRUE, showW = FALSE)
  peaks_bed <- peak_table_list %>%
    sapply(
      \(tab) tab %>% subset(q < 0.1), simplify=F
    ) %>%
    bind_rows(.id = "mark") %>%
    mutate(
      mark = paste0(
        mark,
        " (",
        cut(
          q,
          c(0, 1e-4, 1e-3, 1e-2, 5e-2, 1e-1)
        ) %>%
          fct_recode(
            `****`='(0,0.0001]',
            `***`='(0.0001,0.001]',
            `**`='(0.001,0.01]',
            `*`='(0.01,0.05]',
            `~`='(0.05,0.1]'
          ),
        ")"
      )
    ) %>%
    subset(select=c(chr,start,end,mark))
  with_options(
    list(scipen=100),
    peaks_bed %>%
      write.table(output_path, sep="\t", quote = F, row.names = F, col.names = F)
  )
}