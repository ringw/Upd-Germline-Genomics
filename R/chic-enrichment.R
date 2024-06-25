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
  rles[sapply(rles, length) > bin_size] %>%
    sapply(
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
    n_1 = mod_track %>% attr("n"),
    count_1 = mod_num_rep,
    mu_2 = unbin_track(input_track),
    sd_2 = unbin_track(input_track %>% attr("standard_deviation")),
    n_2 = input_track %>% attr("n"),
    count_2 = input_num_rep
  ) %>%
    rowwise() %>%
    mutate(
      df = WelchSatter(c(sd_1, sd_2), df = c(count_1 - 1, count_2 - 1))$ws.df
    ) %>%
    ungroup() %>%
    mutate(
      sd_test = sqrt(sd_1^2 + sd_2^2),
      t_denom = sqrt(sd_1^2 / n_1 + sd_2^2 / n_2),
      t = (mu_1 - mu_2) / t_denom,
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

chic_track_welch <- function(mod_track, input_tracks, mark_name, bin_size=100) {
  track_list <- append(list(mod=mod_track), input_tracks)
  track_n <- track_list %>% sapply(\(track) track %>% attr("n"))
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
    # t-test grows with sqrt(n) (for the estimate of the mean for each
    # coefficient in the contrast).
    t_denom = sqrt(
      as.numeric(
        sd_values_matrix^2 %*% (ci_values^2 * track_n)
      )
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
          c(0, 1e-4, 1e-3, 1e-2, 5e-2, 1e-1, Inf)
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
    subset(select=c(chr,start,end,mark)) %>%
    arrange(chr, start, mark)
  with_options(
    list(scipen=100),
    peaks_bed %>%
      write.table(output_path, sep="\t", quote = F, row.names = F, col.names = F)
  )
  output_path
}

chic_track_generate_table_by_enrichment <- function(chic_df_list, enrichment_threshold = 1.5) {
  enrichment_track <- chic_df_list %>%
    sapply(\(df) df %>% with(Rle(enrichment, length))) %>%
    RleList
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
      chic_df_list %>% chic_df_list_to_pvalue_list,
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

chic_squeeze_var <- function(rle_lists, sample_size = 50) {
  stopifnot(type_sum(names(rle_lists)) == "chr")
  df <- sapply(rle_lists, \(obj) attr(obj, "n")) - 1

  rle_lengths <- rle_lists[[1]] %>% sapply(length)

  # For chr name, build the list of limma squeezeVar models.
  build_limma_model <- function(chr) {
    limma_chr <- tibble(
      begin = seq(0, rle_lengths[chr], by=sample_size),
      end = pmin(begin + sample_size, rle_lengths[chr]),
      length = end - begin,
      x = ceiling((begin + end) / 2)
    )
    square <- \(x) x^2
    chr_vars <- rle_lists %>% sapply(
      \(rle_list) rle_list[[chr]][limma_chr$x] %>% as.numeric,
      simplify=FALSE
    ) %>%
      do.call(cbind, .) %>%
      square
    limma_chr$limma_model <- apply(
      chr_vars,
      1,
      \(vars) squeezeVar(vars, df),
      simplify=FALSE
    )
    limma_chr
  }
  limma_models <- sapply(names(rle_lengths), build_limma_model, simplify=FALSE)

  var_post_tracks <- sapply(
    names(rle_lists) %>% setNames(seq_along(.), .),
    \(rle_model_ind) sapply(
      names(rle_lengths),
      \(chr) limma_models[[chr]] %>%
        rowwise() %>%
        mutate(
          var.post = limma_model$var.post[rle_model_ind]
        ) %>%
        with(Rle(var.post, length))
    ) %>%
      RleList
  )
  df_track <- sapply(
    limma_models,
    \(chr_model) chr_model %>%
      rowwise() %>%
      mutate(
        df.prior = limma_model$df.prior
      ) %>%
      with(Rle(df.prior, length))
  ) %>%
    RleList
  list(var = var_post_tracks, df = df_track)
}

chic_ttest <- function(
  mod_track,
  input_track,
  sd_track,
  df_track,
  sample_size = 50
) {
  rle_lengths <- sapply(mod_track, length)
  df_residual <- attr(input_track, "n") - 1
  t_test_results <- sapply(
    names(mod_track),
    \(chr)
      tibble(
        begin = seq(0, rle_lengths[chr], by=sample_size),
        end = pmin(begin + sample_size, rle_lengths[chr]),
        length = end - begin,
        x = ceiling((begin + end) / 2),
        enrichment = (mod_track[[chr]][x] / input_track[[chr]][x]) %>%
          as.numeric %>%
          replace(as.logical(input_track[[chr]][x] < 1), NA)
      ) %>%
      mutate(
        p = pt(
          (
            as.numeric(mod_track[[chr]][x])
            - as.numeric(input_track[[chr]][x])
          ) / as.numeric(sd_track[[chr]][x])
          / sqrt(
            1 / attr(mod_track, "n")
            + 1 / attr(input_track, "n")
          ),
          as.numeric(df_track[[chr]][x]) + df_residual,
          lower.tail = FALSE
        )
      ),
    simplify=FALSE
  )
}

chic_df_list_to_pvalue_list <- function(df_list) {
  df_list %>%
    sapply(\(df) df %>% with(Rle(p, length))) %>%
    RleList
}

chic_quantify_broad_peaks <- function(track, track_mask, features, enrichment = 1.5, resolution = 10) {
  track <- (track >= enrichment) * (
    sapply(
      track_mask,
      \(track_mask) as(
        as.logical(track_mask) %>%
          ifelse(1, 0),
        "Rle"
      )
    ) %>%
      RleList
  )
  feature_values <- features %>%
    filter(
      chr %in% names(chr.lengths),
      !is.na(start),
      !is.na(end),
      transcript.length >= 2000
    ) %>%
    split(.$chr) %>%
    mapply(
      \(chr_name, chr_features) {
        feature_names <- rownames(chr_features)
        chr_features <- chr_features %>%
          rowwise %>%
          mutate(all_pos = seq(min(start, end), max(start, end)) %>% list)
        feature_factor <- rep(
          rownames(chr_features),
          sapply(chr_features$all_pos, length)
        ) %>%
          factor(rownames(chr_features))
        feature_lookup <- track[[chr_name]][
          chr_features$all_pos %>% do.call(c, .)
        ]
        sapply(
          feature_lookup %>% split(feature_factor),
          mean
        ) %>%
          setNames(feature_names)
      },
      names(.),
      .,
      SIMPLIFY=FALSE
    ) %>%
    setNames(NULL) %>%
    do.call(c, .)
}

display_peak_location_stats <- function(lst) {
  lst %>%
    sapply(
      \(df) df$region %>% table %>% enframe, simplify=FALSE
    ) %>%
    bind_rows(.id = "mark") %>%
    mutate(
      mark = mark %>% factor(chic.mark.data$mark),
      name = name %>% factor(lst[[1]]$region %>% levels)
    ) %>%
    ggplot(
      aes(x = name, y = value, fill = mark)
    ) + geom_bar(
      stat = "identity", position = "dodge"
    ) + scale_y_continuous(
      expand = expansion(mult = c(0, 0.05))
    ) + scale_fill_viridis_d(
      option = "turbo", begin = 0.3, end = 0.9
    ) + labs(
      x = "Genomic Annotation",
      y = "Number of Peaks (q < 0.05)"
    )
}