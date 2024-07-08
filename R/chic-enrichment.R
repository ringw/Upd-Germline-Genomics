chic_quantify <- function(
  chic_experiment,
  molecule,
  replicate_nums,
  offset,
  global_overdispersion=TRUE,
  test_de=FALSE
) {
  if (length(unique(replicate_nums)) == 1)
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
            as.data.frame %>%
            cbind(
              p_peak=1, L2FC=0
            )
        )
    )
  fit <- glm_gp(
    as.matrix(chic_experiment@elementMetadata),
    # We put "molecule" and "rep" in our tar_map, so create new names. As
    # "molecule" coefs come first, our p_peak contrast will be the log fold
    # change between the two levels of "molecule".
    ~ 0 + mol + R,
    tibble(
      mol = molecule,
      R = replicate_nums %>%
        factor %>%
        `contrasts<-`(value = contr.helmert(length(levels(.))))
    ),
    size_factors = 1,
    offset = offset,
    overdispersion = list(TRUE, "global")[
      c(isFALSE(global_overdispersion), isTRUE(global_overdispersion))][[1]],
    overdispersion_shrinkage = list(TRUE, FALSE)[
      c(isFALSE(global_overdispersion), isTRUE(global_overdispersion))][[1]],
    verbose = TRUE
  )
  GRanges(
    seqnames(chic_experiment),
    ranges(chic_experiment),
    seqlengths = seqlengths(chic_experiment)
  ) %>%
    `elementMetadata<-`(
      value = with(
        if (test_de)
          test_de(
            fit,
            c(-1, 1) %>% c(rep(0, ncol(fit$Beta) - length(.))),
            verbose=T
          )
        else list(pval=1, lfc=0),
        tibble(
          as.data.frame(exp(fit$Beta)) %>%
            replace(. < 1e-8, 0) %>%
            rename_with(~ str_glue("score.{.}")),
          p_peak=pval,
          L2FC=lfc
        )
      )
    ) %>%
    `metadata<-`(
      value = list(
        overdispersions=fit$overdispersions,
        overdispersion_shrinkage_list=fit$overdispersion_shrinkage_list,
        deviances=fit$deviances,
        ridge_penalty=fit$ridge_penalty,
        model_matrix=fit$model_matrix,
        design_formula=fit$design_formula
      )
    )
}

enrich_int_list <- function(p_peak_granges, int_list) {
  sapply(
    int_list,
    \(w) with(
      elementMetadata(gr)[w, ],
      min(
        p_peak[L2FC > 0],
        if (length(p_peak) > 0) 1 else numeric(0),
        na.rm=T
      ) %>%
        replace(!is.finite(.), NA)
    )
  )
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