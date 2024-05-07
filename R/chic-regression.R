chic_perform_regression <- function(df, formula = "fpkm ~ molecule + replicate", contrasts_replicate = contr.helmert) {
  fpkm <- t(sapply(df$density, \(data) data$y))
  molecule <- df$molecule
  replicate <- df$rep
  if (!is.null(contrasts_replicate))
    contrasts(replicate) <- contrasts_replicate(length(levels(replicate)))
  mymanova <- manova(as.formula(formula), data = df)
}

chic_regression_fold_change <- function(mn, df) {
  which_coef_molecule <- grep("^molecule", rownames(mn$coefficients))
  coverage <- mn$coefficients["(Intercept)", ]
  fold_change <- (
    mn$coefficients["(Intercept)", ]
    +
    mn$coefficients[which_coef_molecule, ]
  ) / mn$coefficients["(Intercept)", ]
  reframe(df, chr, x, coverage, fold_change)
}

feature_lengths_to_sliding_windows <- function(feature.lengths, flybase.lengths, width, step) {
  granges <- slidingWindows(
    GRanges(mapply(\(n, l) str_glue("{n}:1-{l}"), names(feature.lengths), feature.lengths)),
    width=width,
    step=step
  )
  sum_features <- which(feature.lengths != flybase.lengths)
  granges[sum_features] <- GRanges(
    mapply(
      \(n, l) str_glue("{n}:1-{l}"),
      names(flybase.lengths)[sum_features],
      flybase.lengths[sum_features]
    )
  ) %>%
    split(
      names(flybase.lengths[sum_features]) %>%
        factor(., .)
    )
  granges
}

track_kde_to_sliding <- function(fold_change, granges) {
  gr <- unlist(granges)
  coverage_chrs <- list()
  fold_change_chrs <- list()
  chrs_to_tile <- setdiff(
    unique(fold_change$chr),
    fold_change$chr[is.na(fold_change$x)]
  )
  coverage_do_not_tile <- fold_change %>%
    subset(is.na(x)) %>%
    pull(coverage, chr)
  fold_change_do_not_tile <- fold_change %>%
    subset(is.na(x)) %>%
    pull(fold_change, chr)
  for (n in names(granges)) {
    if (n %in% chrs_to_tile) {
      fc <- fold_change %>% subset(chr == n)
      loc <- granges[[n]]@ranges@start + granges[[n]]@ranges@width/2
      coverage_chrs <- coverage_chrs %>%
        append(
          list(
            approx(
              fc$x,
              fc$coverage,
              xout = loc
            )$y
          ) %>%
            setNames(n)
        )
      fold_change_chrs <- fold_change_chrs %>%
        append(
          list(
            approx(
              fc$x,
              fc$fold_change,
              xout = loc
            )$y
          ) %>%
            setNames(n)
        )
    }
  }
  coverage_chrs <- coverage_chrs %>%
    append(as.list(coverage_do_not_tile))
  gr$coverage <- unlist(coverage_chrs[names(granges)])
  fold_change_chrs <- fold_change_chrs %>%
    append(as.list(fold_change_do_not_tile))
  gr$fold_change <- unlist(fold_change_chrs[names(granges)])
  gr
}

coverage_interp_obs_granges <- function(gr, predicate, score_col = "fold_change") {
  for (n in names(chr.lengths)) {
    take_metadata <- as.logical(gr@seqnames == n)
    xs <- gr@ranges@start[take_metadata & predicate] + 1/2*gr@ranges@width[take_metadata & predicate]
    ys <- gr@elementMetadata[take_metadata & predicate, score_col]
    xout <- gr@ranges@start[take_metadata] + 1/2*gr@ranges@width[take_metadata]
    approx_out <- approx(
      xs,
      ys,
      xout,
      rule = 2
    )
    gr@elementMetadata[take_metadata, score_col] <- approx_out$y
  }
  gr
}

regress_experiment_replicates <- function(df, ...) {
  models <- tibble(
    replicate = levels(df$rep)
  ) %>%
    rowwise %>%
    mutate(
      model = chic_perform_regression(
        as.list(df) %>%
          append(
            list(
              replicate_subset = {
                rep_relevel <- fct_relevel(.$rep, replicate)
                ctr <- contr.sum(length(levels(rep_relevel)))
                ctr[, 1] <- 0
                contrasts(rep_relevel) <- ctr
                rep_relevel
              }
            )
          ),
        "fpkm ~ molecule + replicate_subset + rep"
      ) %>%
        list
    )
  effect_row <- head(grep("^rep[0-9]", rownames(models$model[[1]]$effects)), 1)
  fpkm <- t(sapply(df$density, \(data) data$y))
  nobs <- nrow(fpkm)
  mmmmmodels <- models
  models$model[[1]]$effects <- models$model[[1]]$effects[, 1000:1003]
  models$model[[2]]$effects <- models$model[[2]]$effects[, 1000:1003]
  models$model[[3]]$effects <- models$model[[3]]$effects[, 1000:1003]
  SSB_fit = sum(models$model[[1]]$effects[2, ]^2) / (nobs - 1)
  square <- \(x) x^2
  SSE = (
    models$model[[1]]$effects %>%
      subset(rownames(.) == "") %>%
      square %>%
      sum
  ) / (nobs - 1)
  result <- tibble(
    replicate = models$replicate,
    SSB_replicate = sapply(models$model, \(m) sum(m$effects[effect_row, ]^2)) / (nobs - 1),
    SSB_fit,
    SSE,
    SST = sum(colVars(fpkm))
  )
}