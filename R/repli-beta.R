# Repliseq Regression Model ----
repli_posterior_bar_colors <- append(
  chic_line_track_colors,
  list(
    Kc167 = "#C27A00",
    S2 = "#4087F4"
  )
)

# Gets SummarizedExperiment ready for Repliseq Regression.
# This is a la creating a DESeqDataSet. For this regression, it initializes
# metadata that will be used in the main computation workhorse.
init_beta_dm_experiment <- function(sce) {
  # We will pretend that the SummarizedExperiment is a SingleCellExperiment
  # - having an optional colData entry named sizeFactor.
  sf <- sce$sizeFactor
  # Initialize to counts per million.
  if (is.null(sf))
    sf <- colSums(assay(sce, "counts")) / 1000 / 1000
  # Col data: replication_value gives some keys to the S-phase quantified FACS
  # sorted fractions. It is most important that the levels represent FACS
  # sorted fractions that start near the beginning of S-phase (for the first
  # fraction) and lastly captures cells near the end of S-phase (for the last
  # factor level). The readout of the FPKM expected value for the entry in the
  # Repliseq expression matrix will depend on the factor level, counting out
  # uniformly spaced quantiles of timing. This is described in the likelihood
  # function for the vector of fragment counts.
  num_fractions <- length(levels(sce$replication_value))
  # Defines uniformly spaced quantiles of timing, these are needed for the
  # likelihood function.
  sce@metadata$beta_regression_cuts <- seq(
    0, 1, length.out=1+num_fractions
  )
  # An offset that is multiplied with the FPKM expected value to produce the
  # expected entry of fragment count in the FPKM data.
  # Zero is not allowed in extraDistr Dirichlet! (For an entry which is observed
  # to always be zero). Infinitesimal value works well, and is so small that it
  # does not influence the optimization process.
  # The size factors matrix must be rectangular (infinitesimal in case of a
  # missing logical spot for one biological replicate x one sorted fraction HTS
  # sample). This size factors matrix is multiplied by a regression expected
  # value FPKM which is a row vector (the same factor for every row of the
  # rectangular Repliseq observations).
  beta_regression_size_factors <- matrix(
    1e-40,
    nrow = length(levels(sce$rep)),
    ncol = num_fractions,
    dimnames = list(levels(sce$rep), levels(sce$replication_value))
  )
  # Linear operator. The row corresponds to an observation (column) in the assay
  # (e.g. Replicate 1 Early, Replicate 1 Early-Mid, ...). The column indices of
  # beta_regression_feed_y count out the entries in a rectangular Repliseq
  # observations matrix (the likelihood function is vector-valued and we handle
  # replicates of the HTS experiment). The column indices count out in
  # column-major order, so the first column encodes whether there is a
  # Replicate 1 Early observation, then the second column encodes whether there
  # is a Replicate 2 Early observation, then the third column (in case of 2
  # replicates) Replicate 1 Early-Mid observation. 
  # No observation (as part of the experimental design) is fine (missing data
  # which always appears as zero in one entry of the rectangular set of Repliseq
  # observed fragment counts). In that case, we manipulate
  # beta_regression_size_factors to drive the regression expected value for that
  # missing HTS sample to be infinitesimally small. Then it is fine that we
  # always observe zero fragments in that particular sample.
  beta_regression_feed_y <- matrix(
    0,
    nrow = ncol(sce),
    ncol = length(beta_regression_size_factors),
    dimnames = list(colnames(sce), NULL)
  )
  for (i in seq(ncol(sce))) {
    rep <- sce$rep[i]
    replication_value <- sce$replication_value[i]
    beta_regression_size_factors[
      as.numeric(rep),
      as.numeric(replication_value)
    ] <- sf[i]
    beta_regression_feed_y[
      i,
      as.numeric(rep)
      + length(levels(sce$rep)) * (as.numeric(replication_value) - 1)
    ] <- 1
  }
  sce@metadata$beta_regression_num_fractions <- num_fractions
  sce@metadata$beta_regression_size_factors <- beta_regression_size_factors
  sce@metadata$beta_regression_feed_y <- beta_regression_feed_y
  sce
}

# Our "plate" of observations enhances the biological replicates by also
# treating non-overlapping consecutive genomic windows (index in the experiment
# rows given by vector "i") as replicates observing the nascent DNA density in
# the region containing all of the windows. Thus, we produce a larger
# rectangular matrix of Repliseq observations than what we would see from the
# above metadata of the SummarizedExperiment.
beta_dm_pull_expr_plate <- function(exper, i)
  # Readout of the assay comes by taking the transpose of the
  # beta_regression_feed_y (as every column of beta_regression_feed_y encodes
  # one entry of the Repliseq observations, the cross product is taking the
  # transpose of beta_regression_feed_y and using matrix multiplication with the
  # counts on the right).
  crossprod(
    exper@metadata[["beta_regression_feed_y"]],
    t(assay(exper, "counts")[i,, drop=F])
  ) %>%
    array(
      dim = c(
        length(levels(exper$rep)),
        length(levels(exper$replication_value)),
        length(i)
      ),
      dimnames = list(
        levels(exper$rep),
        levels(exper$replication_value),
        i
      )
    ) %>%
    aperm(c(1, 3, 2)) %>%
    matrix(
      nrow = length(levels(exper$rep)) * length(i),
      ncol = length(levels(exper$replication_value)),
      dimnames = list(
        with(
          cross_join(tibble(index=i), tibble(rep=levels(exper$rep))),
          str_glue("{rep}@{index}")
        )
      )
    )

# Repliseq Dirichlet-Multinomial likelihood function. The likelihood is
# vector-valued (the observation is one row of the observation matrix Y). The
# Dirichlet-Multinomial distribution is conditioned on the vector sum (which is
# taken to be unvarying compared to the observed total fragment count, as a
# simplification). Crucially, the uniform cuts to regression timing are tested
# against the Beta distribution (parameterized as a "weighted coin toss"
# experiment by the parameters "heads" and "tails"). This defines a bar chart of
# percentage density of nascent DNA in the Repliseq fractions that has a mode RT
# sorted cell fraction for greatest density of fragments, and which should be
# unimodal (and the density declines dramatically on either side). This
# concentration around the mode is parameterized by `heads+tails`. The
# distribution is overdispersed, like DESeq2 and other Negative Binomial
# regression, so we need not only the 2 parameters that give the exact expected
# nascent DNA density in each sorted fraction, but also the overdispersion
# `theta`. Dividing the Beta function distribution of density in the
# fractions by `theta` produces overdispersion similar to what is modeled by
# Negative Binomial regression (by dividing the Negative Binomial parameter by
# `theta`). If we drive the Dirichlet-Multinomial parameters very large by using
# small `theta`, then we approximate a Multinomial distribution. For large
# `theta`, the small concentration parameter to the Dirichlet-Multinomial
# distribution results in very overdispersed predictions. Thus, the likelihood
# is highly comparable in character to the Negative Binomial regression
# likelihood, but is for vectors of sorted cell samples DNA libraries from one
# biological sample.
beta_dm_plain_likelihood_fn <- function(
  Y, beta_regression_cuts, beta_regression_size_factors,
  wts, heads, tails, theta
) {
  n <- rowSums(Y)
  beta_values <- cross_join(
    tibble(cut = beta_regression_cuts),
    tibble(heads = as.numeric(heads), tails = as.numeric(tails))
  )
  pb <- matrix(
    pbeta(beta_values$cut, beta_values$heads, beta_values$tails),
    nrow = length(heads)
  ) %>%
    rowDiffs
  pb <- pmax(pb, 1e-40)
  probs <- purrr::reduce(
    sapply(
      seq(nrow(Y)),
      \(j) {
        y <- Y[j, ]
        sf <- beta_regression_size_factors[j, ]
        mu_unnorm <- t(t(pb) * sf)
        dirparam <- (
          mu_unnorm
          /
          rowSums(mu_unnorm)
          /
          theta
        )
        ddirmnom(
          y,
          n[j],
          dirparam
        )^wts[j]
      },
      simplify=FALSE
    ),
    `*`
  )
}

# TODO comment this
beta_dm_prior_heads_tails_fn <- function(heads, tails) {
  prior <- dnorm(
    sqrt((heads - 1)^2 + (tails - 1)^2),
    mean = 2.5,
    sd = 1
  )
  prior[!(heads >= 1 & tails >= 1)] <- 0
  prior <- prior
}

# The prior on "heads" and "tails" as a separable filter of "X" (angle) and "Y"
# (heads + tails).
beta_dm_prior_filter <- function(repli.polar.coordinates) {
  angle <- repli.polar.coordinates$X
  rho <- repli.polar.coordinates$Y
  heads <- 1 + rho * sin(angle)
  tails <- 1 + rho * cos(angle)
  M <- beta_dm_prior_heads_tails_fn(heads, tails)
  s <- svd(M)
  Xsum <- ncol(M) / (pi/2)
  Xfactor <- Xsum / sum(abs(s$v[,1]))
  list(
    X = abs(s$v[,1]) * Xfactor,
    Y = s$d[1] / Xfactor * abs(s$u[,1])
  )
}

beta_dm_regression_likelihood <- function(
  exper, center_i, wts,
  repli.polar.coordinates,
  theta = NULL,
  xform_scale = 1, xform_center = 0
) {
  angle <- repli.polar.coordinates$X
  rho <- repli.polar.coordinates$Y
  heads <- 1 + rho * sin(angle)
  tails <- 1 + rho * cos(angle)
  beta_regression_cuts <- exper@metadata[["beta_regression_cuts"]]
  beta_regression_size_factors <- exper@metadata[["beta_regression_size_factors"]]
  i <- seq(
    pmax(1, center_i - floor(length(wts) / 2)),
    pmin(nrow(exper), center_i + floor(length(wts) / 2))
  )
  wts <- wts[
    seq(
      1 + floor(length(wts) / 2) + min(i) - center_i,
      1 + floor(length(wts) / 2) + max(i) - center_i
    )
  ]
  Y <- beta_dm_pull_expr_plate(exper, i)
  beta_regression_size_factors <- beta_regression_size_factors[
    rep(seq(nrow(beta_regression_size_factors)), length(i)),
  ]
  wts <- rep(wts, each=length(levels(exper$rep)))

  # Mutate the amount placed on heads or tails according to the xform.
  update_param <- heads > 1.001 & tails > 1.001
  pseudotime_value <- qlogis(
    ((heads - 1) / (heads + tails - 2))[update_param]
  )
  pseudotime_value <- pseudotime_value * xform_scale + xform_center
  updated_param_value <- plogis(pseudotime_value)
  heads[update_param] <- (
    1 + (heads + tails - 2)[update_param] * updated_param_value
  )
  tails[update_param] <- (
    1 + (heads + tails - 2)[update_param] * (1 - updated_param_value)
  )

  if (is.null(theta)) {
    probs <- beta_dm_plain_likelihood_fn(
      Y, beta_regression_cuts, beta_regression_size_factors,
      wts, heads, tails, theta = 0.02
    ) *
      beta_dm_prior_heads_tails_fn(heads, tails)
    param0 <- which.max(probs)
    if (is.matrix(heads)) {
      i <- param0 %% nrow(heads)
      j <- 1 + floor((param0 - 1) / nrow(heads))
      heads0 <- heads[i, j]
      tails0 <- tails[i, j]
    } else {
      heads0 <- heads[param0]
      tails0 <- tails[param0]
    }
    # We cannot actually find the caller of the random function which is leading
    # to "future.seed" warnings. It is likely in the optimize() algorithm. Even
    # then, we cannot readily reproduce the "beta_dm_regression_likelihood"
    # future.seed warning. with_seed will ensure that the result is
    # deterministic even if the "future.seed" warning does not know that we are
    # doing so.
    theta <- with_seed(
      0,
      optimize(
        \(theta) -log(
          beta_dm_plain_likelihood_fn(
            Y, beta_regression_cuts, beta_regression_size_factors,
            wts, heads0, tails0, theta
          )
        ),
        interval = c(0.001, 20)
      )$minimum
    )
  }

  probs <- beta_dm_plain_likelihood_fn(
    Y, beta_regression_cuts, beta_regression_size_factors,
    wts, heads, tails, theta
  )
  probs <- if (is.matrix(heads)) {
    probs %>%
      matrix(
        nrow = nrow(heads), ncol = ncol(heads), dimnames = dimnames(heads)
      )
  } else {
    probs
  }
  probs <- probs %>% `attr<-`("theta", value = theta)
}

beta_dm_regression_post_probability <- function(
  exper, center_i, wts,
  repli.polar.coordinates,
  prior_filter,
  theta = NULL,
  xform_scale = 1, xform_center = 0
) {
  grid <- beta_dm_regression_likelihood(
    exper, center_i, wts,
    repli.polar.coordinates,
    theta,
    xform_scale, xform_center
  )
  angle <- repli.polar.coordinates$X[1, ]
  rho <- repli.polar.coordinates$Y[, 1]
  tibble(
    value = sin(angle) / (sin(angle) + cos(angle)),
    prob = colSums(rho * prior_filter$Y * grid),
    theta = attr(grid, "theta")
  )
}

analyze_repli_experiment <- function(
  repli.experiment,
  repli.sliding.weights,
  repli.polar.coordinates,
  repli.prior.distribution,
  xform_scale,
  xform_center,
  theta = NULL
) {
  if (is.null(theta))
    theta <- numeric(0)
  # For each set of indices into the experiment (each chromosome),
  # parallel-apply the beta fitting function with smooth-weight plate.
  mapply(
    \(inds, thetas) {
      exper <- repli.experiment[inds, ]
      future_lapply(
        seq_along(inds),
        \(lookup) with(
          beta_dm_regression_post_probability(
            exper,
            lookup,
            repli.sliding.weights,
            repli.polar.coordinates,
            repli.prior.distribution,
            theta = if (length(thetas)) thetas[lookup] else NULL,
            xform_scale = xform_scale,
            xform_center = xform_center
          ),
          tibble(
            rowname = inds[lookup],
            value,
            prob,
            theta
          )
        )
      ) %>%
        do.call(rbind, .)
    },
    split(
      seq(nrow(repli.experiment)),
      rowData(repli.experiment)$seqnames
    ),
    split(
      theta,
      rowData(repli.experiment)$seqnames
    ),
    SIMPLIFY = FALSE
  ) %>%
    do.call(rbind, .)
}

quantify_repli_experiment <- function(prob, prior) {
  angle <- seq(0, pi/2, length.out=length(prob))
  timing <- sin(angle)/(sin(angle)+cos(angle))
  d_timing_d_angle <- 1/(1 + sin(2*angle))
  sum(timing * prob * prior * d_timing_d_angle) / sum(prob * prior * d_timing_d_angle)
}

bayes_factor_repli_experiment <- function(elementMetadata, prior) {
  angle <- seq(0, pi/2, length.out=nrow(elementMetadata))
  d_timing_d_angle <- 1/(1 + sin(2*angle))
  full_model <- prod(
    apply(
      elementMetadata,
      2,
      \(v) sum(v * d_timing_d_angle * prior)
    )
  )
  reduced_model <- sum(
    rowProds(as.matrix(elementMetadata)) * d_timing_d_angle * prior
  )
  full_model / reduced_model
}

repli_posterior_bayes_factor_all_pairs <- function(repli.posterior, prior) {
  repli.posterior %>%
    group_by(rowname) %>%
    summarise(
      Germline_Somatic = bayes_factor_repli_experiment(cbind(Germline, Somatic), prior),
      Germline_Kc167 = bayes_factor_repli_experiment(cbind(Germline, Kc167), prior),
      Germline_S2 = bayes_factor_repli_experiment(cbind(Germline, S2), prior),
      Somatic_Kc167 = bayes_factor_repli_experiment(cbind(Somatic, Kc167), prior),
      Somatic_S2 = bayes_factor_repli_experiment(cbind(Somatic, S2), prior),
      Kc167_S2 = bayes_factor_repli_experiment(cbind(Kc167, S2), prior),
      Dynamic_Static_Model = bayes_factor_repli_experiment(cbind(Germline, Somatic, Kc167, S2), prior),
    ) %>%
    subset(select = -rowname)
}

repli_timing_factor_all_pairs <- function(repli_timing) {
  tile_track <- GRanges(
    seqnames(repli_timing[[1]]),
    ranges(repli_timing[[1]]),
    seqlengths = seqlengths(repli_timing[[1]])
  )
  list(
    Germline_Somatic = GRanges(
      tile_track,
      Germline = repli_timing$Germline$score,
      Somatic = repli_timing$Somatic$score
    ),
    Germline_Kc167 = GRanges(
      tile_track,
      Germline = repli_timing$Germline$score,
      Kc167 = repli_timing$Kc167$score
    ),
    Germline_S2 = GRanges(
      tile_track,
      Germline = repli_timing$Germline$score,
      S2 = repli_timing$S2$score
    ),
    Somatic_Kc167 = GRanges(
      tile_track,
      Somatic = repli_timing$Somatic$score,
      Kc167 = repli_timing$Kc167$score
    ),
    Somatic_S2 = GRanges(
      tile_track,
      Somatic = repli_timing$Somatic$score,
      S2 = repli_timing$S2$score
    ),
    Kc167_S2 = GRanges(
      tile_track,
      Kc167 = repli_timing$Kc167$score,
      S2 = repli_timing$S2$score
    ),
    GRanges(
      tile_track,
      Germline = repli_timing$Germline$score,
      Somatic = repli_timing$Somatic$score
    )
  )
}

# Nested peak calling (this version of the function was not integrated back into
# the cell-type differences but is used for static regions). In the region, at
# least 50% of the windows must be lit up, or else the region is excluded.
repli_post_test_nested_peak_calling <- function(wnd) {
  peaks <- GenomicRanges::reduce(wnd)
  peaks <- peaks %>%
    GenomicRanges::resize(width(.) + 10000, fix="center") %>%
    GenomicRanges::reduce() %>%
    GenomicRanges::resize(width(.) - 10000, fix="center") %>%
    subset(width(.) >= 20000)
  mapping <- findOverlaps(peaks, wnd) %>% as("List")
  coverage <- sapply(mapping, \(mapping) sum(width(wnd)[mapping])) / width(peaks)
  peaks[coverage >= 0.5]
}

repli_timing_nested_peak_calling <- function(bayes_factor, timings, hypothesis_testing_factor) {
  tile_track <- GRanges(
    seqnames(timings),
    ranges(timings),
    seqlengths = seqlengths(timings)
  )
  peaks <- GenomicRanges::reduce(
    tile_track[bayes_factor >= hypothesis_testing_factor]
  )
  seqlengths(peaks) <- NA
  peaks <- peaks %>%
    GenomicRanges::resize(width(.) + 10000, fix="center") %>%
    GenomicRanges::reduce() %>%
    GenomicRanges::resize(width(.) - 10000, fix="center") %>%
    subset(width(.) >= 20000)
  diff_timing <- -rowDiffs(as.matrix(elementMetadata(timings)))
  peaks_timing <- findOverlaps(
    peaks,
    tile_track
  ) %>%
    sapply(
      \(inds) diff_timing[inds]
    )
  peaks$NegDiff <- sapply(peaks_timing, min)
  peaks$PosDiff <- sapply(peaks_timing, max)
  names(peaks) <- str_replace(
    make.unique(as.character(seqnames(peaks))),
    "\\.|$",
    paste0(
      ".",
      sign((peaks$NegDiff + peaks$PosDiff) / 2) %>%
        factor(c("1", "-1")) %>%
        `levels<-`(
          value = paste0(colnames(elementMetadata(timings))[1], c("Earlier", "Later"))
        ),
      "."
    )
  ) %>%
    str_replace(
      "\\.$",
      ""
    )
  peaks
}

plot_posterior <- function(repli.posterior, repli.polar.coordinates) {
  angle <- repli.polar.coordinates$X[1, ]
  d_timing_d_angle <- 2/(1 + sin(2*angle))
  timing <- 1 - 2 * (sin(angle)/(sin(angle)+cos(angle)))
  x <- c(
    1,
    1,
    rep(timing[-1] - 0.5 * diff(timing), each = 2),
    -1,
    -1
  )
  polygons <- repli.posterior %>%
    group_by(celltype) %>%
    reframe(
      x,
      y = c(
        0,
        prob %>% `/`(sum(.)*2) %>% rep(each = 2),
        0
      )
    )
  X <- seq(-0.995, 0.995, by=0.01)
  Y <- polygons %>%
    group_by(celltype) %>%
    reframe(
      as_tibble(
        approx(x, y, xout=X)
      )
    )
  unnormalized_posterior <- 0.01 * (split(Y$y, Y$celltype) %>% sapply(sum))
  polygons <- polygons %>%
    group_by(celltype) %>%
    mutate(y = y / unnormalized_posterior[celltype])
  ymax <- max(polygons$y) * 1.05
  ggplot(
    polygons,
    aes(x, y, fill = celltype)
  ) +
    geom_polygon(alpha = 0.75) +
    scale_fill_manual(
      values = repli_posterior_bar_colors %>% setNames(str_to_title(names(.))),
      guide = guide_legend(NULL, override.aes = list(alpha = 1))
    ) +
    scale_x_continuous(labels = as.numeric) +
    coord_cartesian(
      c(1, -1), c(0, ymax), expand = FALSE
    ) +
    labs(
      x = "Timing",
      y = "P(Timing)"
    ) +
    theme(
      aspect.ratio = 1/2,
      legend.position = "bottom",
      plot.margin = margin(0, 5.5, 0, 5.5),
    )
}

plot_fpkm_bayes_mmse <- function(
  fpkm_bars,
  bayes_mmse_param,
  repli.posterior,
  repli.polar.coordinates
) {
  barplot <- ggplot(
    tibble(x = names(repli_level_colors) %>% factor(., .), y = fpkm_bars),
    aes(x, y, fill = x)
  ) +
    geom_bar(stat = "identity") +
    scale_fill_manual("Fraction", values = unlist(repli_level_colors)) +
    coord_cartesian(c(0.4, 4.6), c(0, 1.05 * max(fpkm_bars)), expand = FALSE) +
    labs(
      x = "Fraction",
      y = "FPKM"
    ) +
    theme(legend.position = "bottom")
  predict_beta <- tibble(
    x = seq(0, 1, by = 0.01),
    y = dbeta(x, 1 + bayes_mmse_param[2] * sin(bayes_mmse_param[1]), 1 + bayes_mmse_param[2] * cos(bayes_mmse_param[1]))
  )
  predict_beta <- tibble(
    x = c(0, predict_beta$x),
    y = c(
      0,
      predict_beta$y[2] - 2 * (predict_beta$y[3] - predict_beta$y[2]),
      predict_beta$y[-1]
    )
  )
  predict_point <- tibble(
    x = sin(bayes_mmse_param[1]) / (sin(bayes_mmse_param[1]) + cos(bayes_mmse_param[1])),
    y = approx(predict_beta$x, predict_beta$y, xout = sin(bayes_mmse_param[1]) / (sin(bayes_mmse_param[1]) + cos(bayes_mmse_param[1])))$y
  )
  predict_polygon <- bind_rows(
    list(
      E = rbind(
        subset(predict_beta, x <= 0.25),
        tibble(x = 0.25, y = 0)
      ),
      EM = rbind(
        tibble(x = 0.25, y = 0),
        subset(predict_beta, between(x, 0.25, 0.5)),
        tibble(x = 0.5, y = 0)
      ),
      ML = rbind(
        tibble(x = 0.5, y = 0),
        subset(predict_beta, between(x, 0.5, 0.75)),
        tibble(x = 0.75, y = 0)
      ),
      L = rbind(
        tibble(x = 0.75, y = 0),
        subset(predict_beta, between(x, 0.75, 1))
      )
    ),
    .id = "name"
  )
  fitplot <- ggplot(
    predict_beta,
    aes(x, y)
  ) +
    geom_polygon(aes(fill = name), data = predict_polygon) +
    geom_line() +
    geom_point(data = predict_point, size = 2) +
    scale_x_continuous(
      labels = \(v) 1 - 2*v
    ) +
    scale_fill_manual("Fraction", values = unlist(repli_level_colors)) +
    labs(x = "Timing", y = "Nascent DNA Density") +
    coord_cartesian(c(0, 1), c(0, 1.05 * max(predict_beta$y)), expand = FALSE) +
    theme(legend.position = "bottom")
  w <- 2.5
  h <- 1.25
  cbind(
    set_panel_size(barplot, w = unit(w, "in"), h = unit(h, "in")),
    set_panel_size(fitplot, w = unit(w, "in"), h = unit(h, "in")),
    set_panel_size(plot_posterior(repli.posterior, repli.polar.coordinates), w = unit(w, "in"), h = unit(h, "in"))
  )
}
