beta_prior_draws <- function(n = 200, shape = 4, rate = 0.5) {
  alpha <- rgamma(n * 2, shape, rate) %>%
    subset(. > 1) %>%
    head(n)
  beta <- rgamma(n * 2, shape, rate) %>%
    subset(. > 1) %>%
    head(n)
  x <- seq(0, 1, by = 0.02)
  y <- tibble(alpha, beta) %>%
    rowwise() %>%
    reframe(
      name = str_glue("{alpha}@{beta}"),
      x,
      y = dbeta(x, alpha, beta)
    )
  ggplot(y, aes(x, y, group = name, color = name)) +
    geom_line() +
    scale_color_hue(guide = NULL) +
    theme_cowplot() +
    coord_cartesian(c(0, 1), c(-0.01, 7.5), expand = F)
}

init_beta_dm_experiment <- function(sce) {
  # We will pretend that the SummarizedExperiment is a SingleCellExperiment
  # - having an optional colData entry named sizeFactor.
  sf <- sce$sizeFactor
  if (is.null(sf))
    sf <- colSums(assay(sce, "counts")) / 1000 / 1000
  num_fractions <- length(levels(sce$replication_value))
  sce@metadata$beta_regression_cuts <- seq(
    0, 1, length.out=1+num_fractions
  )
  # Zero is not allowed in extraDistr Dirichlet! (For an entry which is observed
  # to always be zero). Infinitesimal value works well, and is so small that it
  # does not influence the optimization process.
  beta_regression_size_factors <- matrix(
    1e-40,
    nrow = length(levels(sce$rep)),
    ncol = num_fractions,
    dimnames = list(levels(sce$rep), levels(sce$replication_value))
  )
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

# Pulls the "Y" matrix of DM-distributed data (tuples of sample counts). As we
# are interested in a rolling view of "Y", then "i" is a numeric vector of the
# features in the SummarizedExperiment to extract.
beta_dm_pull_expr_plate <- function(exper, i)
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

bayes_factor_repli_experiment <- function(prob1, prob2, prior) {
  angle <- seq(0, pi/2, length.out=length(prob1))
  d_timing_d_angle <- 1/(1 + sin(2*angle))
  full_model <- sum(
    prob1 * d_timing_d_angle * prior
  ) * sum(
    prob2 * d_timing_d_angle * prior
  )
  reduced_model <- sum(prob1 * prob2 * d_timing_d_angle * prior)
  full_model / reduced_model
}
