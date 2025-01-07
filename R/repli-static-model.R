# Static-timing analysis of the Repliseq Joint Distribution.
# Lookup for static timing: Choosing rectangles
# (Angle variable / logit transform) such that with these coordinates into the
# joint distribution, all logistic coefficients (except the intercept) are <=1.
# These form the hyper-rectangles that are used for integrating (rectangle
# method). The integral is compared against the marginal likelihood of the
# original model, which is the product from 4 independent distributions.
build_static_timing_lookup_table <- function(logit_parameter, cell_type_names, cutoff=1) {
  i <- seq_along(logit_parameter)
  rows <- NULL
  for (n in cell_type_names) {
    tbl <- as_tibble(setNames(list(i), n))
    if (is.null(rows))
      rows <- tbl
    else
      rows <- rows %>% cross_join(tbl)
  }
  rows <- tibble(
    rows,
    Intercept = logit_parameter[as.matrix(rows[1:length(cell_type_names)])] %>%
      matrix(nrow=nrow(rows)) %>%
      rowMeans(),
    Effect_Quant = sapply(
      rows[1:length(cell_type_names)],
      \(i) abs(logit_parameter[i] - Intercept) %>% replace(is.na(.) | !is.finite(.), 0)
    ) %>%
      rowMaxs()
  )
  x <- subset(rows[1:4], !is.na(rows$Intercept) & is.finite(rows$Intercept) & rows$Effect_Quant <= cutoff)
  x <- as.matrix(x)
}

# Calculate P(Static) for the 4 cell types. Arguments:
# post: Fitted posterior distributions for each cell type in columns, evaluated
# at certain timing values. The posterior is un-normalized, the columns times
# prior do not need to sum to 1.
# lut: Lookup table defining a region of the joint distribution.
# d_timing_d_angle: The prior is no longer uniform (prior was uniform for angles
# in subset of Quadrant I from origin (1,1) for parameters Alpha, Beta). This is
# an un-normalized prior distribution that has not been multiplied to the
# posterior distribution yet. The implementation of this un-normalized prior is
# a first derivative that represents the change of variables from the angle to
# the Beta distribution mode parameter.
repli_static_model_probability <- function(post, lut, d_timing_d_angle) {
  post_lookup <- lut +
    rep(seq(0, ncol(post) - 1) * nrow(post), each = nrow(lut))
  logNumer <- rowSums(
    log(
      matrix(
        post[post_lookup] * d_timing_d_angle[lut],
        nrow = nrow(post_lookup)
      )
    )
  ) %>%
    logSumExp()
  logDenom <- sum(
    apply(
      post, 2, \(v) log(sum(v * d_timing_d_angle))
    )
  )
  exp(logNumer - logDenom)
}

# Apply repli_static_model_probability.
# repli.post.objs: List of repli.posterior tibbles which represent one cell type
# and all of the entries in the GRanges.
analyze_repli_static_model_probability <- function(repli.post.objs, repli.static.timing.table, lut) {
  npost <- sum(repli.post.objs[[1]]$rowname == 1)
  m <- tail(repli.post.objs[[1]]$rowname, 1)
  inds <- sapply(
    seq(m),
    \(i) repli.static.timing.table$i + (i - 1) * npost,
    simplify = FALSE
  )
  objs <- sapply(
    inds,
    \(i) sapply(
      repli.post.objs,
      \(obj) obj$prob[i]
    ),
    simplify = FALSE
  )
  future_sapply(
    objs,
    repli_static_model_probability,
    lut = lut,
    d_timing_d_angle = repli.static.timing.table$D_Timing_D_Angle
  )
}
