illustrate_coverage_poisson <- function(coverage) {
  dots_x <- sample(length(coverage), 200, rep=T, prob = coverage) + runif(200)
  dots_y <- runif(200)
  ggpubr::ggarrange(
    ggplot(data.frame(x=seq_along(coverage), y=coverage), aes(x, y))
    + geom_line()
    # Red line representing the window
    + geom_line(
      aes(group=g),
      data=tribble(~x, ~y, ~g, 110, -Inf, "a", 110, Inf, "a", 210, -Inf, "b", 210, Inf, "b"),
      color="red"
    )
    + scale_x_continuous(expand=c(0, 0), breaks = c(100, 200), labels = \(x) x * 10)
    + scale_y_continuous(breaks=NULL)
    + theme_cowplot()
    + labs(x=NULL, y=NULL),
    ggplot(data.frame(x=dots_x, y=dots_y), aes(x, y))
    + geom_point()
    + geom_rect(
      aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
      # arbitrary ymin and ymax representing the sequencing depth as a random
      # Poisson-distributed variable
      data.frame(xmin=110, xmax=210, ymin=0.3, ymax=0.9, x=-Inf, y=-Inf),
      color="red",
      fill="transparent"
    )
    + scale_x_continuous(expand=c(0, 0), breaks = c(100, 200), labels = \(x) x * 10)
    + scale_y_continuous(breaks=NULL, expand=rep(0.02, 2))
    + theme_bw()
    + labs(x=NULL, y=NULL),
    ncol = 1,
    heights=c(1,0.5)
  )
}

illustrate_poisson_variable <- function(coverage) {
  dots <- sapply(
    1:3,
    \(i) data.frame(
      x = sample(length(coverage), 200, rep=T, prob = coverage) + runif(200),
      y = runif(200)
    ) %>%
      ggplot(aes(x, y))
    + geom_point()
    + scale_x_continuous(expand=c(0, 0), breaks = c(100, 200), labels = \(x) x * 10)
    + scale_y_continuous(breaks=NULL, expand=rep(0.02, 2))
    + theme_bw()
    + labs(x=NULL, y=NULL),
    simplify=FALSE
  )
  ggarrange(
    plotlist = dots,
    ncol = 1
  )
}

chic_illustrate_mean_variance <- function(chic.replicates, window=10000) {
  chic.table = chic.replicates %>%
    sapply(
      \(raw) raw[1:sum(chr.lengths)] %>%
        smooth_sparse_vector_macs(
          Rle(names(chr.lengths), chr.lengths),
          window,
          sample_size = 200,
          norm=F
        ) %>%
        sapply(list, simplify=F) %>%
        as_tibble,
      simplify=F
    ) %>%
    do.call(rbind, .)
  chic_illustrate_mean_variance_from_rle_table(chic.table)
}

chic_illustrate_mean_variance_gaussian <- function(chic.replicates, bw) {
  chic.table = chic.replicates %>%
    sapply(
      \(raw) raw[1:sum(chr.lengths)] %>%
        smooth_sparse_vector_to_rle_list(
          Rle(names(chr.lengths), chr.lengths),
          bw = bw,
          sample_size = 200
        ) %>%
        sapply(list, simplify=F) %>%
        as_tibble,
      simplify=F
    ) %>%
    do.call(rbind, .)
  chic_illustrate_mean_variance_from_rle_table(
    chic.table,
    # Cannot predict mean-variance relationship of the Gaussian noise.
    slope = NA
  )
}

chic_illustrate_s2_mean_variance <- function(s2_bigwigs, smooth_length=1000, bin_size=200) {
  chic.replicates <- s2_bigwigs %>%
    sapply(
      \(f) f %>% import("bw") %>% with(coverage(., weight=score)) %>%
        `*`(smooth_length / bin_size)
    )
  chic.table = chic.replicates %>%
    sapply(
      \(rle_list) rle_list[names(chr.lengths)] %>%
        sapply(list, simplify=F) %>%
        as_tibble,
      simplify=F
    ) %>%
    do.call(rbind, .)
  chic.table[3,] <- apply(
    chic.table,
    2,
    \(rep_rles) list(
      rle_1 = rep_rles[[1]],
      rle_2 = rep_rles[[2]],
      sample_points = sapply(
        rep_rles,
        \(rle) 1 + cumsum(c(0, rle@lengths[-length(rle@lengths)]))
      ) %>% setNames(NULL) %>% do.call(union, .) %>% sort
    ) %>%
      with(
        Rle(
          values = rpois(
            length(sample_points),
            as.numeric(rle_1[sample_points] + rle_2[sample_points]) / 2
          ),
          lengths = diff(c(sample_points, length(rle_1)))
        )
      ) %>%
      list,
    simplify=FALSE
  )
  chic_illustrate_mean_variance_from_rle_table(chic.table)
}

chic_illustrate_mean_variance_from_rle_table <- function(chic.table, slope=1) {
  chic.mat = apply(
    chic.table,
    2,
    \(rles) rles %>%
      sapply(
        \(rle) rle[sample(length(rle), length(rle) / 5000)] %>% as.numeric,
        simplify=F
      ) %>%
      do.call(cbind, .),
    simplify=FALSE
  ) %>%
    do.call(rbind, .)
  mu_limit <- quantile(rowMeans(chic.mat), 0.9)
  data.frame(mu=chic.mat %>% rowMeans, s2=chic.mat %>% rowVars) %>%
    ggplot(
      aes(mu, s2)
    ) + (# rasterise(
    geom_point(shape=20, size=0.25, alpha=0.1)# ,
    # dpi=30
   ) + geom_smooth() + geom_line(
    data=tribble(~mu, ~s2, 0,0, mu_limit, slope * mu_limit), color="red", linewidth=1
  ) + coord_cartesian(
    c(0, mu_limit),
    c(0, mu_limit^2 * 0.3),
    expand = FALSE
  ) + theme_bw() + theme(
    aspect.ratio = 2
  ) + labs(
    # x bar: bar is symbol for sample mean.
    x = "x\u0305",
    y = bquote(s^"2")
  )
}

plot_chic_rect_window_replicates <- function(chic_list, chr, xs) {
  xs <- range(xs)
  xs <- seq(xs[1], xs[2], by = 50)
  chic_list %>%
    sapply(\(rles) data.frame(x = xs, y = rles[[chr]][xs] %>% as.numeric), simplify=FALSE) %>%
    bind_rows(.id = "sample")
}