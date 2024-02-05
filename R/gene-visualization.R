determine_gene_list <- function(sctransform_quantile, supplemental_bulk_fpkm) {
  gene_lists <- sapply(
    c("germline", "somatic"),
    \(n) union(
      rownames(sctransform_quantile[[n]]) %>% subset(sctransform_quantile[[n]][, "90%"] > 1),
      supplemental_bulk_fpkm %>% subset(percent_expressed > 0.75) %>% rownames
    ),
    simplify = FALSE
  )
}

determine_two_way_venn_area <- function(n1=1500, n2=1000, n_intersect=500) {
  r1 = sqrt(n1 / pi)
  r2 = sqrt(n2 / pi)
  dmax = r1 + r2
  searchsquare = r1 + 2 * r2

  intersection_mc <- function(d) {
    x = runif(100000, min = -searchsquare, max = searchsquare)
    y = runif(100000, min = -searchsquare, max = searchsquare)
    mean(
      (x^2 + y^2) < r1^2
      & ((x-d)^2 + y^2) < r2^2
    ) * (2 * searchsquare)^2
  }

  optim(0, \(d) (intersection_mc(d) - n_intersect)^2, method="Brent", lower=0, upper=dmax)$par
}

plot_gene_lists <- function(gene_lists) {
  circles <- tribble(
    ~ Identity, ~ x0, ~ r,
    "Germline", 0, sqrt(length(gene_lists$germline) / pi),
    "Somatic", 0, sqrt(length(gene_lists$somatic) / pi)
  )
  intersection <- length(purrr::reduce(gene_lists, intersect))
  circles[2, "x0"] <- determine_two_way_venn_area(length(gene_lists$germline), length(gene_lists$somatic), intersection)
  circle_xlim <- c(-circles$r[1] - 3, circles$x0[2] + circles$r[2] + 3)
  circle_ylim <- c(-max(circles$r) - 3, max(circles$r) + 3)
  g <- circles %>%
    ggplot(
      aes(x0=x0, y0=0, r=r, fill=Identity)
    ) + ggforce::geom_circle(
      alpha = 0.3
    ) + scale_x_continuous(
      breaks = scales::breaks_width(5), labels = NULL,
      limits = circle_xlim,
      expand = c(0, 0)
    ) + scale_y_continuous(
      breaks = scales::breaks_width(5), labels = NULL,
      limits = circle_ylim,
      expand = c(0, 0)
    ) + scale_fill_manual(
      values = c(cluster_colors$germline, cluster_colors$somatic)
    ) + theme_bw() + theme(
      panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
      aspect.ratio = diff(circle_ylim) / diff(circle_xlim)
    ) + labs(
      x = NULL, y = NULL
    )
  attr(g, "n") <- list(
    n1 = length(gene_lists$germline),
    n2 = length(gene_lists$somatic),
    intersection = intersection
  )
  g
  # Common genes to annotate:
  # Myc, piwi, sgg, Top1, AGO1
  # Other common genes (add as many as possible)
  # robo1, Su(dx), kel, foi, stai
  # Zir, fwe
  # Somatic genes:
  # tj, roX2, Egfr, E(spl)m3, Timp, Wnt4, so
  # Germline genes:
  # AGO3, nos, vas, Dl, bam, blanks, peb, HP6, otu
  # stet, solo, htl, Appl, moon, Pxt
}