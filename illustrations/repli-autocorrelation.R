tar_load(chic.tile.diameter_1000_chr)
tar_load(
  c(
    repli.timing_Germline_chr,
    repli.timing_Somatic_chr,
    repli.timing_Kc167_chr,
    repli.timing_S2_chr
  )
)
tk <- split(seq_along(chic.tile.diameter_1000_chr), seqnames(chic.tile.diameter_1000_chr)) %>%
  head(7) %>%
  sapply(\(v) c(v, rep(NA, 10000))) %>%
  do.call(c, .)
Y <- sapply(
  list(
    Germline = repli.timing_Germline_chr,
    Somatic = repli.timing_Somatic_chr,
    Kc167 = repli.timing_Kc167_chr,
    S2 = repli.timing_S2_chr
  ),
  \(obj) obj$score[1:max(tk, na.rm=T)]
) %>%
  scale() %>%
  `[`(value = tk, ) %>%
  replace(is.na(.), 0)
Z <- Y %>%
  mvfft() %>%
  `*`(Conj(.)) %>%
  mvfft(inverse=T) %>%
  Re() %>%
  head(5000) %>%
  `/`(rep(.[1,, drop=T], nrow(.)))

library(ggplot2)
apply(
  Z,
  2,
  \(v) tibble(
    x = seq(-4999, 4999) * 1000,
    y = c(rev(v[-1]), v)
  ),
  simplify=FALSE
) %>%
  bind_rows(.id = "celltype") %>%
  ggplot(aes(x, y, color=celltype)) +
  geom_line() +
  scale_x_continuous(
    labels = \(v) str_glue("{v/1000} kb")
  ) +
  coord_cartesian(
    100 * 1000 * c(-1, 1),
    c(0.65, 1),
    expand=F
  ) +
  theme(
    aspect.ratio = 1
  )

apply(
  Z,
  2,
  \(v) tibble(
    x = seq(-9, 9) * 1000,
    y = c(rev(v[2:10]), v[1:10])
  ),
  simplify=FALSE
) %>%
  bind_rows(.id = "celltype") %>%
  ggplot(aes(x, y, color=celltype)) +
  geom_line() +
  scale_x_continuous(
    labels = \(v) str_glue("{v/1000} kb")
  ) +
  theme(
    aspect.ratio = 1
  )

data <- apply(
  Z,
  2,
  \(v) tibble(
    x = seq(-4999, 4999) * 1000,
    y = c(rev(v[-1]), v)
  ),
  simplify=FALSE
) %>%
  bind_rows(.id = "celltype")

data <- data %>%
  subset(between(x, -500 * 1000, 500 * 1000))
data$celltype <- data$celltype %>%
  factor(c("Germline", "Somatic", "Kc167", "S2"))
line_cmap <- c(chic_line_track_colors$germline, chic_line_track_colors$somatic, hcl(50, 102, 57), hcl(255, 102, 57))
labs <- labs(
  x = "Lagged Distance",
  y = "Autocorrelation R",
  color = "Cell Type"
)
data %>%
  ggplot(aes(x, y, color=celltype)) +
  geom_line(
    linewidth = 1 * 25.4 / 72,
  ) +
  annotate(
    "rect",
    xmin = -25000,
    xmax = 25000,
    ymin = 0.75 - 0.25 * 0.05,
    ymax = 1 + 0.25 * 0.05,
    linewidth = 0.5 * 25.4 / 72,
    fill = "transparent",
    color = "black"
  ) +
  scale_x_continuous(
    labels = \(v) str_glue("{v/1000} kb")
  ) +
  scale_color_manual(values = line_cmap) +
  coord_cartesian(
    500 * 1000 * c(-1, 1),
    range(data$y) + c(-0.05, 0.05) * diff(range(data$y)),
    expand=F
  ) +
  labs +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = rel(1)),
    legend.position = "none",
  )
lft <- ggplotGrob(last_plot())
data %>%
  ggplot(aes(x, y, color=celltype)) +
  geom_line(
    linewidth = 1 * 25.4 / 72,
  ) +
  geom_point(
    size = 0.75
  ) +
  scale_x_continuous(
    labels = \(v) str_glue("{v/1000} kb")
  ) +
  scale_color_manual(values = line_cmap) +
  coord_cartesian(
    25 * 1000 * c(-1, 1),
    c(0.75, 1) + 0.25 * 0.05 * c(-1, 1),
    expand=F
  ) +
  labs +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = rel(1)),
  )
rgt <- ggplotGrob(last_plot())

plot(cbind(lft, rgt))
ggsave("figure/Both-Cell-Types/Repli-Autocorrelation-R.pdf", cbind(lft, rgt), w=8, h=3.5)
