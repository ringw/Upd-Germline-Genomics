cbind_features_at_tss <- function(
  chic_raws,
  gene_list,
  metafeatures_path,
  n_obs = 5120,
  tss_region = 4096
) {
  chic_raws <- chic_raws %>%
    split(cut(seq_along(chic_raws), c(0, quantile(seq_along(chic_raws), c(0.34, 0.67, 1)))))
  obs_take <- seq(n_obs) %>%
    split(cut(., seq(0, n_obs, length.out=4))) %>%
    sapply(length)
  obs_genes <- tibble(
    raw_samples = chic_raws,
    obs_take = obs_take
  ) %>%
    rowwise %>%
    reframe(
      raw_samples = list(raw_samples),
      genes_take = sample(gene_list, obs_take) %>%
        split(
          cut(seq_along(.), seq(0, length(.), by=128))
        )
    )
  do.call(
    cbind,
    mapply(
      \(raws, obs_take) cbind_features_at_tss_helper(
        raws, obs_take, metafeatures_path,
        n_obs = length(obs_take), tss_region = tss_region
      ),
      obs_genes$raw_samples,
      obs_genes$genes_take,
      SIMPLIFY=FALSE
    )
  )
}

cbind_features_at_tss_helper <- function(
  chic_raws,
  gene_list,
  metafeatures_path,
  n_obs = 512,
  tss_region = 4096
) {
  features <- metafeatures_path %>%
    read.csv(row.names = 1) %>%
    rownames_to_column %>%
    subset(chr %in% names(chr.lengths) & rowname %in% gene_list) %>%
    mutate(
      read_from = ifelse(strand == "+", start - tss_region/2, end + tss_region/2),
      read_by = ifelse(strand == "+", 1, -1)
    )
  chr_features <- features %>% split(.$chr)
  chr.starts <- cumsum(c(0, chr.lengths[-length(chr.lengths)])) %>%
    setNames(names(chr.lengths))
  chr_lut <- chr_features %>%
    mapply(
      \(df, chr_len, chr_start) (
        rep(
          df$read_from, each = tss_region
        ) + rep(
          df$read_by, each = tss_region
        ) * seq(
          0, tss_region - 1
        )
      ) %>%
        replace(
          !between(., 1, chr_len),
          NA
        ) %>%
        `+`(chr_start) %>%
        # "Extract" base pairs at position 1 in the chromsoome, as we generally
        # did not quantify any fragment ends at the first base pair. This avoids
        # an expensive replace() step on the big matrix.
        replace(is.na(.), 1),
      .,
      chr.lengths[names(.)],
      chr.starts[names(.)],
      SIMPLIFY=FALSE
    )
  chr_data <- cross_join(
    tibble(raw_vec = chic_raws, sample_name = names(chic_raws)),
    tibble(gene_lut = chr_lut, gene_list = chr_features %>% sapply(\(df) df$rowname, simplify=FALSE))
  ) %>%
    mutate(
      obs_counter = split(
        gene_list %>% sapply(length) %>% sum %>% seq,
        Rle(seq(length(gene_list)), gene_list %>% sapply(length))
      ) %>%
        as.list
    )
  obs_take <- sample(
    chr_data$gene_list %>% sapply(length) %>% sum,
    n_obs
  )
  chr_data <- chr_data %>%
    rowwise %>%
    mutate(
      gene_lut = list(
        gene_lut[
          rep(
            obs_counter %in% obs_take,
            each = tss_region
          )
        ]
      ),
      gene_list = list(gene_list[obs_counter %in% obs_take]),
      obs_counter = list(obs_counter[obs_counter %in% obs_take])
    )
  obs_mat <- chr_data %>%
    rowwise %>%
    filter(length(gene_list) > 0) %>%
    mutate(
      obs = matrix(
        raw_vec[gene_lut],
        ncol = length(gene_list),
        dimnames = list(NULL, paste0(sample_name, "_", gene_list))
      ) %>%
        list
    ) %>%
    pull(obs) %>%
    do.call(cbind, .)
}

psd_centered_features <- function(
  obs_mat,
  eval_region = seq(-1500, 1500, by=10)
) {
  # Size of each region centered around the tss
  tss_region <- nrow(obs_mat)

  num_hann <- 1025
  H <- hann(num_hann)
  obs_window_lookup <- cross_join(
    tibble(pos = seq(tss_region)),
    tibble(origin = tss_region / 2 + eval_region)
  ) %>%
    mutate(
      # Look into hann() according to distance of pos from the origin
      I = (pos + ceiling(num_hann / 2) - origin) %>%
        replace(!between(., 1, num_hann), NA)
    ) %>%
    # Grab hann value unless OOB (0)
    mutate(x = H[I] %>% replace(is.na(.), 0)) %>%
    pull(x) %>%
    # Construct a matrix of the sliding window. We take one column and multiply
    # element-wise the observation to window the observation.
    matrix(nrow = tss_region, ncol = length(eval_region), byrow=TRUE)

  square <- \(x) x^2
  psd_data <- tibble(
    x = eval_region,
    obs_lookup = seq_along(eval_region)
  ) %>%
    rowwise %>%
    mutate(
      PSD = (obs_mat * obs_window_lookup[, obs_lookup]) %>%
        mvfft %>%
        Mod %>%
        square %>%
        rowSums %>%
        list,
      .keep = "unused"
    ) %>%
    reframe(
      x,
      fourier_coef = seq(0, length(PSD) - 1),
      wavelength = c(Inf, tss_region / seq(1, length(PSD) - 1)),
      PSD
    )
}

create_psd_gradient <- function(...) scale_fill_gradientn(
  colors = c("white", viridis(7, option = "magma", end = 0.9, direction = -1)),
  values = c(-0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4) %>% rescale,
  ...
)

heatmap_psd_centered_features <- function(
  psd_data, ys = seq(50, 400, by=10)
) {
  tss_region <- max(psd_data$fourier_coef) + 1
  plot_data <- psd_data %>%
    group_by(x) %>%
    reframe(
      wavelength = ys,
      PSD = approx(fourier_coef, PSD, xout=tss_region / ys)$y,
      total_PSD = sum(PSD)
    ) %>%
    group_by(x) %>%
    mutate(PSD = PSD / max(total_PSD))
  ggplot(
    plot_data, aes(x, wavelength, fill=PSD)
  ) + geom_tile() + geom_segment(
    aes(xend=xend, yend=yend, fill=NULL),
    tibble(x=0, wavelength=-Inf, xend=0, yend=Inf),
    color = "darkred"
  ) + create_psd_gradient(
    labels = percent,
    limits = c(0, 0.1),
    oob = squish
  ) + coord_cartesian(
    NULL, c(50, 300),
    expand = FALSE
  ) + scale_x_continuous(
    breaks = c(-1000, -500, 0, 500, 1000),
    labels = c("-1000", "-500", "TSS", "+500", "+1000")
  ) + scale_y_continuous(
    breaks = seq(50, 300, by=25),
    labels = c(50, "", 100, "", 150, "", 200, "", 250, "", 300)
  ) + labs(
    x = "base pairs (center of periodic signal)",
    y = "Periodicity (bp)"
  ) + theme_bw()
}