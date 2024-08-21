autocorrelate_centered_granges <- function(gr, autocor_length=101) {
  square <- \(v) v^2
  data <- tibble(
    granges = as.list(split(gr, seqnames(gr))) %>%
      subset(sapply(., length) > 1),
    fft_length = 2^sapply(granges, \(gr) ceil_discrete_log(pmax(length(gr), 256))+1),
    autocor_est = mapply(
      \(gr, fft_length) c(gr$score, rep(0, fft_length - length(gr))) %>%
        fft %>%
        abs %>%
        square %>%
        fft(inv=TRUE) %>%
        Re %>%
        `/`(fft_length * norm(gr$score, type="2")^2) %>%
        `*`(length(gr)) %>%
        head(autocor_length),
      granges,
      fft_length,
      SIMPLIFY=FALSE
    )
  )
  setNames(
    purrr::reduce(data$autocor_est, `+`) / sum(sapply(data$granges, length)),
    seq(0, autocor_length - 1) * width(gr[10])
  )
}