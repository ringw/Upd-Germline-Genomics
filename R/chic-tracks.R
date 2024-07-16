log2 <- function(n) (log(n) / log(2))

generate_chic_tracks <- function(
  chic.experiment.quantify,
  chic.experiment.quantify.smooth_bw40,
  chic.experiment.quantify.smooth_bw2000,
  granges_to_normalize
) {
  tribble(
    ~filename, ~score, ~score_smooth,
    "Rough_Input",
    list(elementMetadata(chic.experiment.quantify)[, 1]),
    list(elementMetadata(chic.experiment.quantify.smooth_bw2000)[, 1]),
    "Rough_Mark",
    list(elementMetadata(chic.experiment.quantify)[, 2]),
    list(elementMetadata(chic.experiment.quantify.smooth_bw2000)[, 2]),
    "Rough_Mark_L2FC",
    list(
      (
        elementMetadata(chic.experiment.quantify)[, 2]
        / elementMetadata(chic.experiment.quantify)[, 1]
      ) %>%
        `/`(enframe(.) %>% dplyr::slice(granges_to_normalize) %>% deframe %>% median(na.rm=T)) %>%
        pmax(2^-10) %>%
        pmin(2^10) %>%
        log2
    ),
    list(
      (
        elementMetadata(chic.experiment.quantify.smooth_bw2000)[, 2]
        / elementMetadata(chic.experiment.quantify.smooth_bw2000)[, 1]
      ) %>%
        `/`(enframe(.) %>% dplyr::slice(granges_to_normalize) %>% deframe %>% median(na.rm=T)) %>%
        pmax(2^-10) %>%
        pmin(2^10) %>%
        log2
    ),
    "Imputed_Input",
    list(elementMetadata(chic.experiment.quantify)[, 1]),
    list(elementMetadata(chic.experiment.quantify.smooth_bw2000)[, 1]),
    "Imputed_Mark",
    list(elementMetadata(chic.experiment.quantify)[, 2]),
    list(elementMetadata(chic.experiment.quantify.smooth_bw2000)[, 2]),
    "Imputed_Mark_L2FC",
    list(
      (
        elementMetadata(chic.experiment.quantify)[, 2]
        / elementMetadata(chic.experiment.quantify.smooth_bw40)[, 1]
      ) %>%
        `/`(enframe(.) %>% dplyr::slice(granges_to_normalize) %>% deframe %>% median(na.rm=T)) %>%
        pmax(2^-10) %>%
        pmin(2^10) %>%
        log2
    ),
    list(
      (
        elementMetadata(chic.experiment.quantify.smooth_bw2000)[, 2]
        / elementMetadata(chic.experiment.quantify.smooth_bw2000)[, 1]
      ) %>%
        `/`(enframe(.) %>% dplyr::slice(granges_to_normalize) %>% deframe %>% median(na.rm=T)) %>%
        pmax(2^-10) %>%
        pmin(2^10) %>%
        log2
    ),
    "FSeq_Input",
    list(elementMetadata(chic.experiment.quantify.smooth_bw40)[, 1]),
    list(elementMetadata(chic.experiment.quantify.smooth_bw2000)[, 1]),
    "FSeq_Mark",
    list(elementMetadata(chic.experiment.quantify.smooth_bw40)[, 2]),
    list(elementMetadata(chic.experiment.quantify.smooth_bw2000)[, 2]),
    "FSeq_Mark_L2FC",
    list(
      (
        elementMetadata(chic.experiment.quantify.smooth_bw40)[, 2]
        / elementMetadata(chic.experiment.quantify.smooth_bw40)[, 1]
      ) %>%
        `/`(enframe(.) %>% dplyr::slice(granges_to_normalize) %>% deframe %>% median(na.rm=T)) %>%
        pmax(2^-10) %>%
        pmin(2^10) %>%
        log2
    ),
    list(
      (
        elementMetadata(chic.experiment.quantify.smooth_bw2000)[, 2]
        / elementMetadata(chic.experiment.quantify.smooth_bw2000)[, 1]
      ) %>%
        `/`(enframe(.) %>% dplyr::slice(granges_to_normalize) %>% deframe %>% median(na.rm=T)) %>%
        pmax(2^-10) %>%
        pmin(2^10) %>%
        log2
    )
  )
}

generate_chic_tracks_peakcalling <- function(
  chic.experiment.quantify,
  chic.experiment.quantify.smooth_bw40,
  chic.experiment.quantify.smooth_bw2000,
  granges_to_normalize
) {
  if (is.null(chic.experiment.quantify))
    return(tibble(filename=character(0), score=list(), score_smooth=list()))
  p_peak_track <- chic.experiment.quantify$p_peak %>%
    replace(which(chic.experiment.quantify$L2FC < 0), 1)
  filename <- if (max(width(chic.experiment.quantify[1:10])) == 500) "Broad_PValue" else "Rough_PValue"
  tribble(
    ~filename, ~score, ~score_smooth,
    filename, list(p_peak_track), list(p_peak_track)
  )
}