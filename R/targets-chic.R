library(dplyr)
chic.samples = read.csv('chic/chic_samples.csv') %>%
  subset(sample != "" & !sapply(rejected, isTRUE))

chic.fpkm.data <- tribble(
  ~name, ~contrast, ~driver,
  'Germline', c(1,0,0,0,0,0,0), 'nos',
  'Somatic', c(1,1,0,0,0,0,0), 'tj'
)
# We are going to pull a "bed" target in the cross product below, which defines
# the genomic features shown in a heatmap for this particular experiment.
chic.fpkm.data$bed_sym <- rlang::syms(paste0("bed_", chic.fpkm.data$name))

chic.mark.data = tribble(~mark, 'H3K4', 'H3K27', 'H3K9')

# ChIC experiments: For every driver x mark.
chic.experiments <- chic.fpkm.data %>% cross_join(chic.mark.data)
# Pull the "chic" (track) target by name for the driver x mark targets.
chic.experiments <- chic.experiments %>% within({
  experiment_name <- paste(mark, name, sep="_")
  chic_input_sym <- rlang::syms(paste("chic.smooth_250", "input", mark, name, sep="_"))
  chic_mod_sym <- rlang::syms(paste("chic.smooth_25", "mod", mark, name, sep="_"))
  chic_smooth_input_sym <- rlang::syms(paste("chic.smooth", "input", mark, name, sep="_"))
  chic_smooth_250_input_sym <- rlang::syms(paste("chic.smooth_250", "input", mark, name, sep="_"))
  chic_smooth_mod_sym <- rlang::syms(paste("chic.smooth", "mod", mark, name, sep="_"))
  chic_smooth_125_mod_sym <- rlang::syms(paste("chic.smooth_125", "mod", mark, name, sep="_"))
  chic_smooth_250_mod_sym <- rlang::syms(paste("chic.smooth_250", "mod", mark, name, sep="_"))
  chic_input_bam <- rlang::syms(paste("chic.merge.bam", "input", mark, name, sep="_"))
  chic_mod_bam <- rlang::syms(paste("chic.merge.bam", "mod", mark, name, sep="_"))
})

target_chic_load_raw_expr <- function(group, driver) {
  group. <- group
  driver. <- driver
  samples <- chic.samples %>%
    filter(group == group., driver == driver.) %>%
    reframe(molecule, rep, sample, pileup = rlang::syms(str_glue("chic.raw_{sample}")))
  rlang::call2("tibble", rowname = samples$sample, molecule = samples$molecule, rep = call("factor", samples$rep), pileup = samples$pileup)
}

granges_to_normalize_euchromatin <- substitute(
  which(as.logical(seqnames(chic.tile.diameter_40_score) %in% c("2L", "2R", "3L", "3R", "4")))
)
# granges_to_normalize_heterochromatin <- substitute(
#   rep(
#     grep("[23]Cen", as.character(seqnames(chic.tile.diameter_40_score))),
#     seqlengths(seqinfo(chic.tile.diameter_40_score)[grep("[23]Cen", as.character(seqnames(chic.tile.diameter_40_score)), val=T)])
#   )
# )
granges_to_normalize_heterochromatin <- substitute(
  which(
    as.logical(
      seqnames(chic.tile.diameter_40_score) %in%
        grep("[23]Cen", names(chr.lengths), value=T)
    )
  )
)

targets.chic.aligned <- tar_map(
  bowtie.refs,
  names = name,
  tar_map(
    chic.samples,
    names = sample,
    unlist = FALSE,
    tar_file(
      chic.bam,
      with(
        list(name=name, group=group, sample=sample) %>%
          with(
            list(output_path = str_glue("chic/{name}/{group}/{sample}.bam"), batch=batch, sample=sample)
          ),
        {
          run(
            "bash",
            c(
              "-i",
              align_chic_lightfiltering,
              str_replace(bowtie[1], "\\..*", ""),
              paste0(batch, "/", sample, "_R1_001.fastq.gz"),
              paste0(batch, "/", sample, "_R2_001.fastq.gz"),
              output_path
            )
          )
          output_path
        }
      ),
      cue = tar_cue("never"),
      packages = c(
        "dplyr",
        "processx",
        "stringr"
      )
    )
  )
)

targets.chic <- list(
  targets.chic.aligned,

  tar_map(
    mutate(
      bowtie.refs,
      idxstats = rlang::syms(str_glue("bulk_reads_idxstats_chic.bam_GC3772016_S1_L002_{name}"))
    ),
    names = name,
    tar_target(
      chic.tile_chr_granges_list,
      GRanges(idxstats$rname %>% setdiff("*") %>% factor(., .), IRanges(1, width = idxstats$rlength[idxstats$rname != "*"])) %>%
        split(seqnames(.))
    ),
    tar_target(
      chic.tile_all_chr_diameter_20_granges_list,
      slidingWindows(
        unlist(chic.tile_chr_granges_list), w=20L, s=20L
      )
    ),
    tar_target(
      chic.tile_all_chr_diameter_100_granges_list,
      unlist(chic.tile_chr_granges_list, use.names=F) %>%
        slidingWindows(w=100L, s=100L) %>%
        setNames(names(chic.tile_chr_granges_list))
    ),
    tar_target(
      chic.tile_all_chr_diameter_1000_granges_list,
      slidingWindows(
        unlist(chic.tile_chr_granges_list), w=1000L, s=1000L
      )
    ),
    tar_target(
      chic.tile_all_chr_diameter_40_granges_list,
      # Windows should be centered on the diameter-20 windows:
      # First 1-40, then 11-50, 31-70, ...
      tibble(
        granges = slidingWindows(
          unlist(chic.tile_chr_granges_list), w=40L, s=10L
        ) %>%
          sapply(
            \(gr) gr[c(1, seq(2, length(gr), by=2))] %>%
              setNames(NULL) %>%
              append(
                str_glue(
                  as.character(seqnames(gr)[1]),
                  ":",
                  gr[length(gr) - (length(gr) %% 2)]@ranges@start + 20,
                  "-",
                  idxstats$rlength[idxstats$rname == as.character(seqnames(gr)[1])]
                )
              )
          ),
        testgranges = stopifnot(all.equal(sapply(granges, length), sapply(chic.tile_all_chr_diameter_20_granges_list, length)))
      ) %>%
        pull(granges)
    ),
    # Fix to the all chr diameter 100 track! Use resize/restrict which is much
    # simpler than new slidingWindows call/slice into the result to get the
    # desired step size.
    tar_target(
      chic.tile_all_chr_diameter_500_granges_list,
      chic.tile_all_chr_diameter_100_granges_list %>%
        sapply(\(gr) gr %>% resize(500L, "center") %>% restrict(1L, max(end(gr))))
    ),
    tar_target(
      chic.tile.diameter_40,
      sapply(
        levels(seqnames(chic.tile_chr_granges_list[[1]])),
        \(n) if (n %in% names(masked.feature.lengths))
          chic.tile_all_chr_diameter_40_granges_list[[n]]
        else
          chic.tile_chr_granges_list[[n]],
        USE.NAMES = FALSE
      ) %>%
        GRangesList %>%
        unlist %>%
        GRanges(
          seqlengths = idxstats %>%
            subset(
              rname != "*",
              select = c(rname, rlength)
            ) %>%
            deframe
        )
    ),
    tar_target(
      chic.tile.diameter_40_score,
      sapply(
        levels(seqnames(chic.tile_chr_granges_list[[1]])),
        \(n) if (n %in% names(masked.feature.lengths))
          chic.tile_all_chr_diameter_20_granges_list[[n]]
        else
          chic.tile_chr_granges_list[[n]],
        USE.NAMES = FALSE
      ) %>%
        GRangesList %>%
        unlist %>%
        GRanges(
          seqlengths = idxstats %>%
            subset(
              rname != "*",
              select = c(rname, rlength)
            ) %>%
            deframe
        )
    ),
    tar_target(
      chic.tile.diameter_500,
      sapply(
        levels(seqnames(chic.tile_chr_granges_list[[1]])),
        \(n) if (n %in% names(masked.feature.lengths))
          chic.tile_all_chr_diameter_500_granges_list[[n]]
        else
          chic.tile_chr_granges_list[[n]],
        USE.NAMES = FALSE
      ) %>%
        GRangesList %>%
        unlist %>%
        GRanges(
          seqlengths = idxstats %>%
            subset(
              rname != "*",
              select = c(rname, rlength)
            ) %>%
            deframe
        )
    ),
    tar_target(
      chic.tile.diameter_500_score,
      sapply(
        levels(seqnames(chic.tile_chr_granges_list[[1]])),
        \(n) if (n %in% names(masked.feature.lengths))
          chic.tile_all_chr_diameter_100_granges_list[[n]]
        else
          chic.tile_chr_granges_list[[n]],
        USE.NAMES = FALSE
      ) %>%
        GRangesList %>%
        unlist %>%
        GRanges(
          seqlengths = idxstats %>%
            subset(
              rname != "*",
              select = c(rname, rlength)
            ) %>%
            deframe
        )
    ),
    tar_target(
      chic.tile.diameter_1000,
      sapply(
        levels(seqnames(chic.tile_chr_granges_list[[1]])),
        \(n) if (n %in% names(masked.feature.lengths))
          chic.tile_all_chr_diameter_1000_granges_list[[n]]
        else
          chic.tile_chr_granges_list[[n]],
        USE.NAMES = FALSE
      ) %>%
        GRangesList %>%
        unlist %>%
        GRanges(
          seqlengths = idxstats %>%
            subset(
              rname != "*",
              select = c(rname, rlength)
            ) %>%
            deframe
        )
    ),
    tar_target(
      chic.tile.diameter_1000_lookup,
      tibble(
        rname = as.factor(seqnames(chic.tile.diameter_40)),
        ind = seq_along(chic.tile.diameter_40),
        loc = chic.tile.diameter_40@ranges@start + 1/2 * chic.tile.diameter_40@ranges@width
      ) %>%
        group_by(rname) %>%
        reframe(
          dest = list(chic.tile.diameter_1000[as.logical(seqnames(chic.tile.diameter_1000) == rname[1])]),
          lookup = ind[
            findInterval(
              dest[[1]]@ranges@start + 1/2 * dest[[1]]@ranges@width,
              loc
            ) %>%
              pmin(length(loc))
          ]
        ) %>%
        pull(lookup)
    )
  ),

  # ChIC loaded BAM files to KDE and the bandwidth-equivalent rectangular windows.
  tar_map(
    cross_join(
      bowtie.refs,
      chic.samples %>% subset(str_starts(molecule, "H3"))
    ) %>%
      cross_join(
        tribble(
          ~bp_suffix, ~pos_fixer_callable,
          "CN", rlang::sym("paired_end_pos_to_midpoint") # ,
          # "FE", rlang::sym("paired_end_pos_to_5_prime")
        )
      ) %>%
      rowwise %>%
      mutate(
        suffix = str_glue("{sample}_{bp_suffix}_{name}"),
        bulk_reads_misc = rlang::syms(str_glue("bulk_reads_misc_chic.bam_{sample}_{name}")),
        bulk_reads_idxstats = rlang::syms(str_glue("bulk_reads_idxstats_chic.bam_{sample}_{name}")),
        bulk_reads_split = mapply(
          \(ref_name, n) rlang::syms(
            str_glue(
              "bulk_reads_{names(if (name == \"masked\") masked.lengths else chr.lengths)}_{n}"
            )
          ) %>%
            setNames(names(if (name == "masked") masked.lengths else chr.lengths)) %>%
            append(list("list"), .) %>%
            do.call(call, ., quote=T),
          name,
          str_glue("chic.bam_{sample}_{name}"),
          SIMPLIFY=F
        ),
        chic.tile.diameter_40 = rlang::syms(str_glue("chic.tile.diameter_40_{name}")),
        chic.tile.diameter_500 = rlang::syms(str_glue("chic.tile.diameter_500_{name}")),
        max_fragment_length = c(`20240628`=500)[as.character(rep)] %>%
          replace(is.na(.), 200)
      ),
    names = suffix,

    tar_target(
      chic.granges.diameter_40,
      tibble(
        reads = append(bulk_reads_split, list(bulk_reads_misc)) %>%
          sapply(
            \(df) df %>%
              pos_fixer_callable %>%
              filter(between(length, 100, max_fragment_length), mapq >= 20),
            simplify=F
          ) %>%
          setNames(NULL),
        tile_granges = reads %>%
          sapply(
            \(df) chic.tile.diameter_40[
              if (!nrow(df))
                integer(0)
              else if (df$rname[1] %in% names(masked.lengths))
                seqnames(chic.tile.diameter_40) == df$rname[1]
              else
                !(seqnames(chic.tile.diameter_40) %in% names(masked.lengths))
            ]
          ),
        partition_tiles = stopifnot(sum(sapply(tile_granges, length)) == length(chic.tile.diameter_40)),
        granges = mapply(
          count_overlaps_bases, tile_granges, reads
        )
      ) %>%
        with(
          granges %>%
            GRangesList %>%
            unlist %>%
            `metadata<-`(value = list(est_library_size = sum(sapply(reads, nrow))))
        )
    ),
    tar_target(
      chic.granges.peakcalling.diameter_40,
      tibble(
        reads = append(bulk_reads_split, list(bulk_reads_misc)) %>%
          sapply(
            \(df) df %>%
              paired_end_reads_to_fragment_lengths %>%
              filter(between(length, 100, max_fragment_length), mapq >= 20),
            simplify=F
          ) %>%
          setNames(NULL),
        tile_granges = reads %>%
          sapply(
            \(df) chic.tile.diameter_40[
              if (!nrow(df))
                integer(0)
              else if (df$rname[1] %in% names(masked.lengths))
                seqnames(chic.tile.diameter_40) == df$rname[1]
              else
                !(seqnames(chic.tile.diameter_40) %in% names(masked.lengths))
            ]
          ),
        partition_tiles = stopifnot(sum(sapply(tile_granges, length)) == length(chic.tile.diameter_40)),
        granges = mapply(
          \(tile_granges, reads) GRanges(
            tile_granges,
            score = countOverlaps(
              tile_granges,
              reads %>% paired_end_reads_to_granges
            )
          ),
          tile_granges,
          reads
        )
      ) %>%
        with(
          granges %>%
            GRangesList %>%
            unlist %>%
            `metadata<-`(value = list(est_library_size = sum(sapply(reads, nrow))))
        )
    ),
    tar_target(
      chic.granges.peakcalling.diameter_500,
      tibble(
        reads = append(bulk_reads_split, list(bulk_reads_misc)) %>%
          sapply(
            \(df) df %>%
              paired_end_reads_to_fragment_lengths %>%
              filter(between(length, 100, max_fragment_length), mapq >= 20),
            simplify=F
          ) %>%
          setNames(NULL),
        tile_granges = reads %>%
          sapply(
            \(df) chic.tile.diameter_500[
              if (!nrow(df))
                integer(0)
              else if (df$rname[1] %in% names(masked.lengths))
                seqnames(chic.tile.diameter_500) == df$rname[1]
              else
                !(seqnames(chic.tile.diameter_500) %in% names(masked.lengths))
            ]
          ),
        partition_tiles = stopifnot(sum(sapply(tile_granges, length)) == length(chic.tile.diameter_500)),
        granges = mapply(
          \(tile_granges, reads) GRanges(
            tile_granges,
            score = countOverlaps(
              tile_granges,
              reads %>% paired_end_reads_to_granges
            )
          ),
          tile_granges,
          reads
        )
      ) %>%
        with(
          granges %>%
            GRangesList %>%
            unlist %>%
            `metadata<-`(value = list(est_library_size = sum(sapply(reads, nrow))))
        )
    )
  ),

  tar_map(
    chic.experiments %>%
      cross_join(dplyr::rename(bowtie.refs, reference="name")) %>%
      cross_join(tibble(bp_name = c("CN", "peakcalling.sharp", "peakcalling.broad"))) %>%
      mutate(driver = driver) %>%
      rowwise %>%
      mutate(
        chic.tile.diameter_40_score = rlang::syms(str_glue("chic.tile.diameter_40_score_{reference}")),
        chic_sample_names = subset(
          chic.samples$sample,
          chic.samples$driver == driver &
          chic.samples$group == mark &
          str_starts(chic.samples$molecule, "H3")
        ) %>%
          list,
        bp_suffix = bp_name %>% replace(which(grepl("peakcalling", .)), "CN"),
        chic.granges = rlang::syms(
          str_glue(
            "chic.granges",
            if (str_starts(bp_name, "peakcalling")) ".peakcalling" else "",
            ".diameter_{c(CN=40, peakcalling.sharp=40, peakcalling.broad=500)[bp_name]}_{chic_sample_names}_{bp_suffix}_{reference}"
          )
        ) %>%
          list,
        chic.tile = rlang::syms(
          str_glue("chic.tile.diameter_{c(CN=40, peakcalling.sharp=40, peakcalling.broad=500)[bp_name]}_{reference}")
        ),
        chic.tile.score = rlang::syms(
          str_glue("chic.tile.diameter_{c(CN=40, peakcalling.sharp=40, peakcalling.broad=500)[bp_name]}_score_{reference}")
        ),
        chic.tile.diameter_1000 = rlang::syms(str_glue("chic.tile.diameter_1000_{reference}")),
        chic.tile.diameter_1000_lookup = rlang::syms(str_glue("chic.tile.diameter_1000_lookup_{reference}")),
        experimental_group = (
          chic.samples$group[match(chic_sample_names, chic.samples$sample)]
        ) %>%
          list,
        molecule = (
          chic.samples$molecule[match(chic_sample_names, chic.samples$sample)]
        ) %>%
          list,
        model_matrix_rep = (
          chic.samples$rep[match(chic_sample_names, chic.samples$sample)]
        ) %>%
          list,
        chic.experiment.granges.offset.elementMetadata = rlang::syms(
          str_glue(
            "chic.experiment.granges.offset.",
            c(H3K4="euchromatin", H3K27="euchromatin", H3K9="euchromatin")[mark],
            "_{mark}_{name}_{bp_name}_{reference}"
          )
        ),
        granges_to_normalize = list(
          do.call(
            substitute,
            list(
              list(
                H3K4=granges_to_normalize_euchromatin,
                H3K27=granges_to_normalize_euchromatin,
                H3K9=granges_to_normalize_euchromatin
              )[[mark]],
              list(chic.tile.diameter_40_score = chic.tile.diameter_40_score)
            )
          )
        ),
        generate_chic_tracks = rlang::syms(
          c(
            CN="generate_chic_tracks",
            FE="generate_chic_tracks",
            peakcalling.sharp="generate_chic_tracks_peakcalling",
            peakcalling.broad="generate_chic_tracks_peakcalling"
          )
        )[bp_name],
        # global_overdispersion = bp_name != "peakcalling.broad",
        test_de = str_starts(bp_name, "peakcalling")
      ),
    names = mark | name | bp_name | reference,
    tar_target(
      chic.experiment.granges,
      chic.granges[[1]] %>%
        attributes %>%
        append(list(molecule=molecule, rep=model_matrix_rep)) %>%
        with(
          GRanges(
            seqnames,
            ranges,
            seqlengths = seqlengths(chic.granges[[1]]),
            score = sapply(chic.granges, \(gr) gr$score) %>%
              matrix(
                nrow = nrow(.),
                dimnames = list(NULL, str_glue("{molecule}_Rep{rep}"))
              )
          ) %>%
            `metadata<-`(
              value = list(
                est_library_size = sapply(
                  chic.granges,
                  \(gr) gr@metadata$est_library_size
                ) %>%
                  setNames(str_glue("{molecule}_Rep{rep}"))
              )
            )
        )
    ),
    tar_target(
      chic.experiment.granges.offset.euchromatin,
      # Median of 1kb-window coverage scores (rescaled to match the scale of
      # the 40 bp count overlaps). Because each fragment is assigned exactly
      # 1 bp to count overlaps with, we can exactly scale down by the
      # difference in window size using the rolling mean of windows.
      matrix(
        log(
          chic.experiment.granges[
            seqnames(chic.experiment.granges) %in%
            c("2L", "2R", "3L", "3R", "4")
          ] %>%
            elementMetadata %>%
            apply(
              2,
              \(v) v %>% rollmean(50, fill=0, align="center") %>% median
            )
        ) %>%
          rep(each = length(chic.experiment.granges))
        + log(ranges(chic.experiment.granges)@width)
        - log(
          median(
            ranges(chic.experiment.granges)@width[
              as.logical(
                seqnames(chic.experiment.granges) %in%
                c("2L", "2R", "3L", "3R", "4")
              )
            ]
          )
        ),
        nrow = length(chic.experiment.granges),
        ncol = ncol(elementMetadata(chic.experiment.granges)),
        dimnames = list(
          NULL,
          colnames(elementMetadata(chic.experiment.granges))
        )
      )
    ),
    tar_target(
      chic.experiment.granges.offset.heterochromatin,
      matrix(
        log(
          chic.experiment.granges[
            grepl("Cen", seqnames(chic.experiment.granges))
          ] %>%
            attributes %>%
            with(
              apply(
                elementMetadata,
                2,
                \(v) median(v / ranges@width)
              )
            )
        ) %>%
          rep(each = length(chic.experiment.granges))
        + log(ranges(chic.experiment.granges)@width),
        nrow = length(chic.experiment.granges),
        ncol = ncol(elementMetadata(chic.experiment.granges)),
        dimnames = list(
          NULL,
          colnames(elementMetadata(chic.experiment.granges))
        )
      )
    ),
    tar_target(
      chic.experiment.granges.offset,
      GRanges(
        seqnames(chic.experiment.granges),
        ranges(chic.experiment.granges),
        elementMetadata = chic.experiment.granges.offset.elementMetadata
      )
    ),
    tar_target(
      chic.experiment.quantify,
      chic_quantify(
        chic.experiment.granges,
        molecule,
        model_matrix_rep,
        offset = chic.experiment.granges.offset.elementMetadata,
        # global_overdispersion = global_overdispersion,
        test_de = test_de
      )
    ),

    tar_map(
      tibble(
        bw = c(40, 2000),
        bw_name = str_glue("bw{str_trim(format(bw, scientific=F))}")
      ),
      names = bw_name,
      tar_target(
        chic.experiment.quantify.smooth,
        chic.experiment.quantify %>%
          split(seqnames(.)) %>%
          sapply(\(gr) ksmooth_sliding_windows(gr, bw = bw)) %>%
          GRangesList %>%
          unlist
      )
    ),

    tar_file(
      chic.bw.tracks,
      generate_chic_tracks(
        chic.experiment.quantify = chic.experiment.quantify,
        chic.experiment.quantify.smooth_bw40 = chic.experiment.quantify.smooth_bw40,
        chic.experiment.quantify.smooth_bw2000 = chic.experiment.quantify.smooth_bw2000,
        granges_to_normalize = granges_to_normalize
      ) %>%
        rowwise %>%
        summarise(
          filename = with(
            list(mark=mark, name=name, bp_suffix=bp_suffix, reference=reference),
            str_glue("chic/{reference}/{name}_{mark}_{bp_suffix}_{filename}.bw")
          ),
          dir.create = dir.create(dirname(filename), rec=T, showW=F),
          gr = GRanges(
            seqnames(chic.tile.score),
            ranges(chic.tile.score),
            seqlengths = seqlengths(chic.tile.score),
            score = (
              if (grepl("FSeq_[^I]|Imputed", filename))
                score[[1]] %>%
                  replace(
                    which(elementMetadata(chic.experiment.quantify.smooth_bw40)[, 1] < 0.1),
                    score_smooth[[1]][
                      elementMetadata(chic.experiment.quantify.smooth_bw40)[, 1] < 0.1
                    ]
                  )
              else score[[1]]
            ) %>%
              # Now replace actual holes (e.g. in an FBte sequence) with 1 (no enrichment).
              replace(is.na(.), 0)
          ) %>%
            list,
          filename = export(gr, BigWigFile(filename)) %>% as.character
        ) %>%
        pull(filename),
      packages = tar_option_get("packages") %>% c("tidyr")
    ),
    tar_file(
      chic.bw.track.wide,
      if (!grepl("peakcalling", bp_name))
        GRanges(
          chic.tile.diameter_1000,
          score = (
            elementMetadata(chic.experiment.quantify.smooth_bw2000)[, 2]
            / elementMetadata(chic.experiment.quantify.smooth_bw2000)[, 1]
          ) %>%
            `/`(enframe(.) %>% dplyr::slice(granges_to_normalize) %>% deframe %>% median(na.rm=T)) %>%
            `[`(chic.tile.diameter_1000_lookup) %>%
            log2 %>%
            replace(is.na(.), 0)
        ) %>%
          export(
            with(
              list(mark=mark, name=name, bp_suffix=bp_suffix, reference=reference),
              str_glue("chic/{reference}/{name}_{mark}_{bp_suffix}_Wide_Mark_L2FC.bw")
            ) %>%
              BigWigFile
          ) %>%
          as.character
    ),
    tar_target(
      chic.heatmap.tss,
      if (!grepl("peakcalling", bp_name))
        track_to_heatmap(
          grep(
            str_glue(
              if (grepl("H3K9", chic.bw.tracks[1])) "FSeq" else "Imputed",
              "_Mark_L2FC"
            ), chic.bw.tracks, val=T
          ) %>% BigWigFile %>% import %>%
            attributes %>%
            with(
              GRanges(
                seqnames,
                ranges,
                seqinfo = seqinfo,
                score = exp(elementMetadata$score * log(2))
              )
            ),
          as_tibble(read.csv(assay.data.sc)),
          grep("FSeq_Input", chic.bw.tracks, val=T) %>% BigWigFile %>% import,
          mask_threshold = 0
        )
    ),
    tar_target(
      chic.heatmap.tss.nucleosome,
      if (!grepl("peakcalling", bp_name))
        track_to_heatmap(
          grep("Imputed_Input", chic.bw.tracks, val=T) %>% BigWigFile %>% import %>%
            attributes %>%
            with(
              GRanges(
                seqnames,
                ranges,
                seqinfo = seqinfo,
                score = exp(elementMetadata$score * log(2))
              )
            ),
          as_tibble(read.csv(assay.data.sc)),
          grep("FSeq_Input", chic.bw.tracks, val=T) %>% BigWigFile %>% import,
          mask_threshold = 0
        )
    ),
    tar_target(
      chic.heatmap.paneled,
      if (!grepl("peakcalling", bp_name))
        track_to_heatmap_with_panels(
          grep(
            str_glue(
              if (grepl("H3K9", chic.bw.tracks[1])) "FSeq" else "Imputed",
              "_Mark_L2FC"
            ), chic.bw.tracks, val=T
          ) %>% BigWigFile %>% import %>%
            attributes %>%
            with(
              GRanges(
                seqnames,
                ranges,
                seqinfo = seqinfo,
                score = exp(elementMetadata$score * log(2))
              )
            ),
          as_tibble(read.csv(assay.data.sc)),
          mask_track = grep("FSeq_Input", chic.bw.tracks, val=T) %>% BigWigFile %>% import,
          mask_threshold = 0
        )
    )
  ),

  tar_target(
    chic.experiment.nucleosomes,
    chic.tile.diameter_40_score_chr %>%
      `elementMetadata<-`(
        value = cbind(
          elementMetadata(chic.experiment.granges_H3K27_Germline_CN_chr) %>%
            subset(select=c(score.H3_Rep1, score.H3_Rep3, score.H3_Rep20240528)),
          elementMetadata(chic.experiment.granges_H3K27_Somatic_CN_chr) %>%
            subset(select=c(score.H3_Rep1, score.H3_Rep4, score.H3_Rep5))
        ) %>%
          `colnames<-`(
            value = as.character(interaction(rep(c("GSC", "CySC"), each=3), 1:3))
          )
      )
  ),
  tar_target(
    chic.experiment.nucleosomes.offset,
    cbind(
      chic.experiment.granges.offset.euchromatin_H3K27_Germline_CN_chr[
        , c("score.H3_Rep1", "score.H3_Rep3", "score.H3_Rep20240528")
      ],
      chic.experiment.granges.offset.euchromatin_H3K27_Somatic_CN_chr[
        , c("score.H3_Rep1", "score.H3_Rep4", "score.H3_Rep5")
      ]
    ) %>%
      `colnames<-`(value = colnames(chic.experiment.nucleosomes))
  ),
  tar_target(
    chic.experiment.quantify.nucleosomes,
    chic_quantify(
      chic.experiment.nucleosomes,
      structure(
        c(2L, 2L, 2L, 1L, 1L, 1L),
        levels=c("H3.CySC", "H3.GSC"),
        class = "factor"
      ),
      rep(1L, 6),
      chic.experiment.nucleosomes.offset,
      test_de = TRUE,
      test_wald = TRUE
    )
  ),
  tar_target(
    chic.experiment.nuccall.background,
    GRanges(
      chic.tile.diameter_40_score_chr,
      logMu_H3.CySC = log(chic.experiment.quantify_H3K27_Somatic_CN_chr$score.molH3),
      logMu_H3.GSC = log(chic.experiment.quantify_H3K27_Germline_CN_chr$score.molH3)
    ) %>%
      split(seqnames(.)) %>%
      sapply(
        \(gr) if (length(gr) >= 50)
          gr %>%
            `elementMetadata<-`(
              value = gr %>%
                elementMetadata %>%
                apply(
                  2,
                  \(v) v %>%
                    list %>%
                    mapply(
                      \(v, k) rollapply(v, k, median, fill=0, align="center"),
                      .,
                      c(1000/20, 5000/20, 10000/20)
                    ) %>%
                    rowMaxs %>%
                    pmax(0),
                  simplify=FALSE
                ) %>%
                do.call(tibble, .)
            )
        else gr
      ) %>%
      GRangesList %>%
      unlist
  ),
  tar_target(
    chic.test.nucleosomes,
    GRanges(
      chic.tile.diameter_40_score_chr,
      p_CySC_nuc = pnorm(
        log(chic.experiment.quantify.nucleosomes$score.molH3.CySC),
        chic.experiment.nuccall.background$logMu_H3.CySC,
        sd = chic.experiment.quantify.nucleosomes$se_H3.CySC,
        lower.tail = FALSE
      ),
      p_GSC_nuc = pnorm(
        log(chic.experiment.quantify.nucleosomes$score.molH3.GSC),
        chic.experiment.nuccall.background$logMu_H3.GSC,
        sd = chic.experiment.quantify.nucleosomes$se_H3.GSC,
        lower.tail = FALSE
      ),
      p_diff_two_tail = chic.experiment.quantify.nucleosomes$p_peak,
      L2FC = chic.experiment.quantify.nucleosomes$L2FC
    )
  ),
  tar_file(
    chic.test.nucleosomes.bw.tracks,
    tibble(
      name = c("Somatic_Nucleosome", "Germline_Nucleosome", "Germline_Nuc_Diff_Depleted", "Germline_Nuc_Diff_Enriched"),
      trck = list(
        chic.test.nucleosomes$p_CySC_nuc,
        chic.test.nucleosomes$p_GSC_nuc,
        pmin(1, 2 * chic.test.nucleosomes$p_diff_two_tail) %>%
          replace(which(chic.test.nucleosomes$p_CySC_nuc >= 0.05), NA) %>%
          replace(which(is.na(chic.test.nucleosomes$p_CySC_nuc)), NA),
        pmin(1, 2 * (1 - chic.test.nucleosomes$p_diff_two_tail)) %>%
          replace(which(chic.test.nucleosomes$p_GSC_nuc >= 0.05), NA) %>%
          replace(which(is.na(chic.test.nucleosomes$p_GSC_nuc)), NA)
      )
    ) %>%
      rowwise %>%
      reframe(
        write_raw = GRanges(chic.tile.diameter_40_score_chr, score = trck %>% replace_na(1)) %>%
          export(BigWigFile(str_glue("chic/chr/{name}_PValue.bw"))) %>%
          as.character,
        write_indicator = GRanges(chic.tile.diameter_40_score_chr, score = trck %>% `<=`(0.05) %>% as.numeric %>% replace_na(0)) %>%
          export(BigWigFile(str_glue("chic/chr/{name}_Location.bw"))) %>%
          as.character
      ) %>%
      unlist,
    packages = tar_option_get("packages") %>% c("tidyr")
  ),
  tar_target(
    chic.test.nucleosomes.fix_Germline,
    nucleosomes_cleanup_tracks(
      chic.tile.diameter_40_score_chr,
      chic.test.nucleosomes$p_GSC_nuc,
      pmin(1, 2 * (1 - chic.test.nucleosomes$p_diff_two_tail))
    ),
    packages = tar_option_get("packages") %>% c("tidyr")
  ),
  tar_target(
    chic.test.nucleosomes.fix_Somatic,
    nucleosomes_cleanup_tracks(
      chic.tile.diameter_40_score_chr,
      chic.test.nucleosomes$p_CySC_nuc,
      pmin(1, 2 * chic.test.nucleosomes$p_diff_two_tail)
    ),
    packages = tar_option_get("packages") %>% c("tidyr")
  ),
  tar_file(
    chic.test.nucleosomes.fix.bw.tracks,
    tibble(
      name = c("Somatic_Nucleosome", "Germline_Nucleosome", "Germline_Nuc_Diff_Depleted", "Germline_Nuc_Diff_Enriched"),
      gr = list(
        chic.test.nucleosomes.fix_Somatic$Nucleosomes,
        chic.test.nucleosomes.fix_Germline$Nucleosomes,
        chic.test.nucleosomes.fix_Somatic$Diff_Enriched,
        chic.test.nucleosomes.fix_Germline$Diff_Enriched
      )
    ) %>%
      rowwise %>%
      reframe(
        write_fix = gr %>%
          split(seqnames(gr)) %>%
          as.list %>%
          enframe %>%
          rowwise %>%
          reframe(
            value = coverage(value)[[name]] %>%
              attributes %>%
              with(
                if (length(lengths) == 0)
                  GRanges()
                else GRanges(
                  name,
                  IRanges(
                    cumsum(c(1, lengths[-length(lengths)])),
                    width = lengths
                  ),
                  score = as.numeric(values)
                )
              ) %>%
              list
          ) %>%
          unlist(use.names = F) %>%
          GRangesList %>%
          unlist(use.names = F) %>%
          GRanges(
            seqinfo = seqinfo(chic.tile.diameter_40_score_chr)
          ) %>%
          export(BigWigFile(str_glue("chic/chr/{name}_Fix.bw"))) %>%
          as.character
      ) %>%
      unlist
  ),

  # Binning of chic H3 fragments by length and by mapq.
  tar_target(
    chic.nucleosome.fragment.stats.cut.bounds,
    seq(50, 500, by=50)
  ),
  tar_target(
    chic.nucleosome.fragment.stats.cut.levels,
    paste0(
      chic.nucleosome.fragment.stats.cut.bounds %>% head(length(.) - 1),
      "..",
      chic.nucleosome.fragment.stats.cut.bounds %>% tail(length(.) - 1)
    )
  ),
  tar_target(
    chic.mapq.stats.cut.bounds,
    c(0, 1, seq(10, 50, by=10))
  ),
  tar_target(
    chic.mapq.stats.cut.levels,
    c(
      "0",
      "1..10",
      "10..20",
      "20..30",
      "30..40",
      "40..42"
    )
  ),
  tar_map(
    chic.samples %>%
      mutate(celltype = driver %>% forcats::fct_recode(Germline="nos", Somatic="tj") %>% as.character) %>%
      rowwise %>%
      reframe(
        celltype,
        group,
        molecule,
        bulk_reads = rlang::syms(str_glue("bulk_reads_{names(chr.lengths)}_chic.bam_{sample}_chr")),
        chr = names(chr.lengths)
      ) %>%
      filter(group == "H3K4", molecule == "H3") %>%
      group_by(celltype) %>%
      summarise(
        bulk_reads = list(bulk_reads),
        chr = list(chr)
      ),
    names = celltype,
    tar_target(
      chic.nucleosome.fragment.stats,
      left_join(
        cross_join(
          tibble(chr = factor(c("2", "3", "4", "X", "Y"))),
          tibble(fragment_sizes = chic.nucleosome.fragment.stats.cut.levels)
        ),
        sapply(
          bulk_reads,
          \(df) df %>%
            paired_end_reads_to_fragment_lengths %>%
            group_by(cut(length, chic.nucleosome.fragment.stats.cut.bounds)) %>%
            tally %>%
            `colnames<-`(value = c("fragment_sizes", "n")) %>%
            subset(!is.na(fragment_sizes)) %>%
            within(levels(fragment_sizes) <- chic.nucleosome.fragment.stats.cut.levels),
          simplify=F
        ) %>%
          setNames(fct_recode(chr, `2`="2L", `2`="2R", `3`="3L", `3`="3R")) %>%
          bind_rows(.id = "chr"),
        c("chr", "fragment_sizes")
      ) %>%
        mutate(n = n %>% replace_na(0)),
      packages = tar_option_get("packages") %>% c("tidyr")
    ),
    tar_target(
      chic.h3.mapq.stats,
      left_join(
        cross_join(
          tibble(chr = factor(c("2", "3", "4", "X", "Y"))),
          tibble(mapq_range = chic.mapq.stats.cut.levels)
        ),
        sapply(
          bulk_reads,
          \(df) df %>%
            group_by(cut(mapq, chic.mapq.stats.cut.bounds)) %>%
            tally %>%
            `colnames<-`(value = c("mapq_range", "n")) %>%
            subset(!is.na(mapq_range)) %>%
            within(levels(mapq_range) <- chic.mapq.stats.cut.levels),
          simplify=F
        ) %>%
          setNames(fct_recode(chr, `2`="2L", `2`="2R", `3`="3L", `3`="3R")) %>%
          bind_rows(.id = "chr"),
        c("chr", "mapq_range")
      ) %>%
        mutate(n = n %>% replace_na(0)),
      packages = tar_option_get("packages") %>% c("tidyr")
    )
  ),
  
  # Peak calling at TSS using regression of ChIC samples.
  tar_target(
    chic.gene.enrichment,
    tibble(
      read.csv(assay.data.sc),
      chr.test = chr %>% replace(!(. %in% names(chr.lengths)), "*"),
      windows.rough = findOverlaps(
        GRanges(
          chr.test,
          IRanges(
            ifelse(
              chr.test != "*",
              ifelse(
                strand == "+",
                start,
                end - 499
              ),
              1
            ),
            width=500
          ),
          seqlengths = c(seqlengths(chic.tile.diameter_500_score_chr), `*`=1)
        ),
        chic.tile.diameter_40_score_chr
      ) %>%
        as("List") %>%
        as.list,
      windows.broad = findOverlaps(
        GRanges(
          chr.test,
          IRanges(
            ifelse(
              chr.test != "*",
              ifelse(
                strand == "+",
                start + 250,
                end - 250
              ),
              1
            ),
            width=1
          ),
          seqlengths = c(seqlengths(chic.tile.diameter_500_score_chr), `*`=1)
        ),
        chic.tile.diameter_500_score_chr
      ) %>%
        as("List") %>%
        as.list,
      windows.rbind = mapply(
        \(v1, v2) c(v1, length(chic.tile.diameter_40_score_chr) + v2),
        windows.rough,
        windows.broad
      ),
      as_tibble(
        sapply(
          list(
            H3K4_Germline=c(
              chic.experiment.quantify_H3K4_Germline_peakcalling.sharp_chr,
              chic.experiment.quantify_H3K4_Germline_peakcalling.broad_chr
            ),
            H3K27_Germline=c(
              chic.experiment.quantify_H3K27_Germline_peakcalling.sharp_chr,
              chic.experiment.quantify_H3K27_Germline_peakcalling.broad_chr
            ),
            H3K9_Germline=c(
              chic.experiment.quantify_H3K9_Germline_peakcalling.sharp_chr,
              chic.experiment.quantify_H3K9_Germline_peakcalling.broad_chr
            ),
            H3K4_Somatic=c(
              chic.experiment.quantify_H3K4_Somatic_peakcalling.sharp_chr,
              chic.experiment.quantify_H3K4_Somatic_peakcalling.broad_chr
            ),
            H3K27_Somatic=c(
              chic.experiment.quantify_H3K27_Somatic_peakcalling.sharp_chr,
              chic.experiment.quantify_H3K27_Somatic_peakcalling.broad_chr
            ),
            H3K9_Somatic=c(
              chic.experiment.quantify_H3K9_Somatic_peakcalling.sharp_chr,
              chic.experiment.quantify_H3K9_Somatic_peakcalling.broad_chr
            )
          ),
          \(gr) sapply(
            windows.rbind,
            \(w) with(
              elementMetadata(gr)[w, ],
              min(
                p_peak[L2FC > 0],
                if (length(p_peak) > 0) 1 else numeric(0),
                na.rm=T
              ) %>%
                replace(!is.finite(.), NA)
            )
          ),
          simplify=FALSE
        )
      )
    ) %>%
      reframe(
        symbol = X,
        flybase,
        chr,
        start,
        end,
        strand,
        H3K4_Germline,
        H3K27_Germline,
        H3K9_Germline,
        H3K4_Somatic,
        H3K27_Somatic,
        H3K9_Somatic
      ),
    format = "parquet"
  ),
  tar_target(
    chic.gene.enrichment.broad,
    tibble(
      read.csv(assay.data.sc),
      chr.test = chr %>%
        replace(!(. %in% names(chr.lengths)) | is.na(start) | is.na(end), "*"),
      windows.broad = findOverlaps(
        GRanges(
          chr.test,
          IRanges(
            start=start %>% replace(is.na(.), 1),
            end=end %>% replace(is.na(.), 1)
          ),
          seqlengths = c(seqlengths(chic.tile.diameter_500_score_chr), `*`=1)
        ),
        chic.tile.diameter_500_score_chr
      ) %>%
        as("List") %>%
        as.list,
      as_tibble(
        sapply(
          list(
            H3K4_Germline=chic.experiment.quantify_H3K4_Germline_peakcalling.broad_chr,
            H3K27_Germline=chic.experiment.quantify_H3K27_Germline_peakcalling.broad_chr,
            H3K9_Germline=chic.experiment.quantify_H3K9_Germline_peakcalling.broad_chr,
            H3K4_Somatic=chic.experiment.quantify_H3K4_Somatic_peakcalling.broad_chr,
            H3K27_Somatic=chic.experiment.quantify_H3K27_Somatic_peakcalling.broad_chr,
            H3K9_Somatic=chic.experiment.quantify_H3K9_Somatic_peakcalling.broad_chr
          ),
          \(gr) sapply(
            windows.broad,
            \(w) with(
              elementMetadata(gr)[w, ],
              mean(
                p_peak < 0.001 & L2FC > 0,
                na.rm=T
              ) %>%
                replace(!is.finite(.), NA)
            )
          ),
          simplify=FALSE
        )
      )
    ) %>%
      reframe(
        symbol = X,
        flybase,
        chr,
        start,
        end,
        strand,
        H3K4_Germline,
        H3K27_Germline,
        H3K9_Germline,
        H3K4_Somatic,
        H3K27_Somatic,
        H3K9_Somatic
      ),
    format = "parquet"
  ),

  tar_file(
    sd_chic_fragments,
    publish_chic_fragments(
      list(Germline=chic.nucleosome.fragment.stats_Germline, Somatic=chic.nucleosome.fragment.stats_Somatic),
      list(Germline=chic.h3.mapq.stats_Germline, Somatic=chic.h3.mapq.stats_Somatic),
      peaks_publish = tibble(
        chic.gene.enrichment,
        H3K4_Germline_Broad = chic.gene.enrichment.broad$H3K4_Germline,
        H3K27_Germline_Broad = chic.gene.enrichment.broad$H3K27_Germline,
        H3K9_Germline_Broad = chic.gene.enrichment.broad$H3K9_Germline,
        H3K4_Somatic_Broad = chic.gene.enrichment.broad$H3K4_Somatic,
        H3K27_Somatic_Broad = chic.gene.enrichment.broad$H3K27_Somatic,
        H3K9_Somatic_Broad = chic.gene.enrichment.broad$H3K9_Somatic
      ),
      "Supplemental_Data/SD03_Bulk_Sequence_Stats.xlsx"
    )
  ),

  # Profiles of ChIC faceted by quantification
  tar_map(
    experiment.driver %>%
      cross_join(tibble(plot_name=c("TSS", "Paneled"))) %>%
      mutate(
        named_tss_data = mapply(
          \(celltype, plot_name) call(
            "setNames",
            rlang::syms(str_glue("chic.heatmap.{tolower(plot_name)}_{chic.mark.data$mark}_{celltype}_CN_chr")),
            chic.mark.data$mark
          ),
          celltype,
          plot_name,
          SIMPLIFY=F
        ),
        quartile.factor = rlang::syms(str_glue("quartile.factor_{celltype}"))
      ),
    names = celltype | plot_name,
    tar_target(chic.experiment.tss.heatmaps, named_tss_data),
    tar_target(
      facet_genes,
      with(
        read.csv(assay.data.sc),
        tibble(
          facet = chr %>%
            fct_recode(`2`="2L", `2`="2R", `3`="3L", `3`="3R") %>%
            factor(c("X", "2", "3", "4")),
          quant = quartile.factor %>%
            fct_recode(off="Q1", low="Q2", med="Q3", high="Q4"),
          activity = quartile.factor %>%
            fct_recode(off="Q1", active="Q2", active="Q3", active="Q4"),
          gene = X
        ) %>%
          subset(!is.na(facet))
      )
    ),
    tar_target(
      sc_chr_quartile_data,
      mapply(
        chic_heatmap_facet_genes,
        named_tss_data,
        list(subset(facet_genes, select=c(facet, quant, gene))),
        SIMPLIFY=F
      ) %>%
        bind_rows(.id = "mark") %>%
        mutate(mark = factor(mark, chic.mark.data$mark)),
      format = "parquet"
    ),
    tar_target(
      sc_chr_active_data,
      mapply(
        chic_heatmap_facet_genes,
        named_tss_data,
        list(subset(facet_genes, select=c(facet, activity, gene))),
        SIMPLIFY=F
      ) %>%
        bind_rows(.id = "mark") %>%
        mutate(mark = factor(mark, chic.mark.data$mark)),
      format = "parquet"
    ),
    tar_target(
      sc_quartile_data,
      mapply(
        chic_heatmap_facet_genes,
        named_tss_data,
        list(subset(facet_genes, select=c(quant, gene))),
        SIMPLIFY=F
      ) %>%
        bind_rows(.id = "mark") %>%
        mutate(mark = factor(mark, chic.mark.data$mark)),
      format = "parquet"
    ),
    tar_target(
      sc_active_data,
      mapply(
        chic_heatmap_facet_genes,
        named_tss_data,
        list(subset(facet_genes, select=c(activity, gene))),
        SIMPLIFY=F
      ) %>%
        bind_rows(.id = "mark") %>%
        mutate(mark = factor(mark, chic.mark.data$mark)),
      format = "parquet"
    ),
    tar_target(
      sc_extended_data,
      tibble(
        gene_list = names(cpm_gene_lists_extended)[1],
        mapply(
          chic_heatmap_facet_genes,
          named_tss_data,
          list(tibble(gene = cpm_gene_lists_extended[[1]])),
          SIMPLIFY=F
        ) %>%
          bind_rows(.id = "mark") %>%
          mutate(mark = factor(mark, chic.mark.data$mark))
      ),
      pattern = map(cpm_gene_lists_extended),
      format = "parquet"
    )
  ),
  tar_target(
    sc_chr_nucleosome_data_Germline_TSS,
    tibble(
      chic_heatmap_facet_genes(chic.heatmap.tss.nucleosome_H3K27_Germline_CN_chr, subset(facet_genes_Germline_TSS, select=c(facet,activity,gene))),
      mark = ""
    ),
    format = "parquet"
  ),
  tar_target(
    sc_chr_nucleosome_data_Somatic_TSS,
    tibble(
      chic_heatmap_facet_genes(chic.heatmap.tss.nucleosome_H3K27_Somatic_CN_chr, subset(facet_genes_Somatic_TSS, select=c(facet,activity,gene))),
      mark = ""
    ),
    format = "parquet"
  ),

  # Temp targets for loading chic.bw.2 and producing repli graphics
  tar_target(
    chic.track_chr,
    tibble(
      filename = "chic/nos_H3K4.new.FE.bw",
      track = list(import(BigWigFile(filename)))
    ) %>%
      with(
        GRanges(
          track[[1]]@seqnames[!grepl("GAL80|FBte", track[[1]]@seqnames)] %>% droplevels,
          track[[1]]@ranges[!grepl("GAL80|FBte", track[[1]]@seqnames)],
          seqlengths=subset(seqlengths(track[[1]]), !grepl("GAL80|FBte", names(seqlengths(track[[1]]))))
        )
      )
  ),
  tar_target(
    chic.track_masked,
    tibble(
      filename = "chic/nos_H3K4.new.FE.bw",
      track = list(import(BigWigFile(filename)))
    ) %>%
      with(
        GRanges(
          track[[1]]@seqnames[!grepl("GAL80", track[[1]]@seqnames)] %>% droplevels,
          track[[1]]@ranges[!grepl("GAL80", track[[1]]@seqnames)],
          seqlengths=subset(seqlengths(track[[1]]), !grepl("GAL80", names(seqlengths(track[[1]]))))
        )
      )
  ),
  tar_target(chic.results_nos_chr, list(H3K4="chic/nos_H3K4.new.FE.bw", H3K27="chic/nos_H3K27.new.FE.bw", H3K9="chic/nos_H3K9.new.FE.bw")),
  tar_target(chic.results_tj_chr, list(H3K4="chic/tj_H3K4.new.FE.bw", H3K27="chic/tj_H3K27.new.FE.bw", H3K9="chic/tj_H3K9.new.FE.bw")),
  tar_target(chic.results_nos_masked, chic.results_nos_chr),
  tar_target(chic.results_tj_masked, chic.results_tj_chr),

  # Nucleosome positioning analysis.
  tar_file(
    nucleosome_normal_bin,
    tibble(
      make_cmd = processx::run("make", wd="NOrMAL")$status,
      make_test_cmd = processx::run("make", "test", wd="NOrMAL")$status,
      bin = "NOrMAL/NOrMAL"
    )$bin
  ),
  tar_file(
    nucleosome_normal_config, "NOrMAL/config.txt"
  ),
  tar_map(
    tibble(
      cross_join(experiment.driver, tibble(chr = names(chr.lengths))),
      bulk_reads = mapply(
        \(driver, chr) str_glue(
          "bulk_reads_",
          chr,
          "_chic.bam_{chic.samples$sample[chic.samples$driver == driver & chic.samples$molecule == 'H3' & chic.samples$group == 'H3K27']}_chr"
        ) %>%
          rlang::syms(),
        driver,
        chr,
        SIMPLIFY=FALSE
      )
    ),
    names = chr | celltype,
    tar_file(
      nucleosome_analysis,
      tibble(
        do.call(rbind, bulk_reads) %>%
          paired_end_reads_to_fragment_lengths %>%
          subset(between(length, 100, 200)),
        forw_file = tempfile("forw"),
        do_write_forw_file = write.table(data.frame(format(pos, sci=F)), forw_file[1], quote=F, row.names=F, col.names=F),
        rev_file = tempfile("rev"),
        do_write_rev_file = write.table(data.frame(format(fragment_end_crick, sci=F)), rev_file[1], quote=F, row.names=F, col.names=F),
        output_file = str_glue("chic/NOrMAL/", celltype, "/", "Analysis_", chr, ".txt"),
        mk_output = dir.create(dirname(output_file[1]), showW=F, rec=T),
        do_run_normal = processx::run(
          nucleosome_normal_bin,
          c(
            nucleosome_normal_config,
            forw_file[1],
            rev_file[1],
            output_file[1]
          ),
          stdout="",
          stderr=""
        )$status
      )$output_file[1]
    )
  ),
  tar_map(
    tibble(
      experiment.driver,
      nucleosome_analysis = sapply(
        celltype,
        \(celltype) rlang::syms(
          str_glue("nucleosome_analysis_{names(chr.lengths)}_{celltype}")
        ),
        simplify=F
      )
    ),
    names = celltype,
    tar_target(
      nucleosome_analysis_bed,
      nucleosome_normal_to_bed(
        setNames(nucleosome_analysis, names(chr.lengths)),
        str_glue("chic/NOrMAL/", celltype, "_Pos.bed")
      )
    ),
    tar_target(
      nucleosome_analysis_bw,
      nucleosome_analysis_bed %>% bed_to_mark_bigwig(str_glue("chic/NOrMAL/", celltype, "_Pos.bw"))
    )
  ),
  tar_file(
    plot.chic.nucleosome.median,
    save_figures(
      "figure/Both-Cell-Types",
      ".pdf",
      tribble(
        ~rowname, ~figure, ~width, ~height,
        "Median-Nucleosome-Phasing",
        granges_plot_nuc_median_phasing(
          c(
            Germline=nucleosome_analysis_bed_Germline,
            Somatic=nucleosome_analysis_bed_Somatic
          ),
          legend.position = "none"
        ),
        6,
        3
      )
    )
  )
)
