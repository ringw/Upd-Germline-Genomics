library(dplyr)

cell_type_violin_colors <- c(Germline="#96C586", Somatic="#A97AAC")

chic.samples = read.csv('chic/chic_samples.csv') %>%
  subset(sample != "" & !sapply(rejected, isTRUE)) %>%
  mutate(
    sample_glob = sample,
    sample = sample %>%
      replace(
        str_length(ident) != 0,
        ident[str_length(ident) != 0]
      )
  )
chic.samples.dimreduc <- chic.samples %>%
  subset(!grepl("ChIP", group) & grepl("H3", molecule))

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
})

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

targets.chic <- list(
  # Align ChIC BAM script. We do not filter by mapq yet (light filtering), we
  # will do all of that after reading the BAM file into an R data frame.
  tar_map(
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
              list(
                output_path = str_glue("chic/{name}/{group}/{sample}.bam"),
                batch=batch,
                sample=sample
              )
            ),
          {
            run(
              "bash",
              c(
                "-i",
                align_chic_lightfiltering,
                str_replace(bowtie[1], "\\..*", ""),
                paste0(batch, "/", sample_glob, "_R1_001.fastq.gz"),
                paste0(batch, "/", sample_glob, "_R2_001.fastq.gz"),
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
  ),

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
        levels(
          seqnames(chic.tile_chr_granges_list[[1]]) %>%
            as.factor() %>%
            # We are going to split chic.tile.diameter_1000_masked and then
            # expect the first elements in the list to be masked.lengths!
            fct_relevel("2L_Histone_Repeat_Unit", after=7L)
        ),
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
      chic.granges.dinucleosome.diameter_40,
      tibble(
        reads = append(bulk_reads_split, list(bulk_reads_misc)) %>%
          sapply(
            \(df) df %>%
              pos_fixer_callable %>%
              filter(between(length, 200, 500), mapq >= 20),
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

  # H3-GFP ChIC experiment (extract some statistics - fragment size - from some
  # H3 samples which are representative of all of the H3 samples overall).
  tar_map(
    tribble(
      ~celltype, ~bulk_reads,
      "Germline",
      rlang::syms(
        with(
          expand.grid(
            n = c(names(chr.lengths), "misc"),
            sample = c("GC3768007_S7_L001", "GC3768009_S9_L001", "GC76045515_S5_L001")
          ),
          str_glue("bulk_reads_{n}_chic.bam_{sample}_chr")
        )
      ),
      "Somatic",
      rlang::syms(
        with(
          expand.grid(
            n = c(names(chr.lengths), "misc"),
            sample = c("GC3768010_S10_L001", "GC3768013_S13_L001", "GC3768014_S14_L001")
          ),
          str_glue("bulk_reads_{n}_chic.bam_{sample}_chr")
        )
      ),
    ),
    names = celltype,
    tar_file(
      fig.fragment.size,
      save_figures(
        str_glue("figure/", celltype),
        ".pdf",
        tribble(
          ~rowname, ~figure, ~width, ~height,
          "CHIC-H3-Fragment-Size",
          histogram_paired_end_fragment_size(bulk_reads),
          4,
          2.5,
          "CHIC-H3-Fragment-Size-Chr",
          histogram_paired_end_fragment_size(bulk_reads, faceted=TRUE),
          4,
          11,
        )
      )
    )
  ),

  # Chromatin Mark Experiments (H3 and mark ChIP paired samples -> estimate
  # effect i.e. chromatin mark L2FC and apply post-hoc processing i.e. quantify
  # H3 and mark values at each window and then smooth).
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
      chic.experiment.granges.model.frame,
      tibble(
        rowname = paste(
          mark,
          name,
          ifelse(grepl(mark, molecule), "M", molecule),
          model_matrix_rep,
          sep="_"
        ),
        molecule = factor(molecule),
        rep = factor(model_matrix_rep)
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
          rep(each = length(chic.experiment.granges)) +
        log(ranges(chic.experiment.granges)@width) +
        # Update to the GLM offset for male sex. For female, we would probably
        # not introduce an offset (2 X chromosomes). Negative offset will make
        # us pull up the log-quantification by this amount.
        ifelse(
          grepl("[XY]", as.character(seqnames(chic.experiment.granges))),
          -log(2),
          0
        )
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
          as_tibble(read.csv(assay.data.sc)) %>%
            arrange(desc(Upd_cpm[X, tolower(name)])),
          grep("FSeq_Input", chic.bw.tracks, val=T) %>% BigWigFile %>% import,
          mask_threshold = 0
        )
    ),
    tar_target(
      chic.heatmap.tss.nucleosome,
      if (!grepl("peakcalling", bp_name))
        track_to_heatmap(
          grep("Imputed_Input", chic.bw.tracks, val=T) %>% BigWigFile %>% import,
          as_tibble(read.csv(assay.data.sc)),
          grep("FSeq_Input", chic.bw.tracks, val=T) %>% BigWigFile %>% import,
          mask_threshold = -Inf
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
  # Pericentromere / chromosome arms reference. This can be generated from the
  # 5-state chromatin model by searching for classic heterochromatin (GREEN)
  # near the Watson or Crick end (whichever is flanking the centromere). The
  # actual pericentromere of interest is very sparse in coverage of chromatin
  # states. We will start at the centromere end, including that bp, and seek to
  # find the next chromatin state on the chromosome, which should be GREEN. We
  # can fill a short gap (< 1 kb) of another chromatin state near the
  # pericentromere, but otherwise we stop when we find a non-GREEN window.
  # https://flybase.org/reports/FBsf0000579224.htm
  tar_target(
    chromosome_pericetromere_label,
    GRanges(
      c(
        "2L:22192401-23513712",
        "2R:1-5651400",
        "3L:23154101-28110227",
        "3R:1-4229200"
      )
    )
  ),
  tar_target(
    chromosome_arms_diameter_1000,
    {
      overlaps <- as.list(
        findOverlaps(
          chromosome_pericetromere_label, chic.tile.diameter_1000_chr
        )
      )
      GRanges(
        chic.tile.diameter_1000_chr,
        group = factor(
          seqnames(chic.tile.diameter_1000_chr),
          c("2L", "2LC", "2RC", "2R", "3L", "3LC", "3RC", "3R")
        ) %>%
          replace(
            unlist(overlaps),
            rep(
              paste0(seqnames(chromosome_pericetromere_label), "C"),
              sapply(overlaps, length)
            )
          )
      )
    }
  ),
  tar_target(
    chromosome_arms_diameter_40,
    {
      overlaps <- as.list(
        findOverlaps(
          chromosome_pericetromere_label, chic.tile.diameter_40_chr
        )
      )
      GRanges(
        chic.tile.diameter_40_chr,
        group = factor(
          seqnames(chic.tile.diameter_40_chr),
          c("X", "2L", "2LC", "2RC", "2R", "3L", "3LC", "3RC", "3R", "4", "Y")
        ) %>%
          replace(
            unlist(overlaps),
            rep(
              paste0(seqnames(chromosome_pericetromere_label), "C"),
              sapply(overlaps, length)
            )
          )
      )
    }
  ),
  # Tracks which are once per experiment. No parameter options, those have been
  # fixed to particular tracks now!
  tar_map(
    tibble(
      dplyr::rename(chic.experiments, celltype="name"),
      chic.experiment.quantify_peakcalling.sharp = rlang::syms(
        str_glue("chic.experiment.quantify_{mark}_{celltype}_peakcalling.sharp_chr")
      ),
      chic.experiment.quantify_peakcalling.broad = rlang::syms(
        str_glue("chic.experiment.quantify_{mark}_{celltype}_peakcalling.broad_chr")
      ),
      chic.bw.track.wide = rlang::syms(
        str_glue("chic.bw.track.wide_{mark}_{celltype}_CN_chr")
      ),
      chic.experiment.quantify.smooth_bw2000 = rlang::syms(
        str_glue("chic.experiment.quantify.smooth_bw2000_{mark}_{celltype}_CN_chr")
      )
    ),
    names = mark | celltype,
    tar_file(
      chic.bw.track.peaks,
      reduce_peaks_2_tracks(
        chic.experiment.quantify_peakcalling.broad,
        chic.experiment.quantify_peakcalling.sharp
      ) %>%
        peaks_to_genome_coverage() %>%
        export(BigWigFile(with(list(mark=mark, celltype=celltype), str_glue("chic/chr/{celltype}_{mark}_CN_Enrichment_Peaks.bw")))) %>%
        as.character,
      packages = c(
        "dplyr",
        "GenomicRanges",
        "rtracklayer",
        "S4Vectors"
      )
    ),
    tar_target(
      chic.gene.enrichment.head,
      reduce_peaks_2_tracks(
        chic.experiment.quantify_peakcalling.broad,
        chic.experiment.quantify_peakcalling.sharp
      ) %>%
        gtf_granges_extended_from_tss(as_tibble(read.csv(assay.data.sc)))
    ),
    tar_file(
      fig.track,
      save_figures(
        paste0("figure/", celltype),
        ".pdf",
        tribble(
          ~rowname, ~figure, ~width, ~height,
          paste0("Track-Plot-", mark),
          plot_track(
            import(BigWigFile(chic.bw.track.wide)),
            name = "L2FC",
            limits = c(-1, 1),
            breaks = c(-1, 0, 1)
          ),
          5.75,
          4
        )
      ),
      packages = tar_option_get("packages") %>% c("cowplot", "grid", "gtable")
    ),
    tar_target(
      enriched.chromosomes.data,
      tibble(
        label = chromosome_arms_diameter_1000$group %>%
          factor(
            levels = c("X", levels(.), "4", "Y")
          ) %>%
          replace(
            which(seqnames(chic.tile.diameter_1000_chr) %in% c("4", "X", "Y")),
            seqnames(chic.tile.diameter_1000_chr)[
              seqnames(chic.tile.diameter_1000_chr) %in% c("4", "X", "Y")
            ] %>%
              as.character
          ),
        start = start(chic.tile.diameter_1000_chr),
        end = end(chic.tile.diameter_1000_chr),
        L2FC = import(BigWigFile(chic.bw.track.wide))$score,
        # H3 track: Need to re-apply the downsampling here, similar to the
        # "Wide L2FC" bw track that has already been written to disk.
        input = approx_track(
          GRanges(
            chic.tile.diameter_40_score_chr,
            score=chic.experiment.quantify.smooth_bw2000$score.molH3
          ),
          chic.tile.diameter_1000_chr
        )$score,
      ) %>%
        subset(!is.na(label)),
      format = "parquet"
    ),
    tar_target(
      enriched.chromosomes,
      ggplot(
        enriched.chromosomes.data,
        aes(label, L2FC, fill=substr(label, 1, 1))
      ) +
        geom_violin() +
        geom_boxplot(outlier.shape = NA, fill = "transparent") +
        scale_fill_hue(h.start = 75, c = 50, l = 80) +
        coord_cartesian(NULL, c(-1.2, 1.2)) +
        labs(
          x = "Chromosome"
        ) +
        theme(aspect.ratio = 1/2, legend.position = "none")
    ),
    tar_file(
      fig.enriched.chromosomes,
      save_figures(
        paste0("figure/", celltype),
        ".pdf",
        tribble(
          ~rowname, ~figure, ~width, ~height,
          paste0("Enriched-Chromosome-Regions-", mark),
          enriched.chromosomes,
          5.5,
          3.5
        )
      )
    )
  ),
  # SummarizedExperiment which is once per mark. Grabbing a particular track and
  # a particular bowtie reference (chr).
  tar_map(
    tibble(
      chic.mark.data,
      chic.experiment.granges_Germline = rlang::syms(
        str_glue("chic.experiment.granges_{mark}_Germline_CN_chr")
      ),
      chic.experiment.granges_Somatic = rlang::syms(
        str_glue("chic.experiment.granges_{mark}_Somatic_CN_chr")
      ),
      chic.experiment.granges.offset_Germline = rlang::syms(
        str_glue("chic.experiment.granges.offset_{mark}_Germline_CN_chr")
      ),
      chic.experiment.granges.offset_Somatic = rlang::syms(
        str_glue("chic.experiment.granges.offset_{mark}_Somatic_CN_chr")
      ),
      chic.experiment.granges.model.frame_Germline = rlang::syms(
        str_glue("chic.experiment.granges.model.frame_{mark}_Germline_CN_chr")
      ),
      chic.experiment.granges.model.frame_Somatic = rlang::syms(
        str_glue("chic.experiment.granges.model.frame_{mark}_Somatic_CN_chr")
      ),
    ),
    names = mark,
    # Grab counts grouped by chromosome_arms; vars to regress; and size factors.
    tar_target(
      enriched.chromosomes.compare,
      enriched_chromosomes_build_summarized_experiment(
        cbind(
          elementMetadata(chic.experiment.granges_Germline),
          elementMetadata(chic.experiment.granges_Somatic)
        ),
        chromosome_arms_diameter_40,
        bind_rows(
          list(
            GSC=chic.experiment.granges.model.frame_Germline,
            CySC=chic.experiment.granges.model.frame_Somatic
          ),
          .id="celltype"
        ),
        chic.experiment.granges.offset_Germline,
        chic.experiment.granges.offset_Somatic
      )
    ),
    # Model matrix with our contrasts for the batch effect. We have set up one
    # contrast for all of the GSC & CySC batches to be precisely the GSC-CySC
    # cell type difference. So when creating a model matrix which is singular,
    # we remove the column which we created manually which is redundant.
    tar_target(
      enriched.chromosomes.compare.model.matrix,
      model.matrix(
        ~ molecule * celltype + rep,
        colData(enriched.chromosomes.compare)
      ) %>%
        subset(select = -repCelltype)
    ),
    # Fit GLM, reduced model, and test, yielding p-value, for celltype - L2FC
    # interaction coefficient.
    tar_target(
      enriched.chromosomes.compare.test,
      enriched.chromosomes.compare %>%
        glm_gp(
          enriched.chromosomes.compare.model.matrix,
          size_factors = .$sf
        ) %>%
        test_de(
          replace(
            rep(0, ncol(.$Beta)),
            ncol(.$Beta),
            1
          )
        )
    )
  ),
  tar_target(
    enriched.chromosomes.compare.signif,
    enriched_chromosomes_compare_signif(
      list(
        H3K4 = enriched.chromosomes.compare.test_H3K4,
        H3K27 = enriched.chromosomes.compare.test_H3K27,
        H3K9 = enriched.chromosomes.compare.test_H3K9
      )
    )
  ),
  # Enriched Chromosomes significance marks, to be added to the plot later.
  tar_target(
    enriched.chromosomes.both.data,
    rbind(
      tibble(celltype = "Germline", mark = "H3K4me3", enriched.chromosomes.data_H3K4_Germline),
      tibble(celltype = "Germline", mark = "H3K27me3", enriched.chromosomes.data_H3K27_Germline),
      tibble(celltype = "Germline", mark = "H3K9me3", enriched.chromosomes.data_H3K9_Germline),
      tibble(celltype = "Somatic", mark = "H3K4me3", enriched.chromosomes.data_H3K4_Somatic),
      tibble(celltype = "Somatic", mark = "H3K27me3", enriched.chromosomes.data_H3K27_Somatic),
      tibble(celltype = "Somatic", mark = "H3K9me3", enriched.chromosomes.data_H3K9_Somatic)
    ) %>%
      group_by(label, start) %>%
      filter(all(input > 0.1)),
    format = "parquet"
  ),
  tar_target(
    enriched.chromosomes.both,
    sapply(
      c("X", "2", "3", "4", "Y"),
      purrr::partial(plot_enriched_chromosomes_combined, enriched.chromosomes.both.data),
      simplify=FALSE
    )
  ),
  tar_file(
    fig.enriched.chromosomes.both,
    save_figures(
      "figure/Both-Cell-Types",
      ".pdf",
      tribble(
        ~rowname, ~figure, ~width, ~height,
        "Enriched-Chromosome-Regions-X",
        enriched.chromosomes.both$X %>% annotate_enriched_chromosomes_combined(),
        3, 3,
        "Enriched-Chromosome-Regions-2",
        enriched.chromosomes.both$`2` %>% annotate_enriched_chromosomes_combined(),
        8, 3,
        "Enriched-Chromosome-Regions-3",
        enriched.chromosomes.both$`3` %>% annotate_enriched_chromosomes_combined(),
        8, 3,
        "Enriched-Chromosome-Regions-4",
        enriched.chromosomes.both$`4` %>% annotate_enriched_chromosomes_combined(),
        3, 3,
        "Enriched-Chromosome-Regions-Y",
        enriched.chromosomes.both$Y %>% annotate_enriched_chromosomes_combined(),
        3, 3,
      )
    )
  ),

  tar_target(
    chic.experiment.nucleosomes,
    chic.tile.diameter_40_score_chr %>%
      `elementMetadata<-`(
        value = cbind(
          elementMetadata(chic.experiment.granges_H3K4_Germline_CN_chr) %>%
            subset(select=c(score.H3_Rep1, score.H3_Rep2)),
          elementMetadata(chic.experiment.granges_H3K4_Somatic_CN_chr) %>%
            subset(select=c(score.H3_Rep1, score.H3_Rep3))
        ) %>%
          `colnames<-`(
            value = as.character(interaction(rep(c("GSC", "CySC"), each=2), 1:2))
          )
      )
  ),
  tar_target(
    chic.experiment.dinucleosomes,
    chic.tile.diameter_40_score_chr %>%
      `elementMetadata<-`(
        # These sample ids must match the above target & rep number found in chic_samples.csv!
        value = data.frame(
          GSC.1 = chic.granges.dinucleosome.diameter_40_GC2894001_S1_L001_CN_chr$score,
          GSC.2 = chic.granges.dinucleosome.diameter_40_GC2894002_S2_L001_CN_chr$score,
          CySC.1 = chic.granges.dinucleosome.diameter_40_GC3772016_S1_L002_CN_chr$score,
          CySC.2 = chic.granges.dinucleosome.diameter_40_GC3772018_S3_L002_CN_chr$score
        )
      )
  ),
  tar_target(
    chic.experiment.nucleosomes.offset,
    matrix(
      -log(2) * as.numeric(grepl("[XY]", seqnames(chic.tile.diameter_40_chr))),
      nrow = length(chic.tile.diameter_40_chr),
      ncol = 4,
      dimnames = list(NULL, colnames(elementMetadata(chic.experiment.nucleosomes)))
    ) +
      list(
        GSC.1 = chic.granges.diameter_40_GC2894001_S1_L001_CN_chr$score,
        GSC.2 = chic.granges.diameter_40_GC2894002_S2_L001_CN_chr$score,
        CySC.1 = chic.granges.diameter_40_GC3772016_S1_L002_CN_chr$score,
        CySC.2 = chic.granges.diameter_40_GC3772018_S3_L002_CN_chr$score
      ) %>%
        sapply(
          \(tr) interp.median(
            tr[as.logical(seqnames(chic.tile.diameter_40_chr) %in% c("2L", "2R", "3L", "3R", "4"))],
            na.rm = T
          )
        ) %>%
        log %>%
        rep(each = length(chic.tile.diameter_40_chr)),
    packages = tar_option_get("packages") %>% c("psych")
  ),
  tar_target(
    chic.experiment.quantify.nucleosomes,
    chic_quantify(
      chic.experiment.nucleosomes,
      structure(
        c(2L, 2L, 1L, 1L),
        levels=c("H3.CySC", "H3.GSC"),
        class = "factor"
      ),
      rep(1L, 4),
      chic.experiment.nucleosomes.offset,
      test_de = TRUE,
      test_wald = TRUE
    )
  ),
  tar_target(
    chic.experiment.quantify.dinucleosomes,
    chic_quantify(
      chic.experiment.dinucleosomes,
      structure(
        c(2L, 2L, 1L, 1L),
        levels=c("H3.CySC", "H3.GSC"),
        class = "factor"
      ),
      rep(1L, 4),
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
  tar_map(
    tribble(
      ~celltype, ~chic.test.nucleosomes.fix, ~score_column,
      "Germline", rlang::sym("chic.test.nucleosomes.fix_Germline"), "score.molH3.GSC",
      "Somatic", rlang::sym("chic.test.nucleosomes.fix_Somatic"), "score.molH3.CySC",
    ),
    names = celltype,
    tar_target(
      chic.fix.quantify.monosomes,
      granges_approx(
        chic.test.nucleosomes.fix$Nucleosomes,
        GRanges(
          chic.tile.diameter_40_chr,
          score = elementMetadata(chic.experiment.quantify.nucleosomes)[, score_column]
        )
      )
    ),
    tar_target(
      chic.fix.quantify.dinucleosomes,
      granges_approx(
        chic.test.nucleosomes.fix$Nucleosomes,
        GRanges(
          chic.tile.diameter_40_chr,
          score = elementMetadata(chic.experiment.quantify.dinucleosomes)[, score_column]
        )
      )
    ),
    tar_target(
      plot.chic.nuc.boxplot.monosomes,
      chic.fix.quantify.monosomes$score %>%
        `/`(
          chic.fix.quantify.monosomes$score[
            as.logical(seqnames(chic.fix.quantify.monosomes) %in% c("2L", "2R", "3L", "3R", "4"))
          ] %>%
            median
        ) %>%
        GRanges(chic.test.nucleosomes.fix$Nucleosomes, score=.) %>%
        chic_stats_boxplot() +
        theme(aspect.ratio = 1/2) +
        labs(title = "Monosomes")
    ),
    tar_target(
      plot.chic.nuc.boxplot.dinucleosomes,
      chic.fix.quantify.dinucleosomes$score %>%
        `/`(
          chic.fix.quantify.monosomes$score[
            as.logical(seqnames(chic.fix.quantify.monosomes) %in% c("2L", "2R", "3L", "3R", "4"))
          ] %>%
            median
        ) %>%
        GRanges(chic.test.nucleosomes.fix$Nucleosomes, score=.) %>%
        chic_stats_boxplot() +
        theme(aspect.ratio = 1/2) +
        labs(title = "Dinucleosomes")
    ),
    tar_file(
      fig.chic.nuc.boxplot,
      save_figures(
        str_glue("figure/", celltype),
        ".pdf",
        tribble(
          ~rowname, ~figure, ~width, ~height,
          "CHIC-H3-Box-Plot",
          plot_grid(
            plot.chic.nuc.boxplot.monosomes,
            plot.chic.nuc.boxplot.dinucleosomes,
            nrow = 2
          ),
          6,
          6.75
        )
      ),
      packages = tar_option_get("packages") %>% c("cowplot")
    )
  ),
  tar_file(
    fig.chic.nuc.boxplot.both,
    save_figures(
      "figure/Both-Cell-Types",
      ".pdf",
      tribble(
        ~rowname, ~figure, ~width, ~height,
        "CHIC-H3-Box-Plot",
        plot_grid(
          list(
            Germline=plot.chic.nuc.boxplot.monosomes_Germline$data,
            Somatic=plot.chic.nuc.boxplot.monosomes_Somatic$data
          ) %>%
            bind_rows(.id = "celltype") %>%
            subset(x != "rDNA") %>%
            ggplot(aes(x, y, fill=celltype)) +
            geom_boxplot(fatten=3.5, outlier.shape = NA) +
            scale_y_continuous(
              trans='log',
              breaks=c(0.5, 1, 2, 4)
            ) +
            coord_cartesian(c(0.36, 5.64), c(0.25, 5), expand=F) +
            scale_fill_manual(values=setNames(unlist(chic_line_track_colors), NULL)) +
            theme_bw() +
            theme(aspect.ratio=1/2) +
            labs(title = "Monosomes", x = NULL, y = "Enrichment (vs Auto Monosome Median)"),
          list(
            Germline=plot.chic.nuc.boxplot.dinucleosomes_Germline$data,
            Somatic=plot.chic.nuc.boxplot.dinucleosomes_Somatic$data
          ) %>%
            bind_rows(.id = "celltype") %>%
            subset(x != "rDNA") %>%
            ggplot(aes(x, y, fill=celltype)) +
            geom_boxplot(fatten=3.5, outlier.shape = NA) +
            scale_y_continuous(
              trans='log',
              breaks=c(0.25, 0.5, 1, 2),
              minor_breaks=0.125
            ) +
            coord_cartesian(c(0.36, 5.64), c(0.1, 3), expand=F) +
            scale_fill_manual(values=setNames(unlist(chic_line_track_colors), NULL)) +
            theme_bw() +
            theme(aspect.ratio=1/2) +
            labs(title = "Dinucleosomes", x = NULL, y = "Enrichment (vs Auto Monosome Median)"),
          nrow = 2
        ),
        8,
        6.75
      )
    )
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
        write_fix = GRanges(gr, seqinfo=seqinfo(chic.tile.diameter_40_score_chr)) %>%
          peaks_to_genome_coverage() %>%
          export(BigWigFile(str_glue("chic/chr/{name}_Fix.bw"))) %>%
          as.character
      ) %>%
      unlist
  ),
  tar_target(
    chic.nucleosome.est.query,
    bulk_reads_idxstats_chic.bam_GC2894001_S1_L001_chr %>%
      subset(rname %in% c(names(chr.lengths), "rDNA")) %$%
      GRanges(rname, IRanges(1, width=rlength)) %>%
      slidingWindows(100000L, 100000L) %>%
      sapply(
        # Now extend each window to a good size for estimation.
        \(gr) GRanges(
          seqnames(gr),
          gr %>%
            ranges %>%
            resize(200000L, fix="center") %>%
            restrict(
              1L,
              if (length(gr))
                seqlengths(gr)[as.character(seqnames(gr)[1])]
              else
                1L
            ),
          seqinfo = seqinfo(gr)
        )
      ) %>%
      GRangesList() %>%
      unlist()
  ),
  # Map over cell types for Nucleosome Spacing Estimation.
  tar_map(
    tribble(
      ~celltype, ~p_value_column,
      "Germline", "p_GSC_nuc",
      "Somatic", "p_CySC_nuc"
    ),
    names = celltype,
    # Reduced GRanges of the rough locations of every "nucleosome" (even ones
    # that are actually improbably low-confidence).
    tar_target(
      chic.nucleosome.unprocessed.loc,
      chic.test.nucleosomes[
        which(
          pull(as.data.frame(elementMetadata(chic.test.nucleosomes)), p_value_column) < 0.05
        )
      ] %>%
        split(seqnames(.)) %>%
        sapply(\(gr) if (length(gr)) GRanges(seqnames(gr)[1], IRanges::reduce(ranges(gr)), seqinfo=seqinfo(gr)) else GRanges(seqinfo=seqinfo(gr))) %>%
        GRangesList %>%
        unlist
    ),
    # Mode Nucleosome Spacing estimator using KDE.
    tar_target(
      chic.nucleosome.est,
      GRanges(
        chic.nucleosome.est.query,
        score = findOverlaps(chic.nucleosome.est.query, chic.nucleosome.unprocessed.loc) %>%
          as("List") %>%
          mapply(
            \(inds, nucloc) nucloc[inds] %>%
              diff %>%
              subset(. < 1000) %>%
              (
                function(v) if (length(v) >= 10)
                  density(v, n=16*1024, adjust=0.25) %$%
                    x[which.max(y)]
                else NA
              ),
            .,
            list(mid(chic.nucleosome.unprocessed.loc))
          )
      )
    )
  ),
  # Figure of nucleosome spacing estimate track.
  tar_file(
    fig.chic.nucleosome.est,
    save_figures(
      "figure/Both-Cell-Types",
      ".pdf",
      tribble(
        ~rowname, ~figure, ~width, ~height,
        "Nucleosome-Repeat-Length-Periodicity",
        plot_track_2score(
          GRanges(
            seqnames(nucleosome.repeat.length_Germline),
            ranges(nucleosome.repeat.length_Germline),
            score_1 = nucleosome.repeat.length_Germline$score %>%
              replace(
                nucleosome.repeat.length_Germline$n == 1 &
                  !as.logical(seqnames(nucleosome.repeat.length_Germline) %in% c("4", "Y")),
                NA
              ),
            score_2 = nucleosome.repeat.length_Somatic$score %>%
              replace(
                nucleosome.repeat.length_Somatic$n == 1 &
                  !as.logical(seqnames(nucleosome.repeat.length_Somatic) %in% c("4", "Y")),
                NA
              )
          ),
          name = "bp",
          limits = c(140, 220),
          breaks = c(140, 180, 220)
        ),
        5.75,
        4
      )
    ),
    packages = tar_option_get("packages") %>% c("cowplot", "grid", "gtable")
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
        H3K9_Somatic_Broad = chic.gene.enrichment.broad$H3K9_Somatic,
        H3K4_Germline_Head = chic.gene.enrichment.head_H3K4_Germline,
        H3K27_Germline_Head = chic.gene.enrichment.head_H3K27_Germline,
        H3K9_Germline_Head = chic.gene.enrichment.head_H3K9_Germline,
        H3K4_Somatic_Head = chic.gene.enrichment.head_H3K4_Somatic,
        H3K27_Somatic_Head = chic.gene.enrichment.head_H3K27_Somatic,
        H3K9_Somatic_Head = chic.gene.enrichment.head_H3K9_Somatic,
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
            str_glue("{chic.mark.data$mark}me3")
          ),
          celltype,
          plot_name,
          SIMPLIFY=F
        ),
        chic.heatmap.tss.nucleosome = mapply(
          \(celltype, plot_name)
            if (plot_name == "TSS") {
              rlang::sym(str_glue("chic.heatmap.tss.nucleosome_H3K27_{celltype}_CN_chr"))
            } else {
              numeric(0)
            },
          celltype,
          plot_name
        ),
        quartile.factor = rlang::syms(str_glue("quartile.factor_{celltype}")),
        repli.factor = rlang::syms(str_glue("repli.gene_{celltype}")),
        repli.gene.active = rlang::syms(str_glue("repli.gene.active_{celltype}")),
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
          repli = repli.factor[X],
          repli.active = repli.gene.active[X],
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
        mutate(mark = factor(mark, str_glue("{chic.mark.data$mark}me3"))),
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
        mutate(mark = factor(mark, str_glue("{chic.mark.data$mark}me3"))),
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
        mutate(mark = factor(mark, str_glue("{chic.mark.data$mark}me3"))),
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
        mutate(mark = factor(mark, str_glue("{chic.mark.data$mark}me3"))),
      format = "parquet"
    ),
    tar_target(
      repli_quartile_data,
      mapply(
        chic_heatmap_facet_genes,
        named_tss_data,
        list(subset(facet_genes, select=c(repli, gene))),
        SIMPLIFY=F
      ) %>%
        bind_rows(.id = "mark") %>%
        mutate(mark = factor(mark, str_glue("{chic.mark.data$mark}me3"))),
      format = "parquet"
    ),
    tar_target(
      repli_quartile_active_data,
      mapply(
        chic_heatmap_facet_genes,
        named_tss_data,
        list(subset(facet_genes, select=c(repli.active, gene))),
        SIMPLIFY=F
      ) %>%
        bind_rows(.id = "mark") %>%
        mutate(mark = factor(mark, str_glue("{chic.mark.data$mark}me3"))),
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
          mutate(mark = factor(mark, str_glue("{chic.mark.data$mark}me3")))
      ),
      pattern = map(cpm_gene_lists_extended),
      format = "parquet"
    ),
    tar_target(
      sc_bivalency_data,
      mapply(
        chic_heatmap_facet_genes,
        named_tss_data,
        tibble(
          activity = Upd_cpm[, tolower(celltype)] %>%
            `>=`(5) %>%
            as.factor %>%
            fct_recode(off="FALSE", active="TRUE"),
          H3K4 = pull(chic.gene.enrichment, str_glue("H3K4_", str_to_sentence(celltype))) %>%
            `<`(0.001) %>%
            as.factor %>%
            fct_recode(`~`="FALSE", `+`="TRUE"),
          H3K27 = pull(chic.gene.enrichment, str_glue("H3K27_", str_to_sentence(celltype))) %>%
            `<`(0.001) %>%
            as.factor %>%
            fct_recode(`~`="FALSE", `+`="TRUE"),
          gene = chic.gene.enrichment$symbol
        ) %>%
          filter(!is.na(H3K4), !is.na(H3K27)) %>%
          list,
        SIMPLIFY=F
      ) %>%
        bind_rows(.id = "mark") %>%
        mutate(mark = factor(mark, str_glue("{chic.mark.data$mark}me3"))),
      format = "parquet"
    ),
    tar_target(
      sc_chr_bivalency_data,
      mapply(
        chic_heatmap_facet_genes,
        named_tss_data,
        tibble(
          facet = read.csv(assay.data.sc)$chr %>%
            fct_recode(`2`="2L", `2`="2R", `3`="3L", `3`="3R") %>%
            factor(c("X", "2", "3", "4")),
          activity = Upd_cpm[, tolower(celltype)] %>%
            `>=`(5) %>%
            as.factor %>%
            fct_recode(off="FALSE", active="TRUE"),
          H3K4 = pull(chic.gene.enrichment, str_glue("H3K4_", str_to_sentence(celltype))) %>%
            `<`(0.001) %>%
            as.factor %>%
            fct_recode(`~`="FALSE", `+`="TRUE"),
          H3K27 = pull(chic.gene.enrichment, str_glue("H3K27_", str_to_sentence(celltype))) %>%
            `<`(0.001) %>%
            as.factor %>%
            fct_recode(`~`="FALSE", `+`="TRUE"),
          gene = chic.gene.enrichment$symbol
        ) %>%
          filter(!is.na(facet), !is.na(H3K4), !is.na(H3K27)) %>%
          list,
        SIMPLIFY=F
      ) %>%
        bind_rows(.id = "mark") %>%
        mutate(mark = factor(mark, str_glue("{chic.mark.data$mark}me3"))),
      format = "parquet"
    ),
    tar_target(
      sc_h3k4_data,
      mapply(
        chic_heatmap_facet_genes,
        named_tss_data,
        tibble(
          activity = Upd_cpm[, tolower(celltype)] %>%
            `>=`(5) %>%
            as.factor %>%
            fct_recode(off="FALSE", active="TRUE"),
          GermlineH3K4 = chic.gene.enrichment$H3K4_Germline %>%
            `<`(0.001) %>%
            as.factor %>%
            fct_recode(`~`="FALSE", `+`="TRUE"),
          SomaticH3K4 = chic.gene.enrichment$H3K4_Somatic %>%
            `<`(0.001) %>%
            as.factor %>%
            fct_recode(`~`="FALSE", `+`="TRUE"),
          gene = chic.gene.enrichment$symbol
        ) %>%
          filter(!is.na(GermlineH3K4), !is.na(SomaticH3K4)) %>%
          list,
        SIMPLIFY=F
      ) %>%
        bind_rows(.id = "mark") %>%
        mutate(mark = factor(mark, str_glue("{chic.mark.data$mark}me3"))),
      format = "parquet"
    ),
    # H3 profile. By gene classification (quant).
    tar_target(
      sc_nucleosome_quartile_data,
      chic_heatmap_facet_genes(
        chic.heatmap.tss.nucleosome,
        subset(facet_genes, select=c(quant, gene))
      ),
      format = "parquet"
    ),
    # H3 profile. Filtered so that the other cell type's gene classification
    # must be the "off" level.
    tar_target(
      sc_nucleosome_exclusive_quartile_data,
      chic_heatmap_facet_genes(
        chic.heatmap.tss.nucleosome,
        subset(
          facet_genes,
          Upd_cpm[
            facet_genes$gene,
            c(Germline="somatic", Somatic="germline")[celltype]
          ] < 5,
          select=c(quant, gene)
        )
      ),
      format = "parquet"
    )
  ),
  tar_target(
    sc_chr_nucleosome_data_Germline_TSS,
    tibble(
      chic.heatmap.tss.nucleosome_H3K27_Germline_CN_chr %>%
        chic_heatmap_facet_genes(subset(facet_genes_Germline_TSS, select=c(facet,activity,gene))),
      mark = ""
    ),
    format = "parquet"
  ),
  tar_target(
    sc_chr_nucleosome_data_Somatic_TSS,
    tibble(
      chic.heatmap.tss.nucleosome_H3K27_Somatic_CN_chr %>%
        chic_heatmap_facet_genes(subset(facet_genes_Somatic_TSS, select=c(facet,activity,gene))),
      mark = ""
    ),
    format = "parquet"
  ),
  tar_target(
    sc_chr_nucleosome_data.ExclusiveGermlineGenes_Germline_TSS,
    tibble(
      chic.heatmap.tss.nucleosome_H3K27_Germline_CN_chr %>%
        chic_heatmap_facet_genes(
          facet_genes_Germline_TSS %>%
            subset(
              gene %in% c(
                cpm_gene_lists_extended$ExclusiveGermlineGenes,
                cpm_gene_lists_extended$OffGenes
              ),
              select=c(facet,activity,gene)
            )
        ),
      mark = ""
    ),
    format = "parquet"
  ),
  tar_target(
    sc_chr_nucleosome_data.ExclusiveSomaticGenes_Somatic_TSS,
    tibble(
      chic.heatmap.tss.nucleosome_H3K27_Somatic_CN_chr %>%
        chic_heatmap_facet_genes(
          facet_genes_Somatic_TSS %>%
            subset(
              gene %in% c(
                cpm_gene_lists_extended$ExclusiveSomaticGenes,
                cpm_gene_lists_extended$OffGenes
              ),
              select=c(facet,activity,gene)
            )
        ),
      mark = ""
    ),
    format = "parquet"
  ),

  # ChIC-seq Profile Graphics.
  tar_file(
    fig.bivalent,
    save_figures(
      "figure/Both-Cell-Types",
      ".pdf",
      tribble(
        ~rowname, ~figure, ~width, ~height,
        "CHIC-TSS-AllMarks-Valency",
        chic_panel_gtable_binary(
          list(Germline=sc_bivalency_data_Germline_TSS, Somatic=sc_bivalency_data_Somatic_TSS) %>%
            bind_rows(.id = "genes") %>%
            rev %>%
            tibble(facet = 1) %>%
            subset(between(as.numeric(pos), match("-500", levels(pos)), match("500", levels(pos)))),
          \(...) chic_plot_average_profiles_facet_grid(...) +
            theme(panel.spacing = unit(1.25, "lines"), plot.margin = margin(5.5, 10, 5.5, 5.5)),
          # Pull tlhe chic.gene.enrichment column. The column name is the name
          # of the lysine residue being tested. The modification being tested is
          # always me3.
          c(CelltypeK4="H3K4", CelltypeK27="H3K27"),
          6, 2.5,
          unlist(chic_line_track_colors) %>% setNames(NULL),
          linewidth = rep(0.66, 2)
        ),
        12.5,
        6.25,
        "CHIC-AllMarks-Valency",
        chic_panel_gtable_binary(
          list(Germline=sc_bivalency_data_Germline_Paneled, Somatic=sc_bivalency_data_Somatic_Paneled) %>%
            bind_rows(.id = "genes") %>%
            rev %>%
            tibble(facet = 1),
          chic_plot_paneled_profiles_facet_grid,
          c(CelltypeK4="H3K4", CelltypeK27="H3K27"),
          8.5, 3,
          unlist(chic_line_track_colors) %>% setNames(NULL),
          linewidth = rep(0.66, 2)
        ),
        18,
        7
      )
    ),
    packages = tar_option_get("packages") %>% c("grid", "gtable")
  ),
  tar_file(
    fig.celltype.k4.tss,
    save_figures(
      "figure/Both-Cell-Types",
      ".pdf",
      tribble(
        ~rowname, ~figure, ~width, ~height,
        "CHIC-TSS-H3K4-Compared",
        chic_panel_gtable_binary(
          list(Germline=sc_h3k4_data_Germline_TSS, Somatic=sc_h3k4_data_Somatic_TSS) %>%
            bind_rows(.id = "genes") %>%
            rev %>%
            tibble(facet = 1) %>%
            subset(between(as.numeric(pos), match("-500", levels(pos)), match("500", levels(pos)))),
          \(...) chic_plot_average_profiles_facet_grid(...) +
            theme(panel.spacing = unit(1.25, "lines"), plot.margin = margin(5.5, 10, 5.5, 5.5)),
          c(GermlineK4="GermlineH3K4", SomaticK4="SomaticH3K4"),
          6, 2.5,
          unlist(chic_line_track_colors) %>% setNames(NULL),
          linewidth = rep(0.66, 2)
        ),
        12.5,
        6.25
      )
    ),
    packages = tar_option_get("packages") %>% c("grid", "gtable")
  ),
  tar_map(
    tribble(
      ~celltype, ~sc_bivalency_data, ~plot_color,
      "Germline", rlang::sym("sc_bivalency_data_Germline_TSS"), chic_line_track_colors$germline,
      "Somatic", rlang::sym("sc_bivalency_data_Somatic_TSS"), chic_line_track_colors$somatic,
    ),
    names = celltype,
    tar_file(
      fig.bivalent,
      save_figures(
        str_glue("figure/", celltype),
        ".pdf",
        tribble(
          ~rowname, ~figure, ~width, ~height,
          "CHIC-TSS-K4-K27-Valency",
          sc_bivalency_data %>%
            filter(
              mark %in% c("H3K4me3", "H3K27me3"),
              activity == "active",
              between(
                as.numeric(pos),
                match("-500", levels(pos)),
                match("500", levels(pos))
              )
            ) %>%
            plot_valency_k4_k27(plot_color),
          5.25,
          4.5,
        )
      )
    )
  ),

  # Nucleosome repeat length analysis.
  tar_map(
    tibble(
      cross_join(
        subset(
          chic.samples,
          group == "H3K27" &
            molecule == "H3",
          select = c(sample, driver)
        ),
        tibble(chr = names(chr.lengths))
      ),
      celltype = c(nos="Germline", tj="Somatic")[driver],
      bulk_reads = rlang::syms(
        str_glue("bulk_reads_{chr}_chic.bam_{sample}_chr")
      ),
      bulk_reads_idxstats = rlang::syms(
        "bulk_reads_idxstats_chic.bam_GC3768007_S7_L001_chr"
      )
    ),
    names = celltype | sample | chr,
    tar_target(
      nucleosome.repeat.length.chr,
      bulk_reads %>%
        subset(mapq >= 20) %>%
        paired_end_reads_to_granges() %>%
        GRanges(
          seqlengths = pull(bulk_reads_idxstats, rlength, rname)
        ) %>%
        sliding_nucleosome_qc_for_histogram_usability() %>%
        sliding_nucleosome_repeat_length_analysis(
          width = 2.5 * 1000 * 1000,
          step = 0.5 * 1000 * 1000
        )
    )
  ),
  tar_target(
    nucleosome.repeat.length.list_Germline,
    list(
      c(
        nucleosome.repeat.length.chr_Germline_GC3768007_S7_L001_2L,
        nucleosome.repeat.length.chr_Germline_GC3768007_S7_L001_2R,
        nucleosome.repeat.length.chr_Germline_GC3768007_S7_L001_3L,
        nucleosome.repeat.length.chr_Germline_GC3768007_S7_L001_3R,
        nucleosome.repeat.length.chr_Germline_GC3768007_S7_L001_4,
        nucleosome.repeat.length.chr_Germline_GC3768007_S7_L001_X,
        nucleosome.repeat.length.chr_Germline_GC3768007_S7_L001_Y
      ),
      c(
        nucleosome.repeat.length.chr_Germline_GC3768009_S9_L001_2L,
        nucleosome.repeat.length.chr_Germline_GC3768009_S9_L001_2R,
        nucleosome.repeat.length.chr_Germline_GC3768009_S9_L001_3L,
        nucleosome.repeat.length.chr_Germline_GC3768009_S9_L001_3R,
        nucleosome.repeat.length.chr_Germline_GC3768009_S9_L001_4,
        nucleosome.repeat.length.chr_Germline_GC3768009_S9_L001_X,
        nucleosome.repeat.length.chr_Germline_GC3768009_S9_L001_Y
      ),
      c(
        nucleosome.repeat.length.chr_Germline_GC76045515_S5_L001_2L,
        nucleosome.repeat.length.chr_Germline_GC76045515_S5_L001_2R,
        nucleosome.repeat.length.chr_Germline_GC76045515_S5_L001_3L,
        nucleosome.repeat.length.chr_Germline_GC76045515_S5_L001_3R,
        nucleosome.repeat.length.chr_Germline_GC76045515_S5_L001_4,
        nucleosome.repeat.length.chr_Germline_GC76045515_S5_L001_X,
        nucleosome.repeat.length.chr_Germline_GC76045515_S5_L001_Y
      )
    )
  ),
  tar_target(
    nucleosome.repeat.length_Germline,
    nucleosome.repeat.length.list_Germline %>%
      nucleosome_repeat_length_calling()
  ),
  tar_target(
    nucleosome.repeat.length.peaks_Germline,
    nucleosome.repeat.length.list_Germline %>%
      nucleosome_repeat_length_peak_analysis() %>%
      GenomicRanges::resize(500000, fix="center") %>%
      GenomicRanges::reduce()
  ),
  tar_target(
    nucleosome.repeat.length.list_Somatic,
    list(
      c(
        nucleosome.repeat.length.chr_Somatic_GC3768010_S10_L001_2L,
        nucleosome.repeat.length.chr_Somatic_GC3768010_S10_L001_2R,
        nucleosome.repeat.length.chr_Somatic_GC3768010_S10_L001_3L,
        nucleosome.repeat.length.chr_Somatic_GC3768010_S10_L001_3R,
        nucleosome.repeat.length.chr_Somatic_GC3768010_S10_L001_4,
        nucleosome.repeat.length.chr_Somatic_GC3768010_S10_L001_X,
        nucleosome.repeat.length.chr_Somatic_GC3768010_S10_L001_Y
      ),
      c(
        nucleosome.repeat.length.chr_Somatic_GC3768013_S13_L001_2L,
        nucleosome.repeat.length.chr_Somatic_GC3768013_S13_L001_2R,
        nucleosome.repeat.length.chr_Somatic_GC3768013_S13_L001_3L,
        nucleosome.repeat.length.chr_Somatic_GC3768013_S13_L001_3R,
        nucleosome.repeat.length.chr_Somatic_GC3768013_S13_L001_4,
        nucleosome.repeat.length.chr_Somatic_GC3768013_S13_L001_X,
        nucleosome.repeat.length.chr_Somatic_GC3768013_S13_L001_Y
      ),
      c(
        nucleosome.repeat.length.chr_Somatic_GC3768014_S14_L001_2L,
        nucleosome.repeat.length.chr_Somatic_GC3768014_S14_L001_2R,
        nucleosome.repeat.length.chr_Somatic_GC3768014_S14_L001_3L,
        nucleosome.repeat.length.chr_Somatic_GC3768014_S14_L001_3R,
        nucleosome.repeat.length.chr_Somatic_GC3768014_S14_L001_4,
        nucleosome.repeat.length.chr_Somatic_GC3768014_S14_L001_X,
        nucleosome.repeat.length.chr_Somatic_GC3768014_S14_L001_Y
      )
    )
  ),
  tar_target(
    nucleosome.repeat.length_Somatic,
    nucleosome.repeat.length.list_Somatic %>%
      nucleosome_repeat_length_calling()
  ),
  tar_target(
    nucleosome.repeat.length.peaks_Somatic,
    nucleosome.repeat.length.list_Somatic %>%
      nucleosome_repeat_length_peak_analysis() %>%
      GenomicRanges::resize(500000, fix="center") %>%
      GenomicRanges::reduce()
  ),

  # Dimension Reduction.
  tar_target_raw(
    "chic.dimension.reduction",
    call(
      "prcomp_irlba",
      call(
        "t",
        call(
          "devianceResiduals",
          call(
            "sapply",
            do.call(
              call,
              append(
                list("list"),
                setNames(
                  rlang::syms(str_glue("chic.granges.peakcalling.diameter_500_{chic.samples.dimreduc$sample}_CN_chr")),
                  chic.samples.dimreduc$sample
                )
              ),
              quote=T
            ),
            quote(\(granges) granges$score)
          )
        )
      ),
      n=10,
      center=FALSE
    )
  ),
  tar_file(
    fig.chic.dimension.reduction,
    save_figures(
      "figure/Both-Cell-Types",
      ".pdf",
      tribble(
        ~rowname, ~figure, ~width, ~height,
        "CHIC-PCA",
        ggarrange(
          ggplot(
            tibble(
              as.data.frame(chic.dimension.reduction$x),
              mutate(
                chic.samples.dimreduc,
                group = relevel(factor(group), "H3K4"),
                molecule = fct_relevel(factor(molecule), c("H3", "H3K4me3"))
              )
            ),
            aes(-PC1, PC2, size=driver, color=molecule, shape=group)
          ) +
            geom_point() +
            scale_size_manual(values=c(1.25, 1.75), guide=guide_legend(override.aes=list(size=c(1.25, 2), color=muted(unlist(chic_line_track_colors), c=80, l=65)))) +
            scale_shape_manual(values=c(18, 19, 17), guide=guide_legend(override.aes=list(size=3))) +
            scale_color_manual(values=c("#17d2d5", "#527a36", "#cfdc1d", "#d53840"), guide=guide_legend(override.aes=list(size=3, shape=15))) +
            theme_bw() +
            theme(aspect.ratio = 1) +
            labs(x = "PC1"),
          ggplot(
            tibble(
              as.data.frame(chic.dimension.reduction$x),
              mutate(
                chic.samples.dimreduc,
                group = relevel(factor(group), "H3K4"),
                molecule = fct_relevel(factor(molecule), c("H3", "H3K4me3"))
              )
            ),
            aes(-PC1, PC3, size=driver, color=molecule, shape=group)
          ) +
            geom_point() +
            scale_size_manual(values=c(1.25, 1.75)) +
            scale_shape_manual(values=c(18, 19, 17), guide=guide_legend(override.aes=list(size=3))) +
            scale_color_manual(values=c("#17d2d5", "#527a36", "#cfdc1d", "#d53840"), guide=guide_legend(override.aes=list(size=3))) +
            theme_bw() +
            theme(aspect.ratio = 1) +
            labs(x = "PC1"),
          common.legend = TRUE,
          legend = "right"
        ),
        8,
        3.75
      )
    )
  )
)
