library(dplyr)

cell_type_violin_colors <- c(Germline="#96C586", Somatic="#A97AAC")

# ChIC Sample Sheet after Demultiplexing of Reads ----
# "ident": An R symbol-like name to be used in target names for the sample. It
# may be empty and defaults to the same as "sample".
# "sample": Either R symbol-like filename/accession name, or a glob of such.
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
# Which chic.samples to analyze in PCA. We will elide our alternative "ChIP"
# samples, including the IgG or whole chromatin samples. The molecules that we
# are interested in are: H3, H3K4me3, H3K27me3, H3K9me3.
chic.samples.dimreduc <- chic.samples %>%
  subset(!grepl("ChIP", group) & grepl("H3", molecule))

# The cell types.
chic.fpkm.data <- tibble(name=c("Germline", "Somatic"), driver=c("nos", "tj"))

# The ChIC experiment conditions (antibody to H3 with mark).
chic.mark.data = tribble(~mark, "H3K4", "H3K27", "H3K9")

# ChIC experiments: For every driver x mark.
chic.experiments <- chic.fpkm.data %>%
  cross_join(chic.mark.data) %>%
  # Pull the "chic" (track) target by name for the driver x mark targets.
  mutate(experiment_name = paste(mark, name, sep="_"))

# Alignment Job, producing all of the mapped fragments from bowtie2 ----
targets.chic.align <- tar_map(
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
)

# Sliding Windows based on the chromosomes or masked chr lengths ----
# First, we create a GRanges list covering the ref seqs from 1 to the length,
# per the "samtools idxstats" command result from any one of the bam files.
# Second, we use a variety of slidingWindows which will be applied later.
targets.chic.sliding <- tar_map(
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
)

# Coverage proceeds using GRanges for each sample independently ----
targets.chic.coverage <- tar_map(
  # For each chic sample and reference:
  # Generate the names of targets split by 7 chrs and "misc".
  # Generate the name of the "sharp" GRanges (40bp diameter, 20bp step size)
  # and "broad" GRanges (500bp diameter, 100bp step size).
  # Parameters for ChIC: Fragment length including both 5' end base pairs is to
  # be in [100, 200] bp (this is our monosome length, generous to the DNA being
  # digested a bit too much).
  # For ChIC, we are only to look at fragments filtered by MAPQ >= 20.
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
)

# Regression, Gaussian kernel, and L2FC computation ----
targets.chic.tracks <- list(
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
        test_de = test_de,
        test_wald = test_de
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
        granges_to_normalize = (
          seqnames(chic.tile.diameter_40_score) %in%
            c("2L", "2R", "3L", "3R", "4")
        ) %>%
          as.logical() %>%
          which()
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
            `/`(
              enframe(.) %>%
                dplyr::slice(
                  (
                    seqnames(chic.tile.diameter_40_score) %in%
                      c("2L", "2R", "3L", "3R", "4")
                  ) %>%
                    as.logical() %>%
                    which()
                ) %>%
                deframe() %>%
                median(na.rm=T)
            ) %>%
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
  # For each experiment - plot all of the 7 chromosome arm / spanning the
  # chromosome into both telomeres dm6 sequences.
  tar_map(
    dplyr::rename(chic.experiments, celltype="name") %>%
      mutate(
        chic.bw.track.wide = rlang::syms(
          str_glue("chic.bw.track.wide_{experiment_name}_CN_chr")
        )
      ),
    names = mark | celltype,
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
            breaks = c(-1, 0, 1),
            color = classification_colors_fig4[[tolower(celltype)]]
          ),
          5.75,
          4
        )
      ),
      packages = tar_option_get("packages") %>% c("cowplot", "grid", "gtable")
    )
  ),
  # Collate the report: Each experiment and the distribution on the dm6
  # sequences cut by our definition of the pericentromere.
  tar_file(
    fig.report.enriched.chromosomes,
    save_figures(
      "figure/Both-Cell-Types",
      ".pdf",
      tribble(
        ~rowname, ~figure, ~width, ~height,
        "Enriched-Chromosome-Regions",
        report_chic_l2fe(
          tribble(
            ~celltype, ~mark, ~track,
            "Germline", "H3K4me3", chic.experiment.quantify_H3K4_Germline_peakcalling.broad_chr,
            "Germline", "H3K27me3", chic.experiment.quantify_H3K27_Germline_peakcalling.broad_chr,
            "Germline", "H3K9me3", chic.experiment.quantify_H3K9_Germline_peakcalling.broad_chr,
            "Somatic", "H3K4me3", chic.experiment.quantify_H3K4_Somatic_peakcalling.broad_chr,
            "Somatic", "H3K27me3", chic.experiment.quantify_H3K27_Somatic_peakcalling.broad_chr,
            "Somatic", "H3K9me3", chic.experiment.quantify_H3K9_Somatic_peakcalling.broad_chr,
          ),
          chromosome_pericetromere_label
        ),
        8.5, 11,
      )
    ),
    packages = tar_option_get("packages") %>% c("egg", "grid", "gtable")
  )
)

# Histone Code Summary Targeting Non-Coding Pericentromere ----
# Most interested in the 5-state chromatin model to identify heterochromatin,
# then extend the large heterochromatin feature all the way in the centromere
# direction (where the limitation was just sparse coverage; the newer 9-state
# chromatin model may be even finer/less covered).
# For a statistical test, quickly further condense our 40-bp sliding window 
# samples (we know that they have 50% overlap in windows; double-counting) and
# use those GRanges to get new GRanges for chromosome arms windows.
# Then test DE (glmGamPoi).
targets.chic.pericentromere <- list(
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
  tar_target(chromosome_arms_diameter_1000_chr, chromosome_arms_diameter_1000),
  tar_target(
    chromosome_arms_diameter_1000_masked,
    {
      overlaps <- as.list(
        findOverlaps(
          chromosome_pericetromere_label, chic.tile.diameter_1000_masked
        )
      )
      GRanges(
        chic.tile.diameter_1000_masked,
        group = factor(
          seqnames(chic.tile.diameter_1000_masked),
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
    chromosome_arms_diameter_500_score_chr,
    {
      overlaps <- as.list(
        findOverlaps(
          chromosome_pericetromere_label, chic.tile.diameter_500_score_chr
        )
      )
      GRanges(
        chic.tile.diameter_500_score_chr,
        group = factor(
          seqnames(chic.tile.diameter_500_score_chr),
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
    chromosome_arms_diameter_500_score_masked,
    {
      overlaps <- as.list(
        findOverlaps(
          chromosome_pericetromere_label, chic.tile.diameter_500_score_masked
        )
      )
      GRanges(
        chic.tile.diameter_500_score_masked,
        group = factor(
          seqnames(chic.tile.diameter_500_score_masked),
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
  )
)

# Transposon / Masked FASTA Reference Analysis ----
targets.chic.transposon <- list(
  tar_target(
    enriched.transposable.elements.peakcalling.broad.masked,
    seqnames(chic.experiment.quantify_H3K4_Germline_peakcalling.broad_masked) %>%
      grep("FBte[0-9]{7}", .) %>%
      as.integer()
  ),
  tar_file(
    transposable.elements.metadata,
    "references/transposon_sequence_set_metadata.txt"
  ),
  tar_target(
    enriched.transposable.elements.factor,
    transposable.elements.metadata %>%
      read.table(sep="\t", quote="", skip=1, comment="") %>%
      as_tibble() %>%
      factor_transposable_elements(
        seqnames(chic.experiment.quantify_H3K4_Germline_peakcalling.broad_masked)[
          enriched.transposable.elements.peakcalling.broad.masked
        ] %>%
          as.character(),
        .
      )
  ),
  tar_file(
    fig.enriched.transposable.elements,
    save_figures(
      "figure/Both-Cell-Types",
      ".pdf",
      tribble(
        ~rowname, ~figure, ~width, ~height,
        "Enriched-Transposable-Elements",
        plot_transposable_element_enrich_chromatin(
          chic.experiment.quantify_H3K4_Germline_peakcalling.broad_masked[
            enriched.transposable.elements.peakcalling.broad.masked
          ],
          chic.experiment.quantify_H3K4_Somatic_peakcalling.broad_masked[
            enriched.transposable.elements.peakcalling.broad.masked
          ],
          chic.experiment.quantify_H3K27_Germline_peakcalling.broad_masked[
            enriched.transposable.elements.peakcalling.broad.masked
          ],
          chic.experiment.quantify_H3K27_Somatic_peakcalling.broad_masked[
            enriched.transposable.elements.peakcalling.broad.masked
          ],
          chic.experiment.quantify_H3K9_Germline_peakcalling.broad_masked[
            enriched.transposable.elements.peakcalling.broad.masked
          ],
          chic.experiment.quantify_H3K9_Somatic_peakcalling.broad_masked[
            enriched.transposable.elements.peakcalling.broad.masked
          ]
        ) %>%
          `+`(
            theme(plot.margin = margin(5.5, 130, 5.5, 5.5))
          ) %>%
          set_panel_size(w = unit(2.5, "in"), h = unit(3.75, "in")),
        5, 5,
        "Enriched-Transposable-Elements-DNA",
        plot_transposable_element_enrich_chromatin(
          chic.experiment.quantify_H3K4_Germline_peakcalling.broad_masked[
            enriched.transposable.elements.peakcalling.broad.masked[
              which(enriched.transposable.elements.factor == "DNA transposon")
            ]
          ],
          chic.experiment.quantify_H3K4_Somatic_peakcalling.broad_masked[
            enriched.transposable.elements.peakcalling.broad.masked[
              which(enriched.transposable.elements.factor == "DNA transposon")
            ]
          ],
          chic.experiment.quantify_H3K27_Germline_peakcalling.broad_masked[
            enriched.transposable.elements.peakcalling.broad.masked[
              which(enriched.transposable.elements.factor == "DNA transposon")
            ]
          ],
          chic.experiment.quantify_H3K27_Somatic_peakcalling.broad_masked[
            enriched.transposable.elements.peakcalling.broad.masked[
              which(enriched.transposable.elements.factor == "DNA transposon")
            ]
          ],
          chic.experiment.quantify_H3K9_Germline_peakcalling.broad_masked[
            enriched.transposable.elements.peakcalling.broad.masked[
              which(enriched.transposable.elements.factor == "DNA transposon")
            ]
          ],
          chic.experiment.quantify_H3K9_Somatic_peakcalling.broad_masked[
            enriched.transposable.elements.peakcalling.broad.masked[
              which(enriched.transposable.elements.factor == "DNA transposon")
            ]
          ]
        ) %>%
          `+`(
            theme(plot.margin = margin(5.5, 130, 5.5, 5.5))
          ) %>%
          set_panel_size(w = unit(2.5, "in"), h = unit(3.75, "in")),
        5, 5,
        "Enriched-Transposable-Elements-Retrotransposon",
        plot_transposable_element_enrich_chromatin(
          chic.experiment.quantify_H3K4_Germline_peakcalling.broad_masked[
            enriched.transposable.elements.peakcalling.broad.masked[
              which(enriched.transposable.elements.factor == "retrotransposon")
            ]
          ],
          chic.experiment.quantify_H3K4_Somatic_peakcalling.broad_masked[
            enriched.transposable.elements.peakcalling.broad.masked[
              which(enriched.transposable.elements.factor == "retrotransposon")
            ]
          ],
          chic.experiment.quantify_H3K27_Germline_peakcalling.broad_masked[
            enriched.transposable.elements.peakcalling.broad.masked[
              which(enriched.transposable.elements.factor == "retrotransposon")
            ]
          ],
          chic.experiment.quantify_H3K27_Somatic_peakcalling.broad_masked[
            enriched.transposable.elements.peakcalling.broad.masked[
              which(enriched.transposable.elements.factor == "retrotransposon")
            ]
          ],
          chic.experiment.quantify_H3K9_Germline_peakcalling.broad_masked[
            enriched.transposable.elements.peakcalling.broad.masked[
              which(enriched.transposable.elements.factor == "retrotransposon")
            ]
          ],
          chic.experiment.quantify_H3K9_Somatic_peakcalling.broad_masked[
            enriched.transposable.elements.peakcalling.broad.masked[
              which(enriched.transposable.elements.factor == "retrotransposon")
            ]
          ]
        ) %>%
          `+`(
            theme(plot.margin = margin(5.5, 130, 5.5, 5.5))
          ) %>%
          set_panel_size(w = unit(2.5, "in"), h = unit(3.75, "in")),
        5, 5,
      )
    ),
    packages = tar_option_get("packages") %>% c("egg")
  )
)

# Analysis of an H3 experiment (coverage analysis of monosomes). ----
targets.chic.h3 <- list(
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
  )
)

# Peak calling at TSS using regression of ChIC samples. ----
targets.chic.peakcalling <- list(
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
  # L2FE pulled from the tracks at every gene's TSS.
  tar_target(
    chic.gene.enrichment.l2fc,
    tibble(
      read.csv(assay.data.sc),
      chr.test = chr %>% replace(!(. %in% names(chr.lengths)), "*"),
      windows.broad = {
        ov <- findOverlaps(
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
        )
        sparseVector(i = from(ov), x = to(ov), length = length(quartile.factor_Germline)) %>%
          as.numeric() %>%
          replace(. == 0, NA)
      },
      as_tibble(
        sapply(
          list(
            H3K4_Germline=c(
              chic.experiment.quantify_H3K4_Germline_peakcalling.broad_chr
            ),
            H3K27_Germline=c(
              chic.experiment.quantify_H3K27_Germline_peakcalling.broad_chr
            ),
            H3K9_Germline=c(
              chic.experiment.quantify_H3K9_Germline_peakcalling.broad_chr
            ),
            H3K4_Somatic=c(
              chic.experiment.quantify_H3K4_Somatic_peakcalling.broad_chr
            ),
            H3K27_Somatic=c(
              chic.experiment.quantify_H3K27_Somatic_peakcalling.broad_chr
            ),
            H3K9_Somatic=c(
              chic.experiment.quantify_H3K9_Somatic_peakcalling.broad_chr
            )
          ),
          \(gr) gr$L2FC[windows.broad] %>% replace(!is.finite(.), NA),
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
  )
)

# ChIC-seq querying track by gene, grouping genes, compute mean line plot ----
targets.chic.lineplot <- list(
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
      facet_genes_bivalent,
      reframe(
        chic.gene.enrichment,
        gene = symbol,
        bivalent = interaction(
          cbind(H3K4_Germline, H3K27_Germline) %>%
            `<`(1e-3) %>%
            replace(is.na(.), FALSE) %>%
            rowAlls() %>%
            `&`(quartile.factor_Germline != "Q1"),
          cbind(H3K4_Somatic, H3K27_Somatic) %>%
            `<`(1e-3) %>%
            replace(is.na(.), FALSE) %>%
            rowAlls() %>%
            `&`(quartile.factor_Somatic != "Q1")
        ) %>%
          `levels<-`(
            value = c(
              "off",
              "GSC",
              "CySC",
              "both"
            )
          )
      )
    ),
    tar_target(
      facet_genes_valency,
      reframe(
        chic.gene.enrichment,
        gene = symbol,
        valency = sapply(
          c("H3K4", "H3K27"),
          \(n) get(paste(n, celltype, sep="_")) %>%
            `<`(1e-3) %>%
            replace(is.na(.), FALSE),
          simplify = FALSE
        ) %>%
          do.call(interaction, .) %>%
          replace(
            list(
              Germline = quartile.factor_Germline,
              Somatic = quartile.factor_Somatic
            )[[celltype]] ==
              "Q1",
            "FALSE.FALSE"
          ) %>%
          `levels<-`(
            value = c(
              "~",
              "H3K4",
              "H3K27",
              "bivalent"
            )
          )
      )
    ),
    tar_target(
      facet_genes_diff_timing,
      tibble(
        gene = names(diff.timing.factor),
        timing = diff.timing.factor,
        # activity = diff.timing.isactive,
        activity = quartile.factor %>%
          fct_recode(off="Q1", active="Q2", active="Q3", active="Q4")
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
      repli_diff_timing_data,
      mapply(
        chic_heatmap_facet_genes,
        named_tss_data,
        list(facet_genes_diff_timing),
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
    # H3 profile. By ChIC TSS enrichment K4 & K27.
    tar_target(
      sc_nucleosome_valency_data,
      chic_heatmap_facet_genes(
        chic.heatmap.tss.nucleosome,
        facet_genes_valency
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
            coord_cartesian(NULL, chic_average_profile_limits, ex=F) +
            theme(panel.spacing = unit(1.25, "lines"), plot.margin = margin(5.5, 10, 5.5, 5.5)),
          c(GermlineK4="GermlineH3K4", SomaticK4="SomaticH3K4"),
          6, 2.5,
          unlist(chic_line_track_colors) %>% setNames(NULL),
          linewidth = rep(0.66, 2),
          chic_average_profile_limits = c(0.35, 8.5)
        ),
        12.5,
        6.25
      )
    ),
    packages = tar_option_get("packages") %>% c("grid", "gtable")
  )
)

targets.chic <- list(
  targets.chic.align,
  targets.chic.sliding,
  targets.chic.coverage,
  targets.chic.tracks,
  targets.chic.pericentromere,
  targets.chic.transposon,
  targets.chic.h3,
  targets.chic.peakcalling,
  targets.chic.lineplot
)
