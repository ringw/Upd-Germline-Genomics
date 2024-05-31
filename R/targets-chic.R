library(dplyr)
chic.samples = read.csv('chic/chic_samples.csv') %>%
  subset(sample != "" & !sapply(rejected, isTRUE))

chic.fpkm.data <- tribble(
  ~name, ~contrast, ~driver,
  'Germline', c(1,0,0,0,0,0,0), 'Nos',
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
    ),
    tar_target(
      chic.readcount,
      as.integer(
        c(
          run(
            "bash",
            c(
              "-c",
              paste0(
                "rclone cat sharepoint:'Bio Data'/",
                batch,
                "/",
                sample,
                "_R1_001.fastq.gz",
                " | gunzip -c | wc -l"
              )
            )
          )$stdout,
          run(
            "bash",
            c(
              "-c",
              paste0(
                "rclone cat sharepoint:'Bio Data'/",
                batch,
                "/",
                sample,
                "_R2_001.fastq.gz",
                " | gunzip -c | wc -l"
              )
            )
          )$stdout
        )
      ),
      packages = "processx"
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
      chic.tile_all_chr_diameter_60_granges_list,
      slidingWindows(
        unlist(chic.tile_chr_granges_list), w=60L, s=20L
      )
    ),
    tar_target(
      chic.tile_all_chr_diameter_20_granges_list,
      slidingWindows(
        unlist(chic.tile_chr_granges_list), w=20L, s=20L
      )
    ),
    tar_target(
      chic.tile.diameter_60,
      sapply(
        levels(seqnames(chic.tile_chr_granges_list[[1]])),
        \(n) if (n %in% names(masked.lengths))
          chic.tile_all_chr_diameter_60_granges_list[[n]]
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
      chic.tile.diameter_60_score,
      sapply(
        levels(seqnames(chic.tile_chr_granges_list[[1]])),
        \(n) if (n %in% names(masked.lengths))
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
      chic.tile.diameter_60_score_lookup,
      mapply(
        \(tile_gr, score_gr, lookup_vec) if (length(lookup_vec) > 1)
          c(
            NA,
            lookup_vec,
            if (length(score_gr) >= length(tile_gr) + 2) NA else NULL
          )
        else lookup_vec,
        chic.tile.diameter_60 %>% split(seqnames(.)),
        chic.tile.diameter_60_score %>% split(seqnames(.)),
        seq(length(chic.tile.diameter_60)) %>% split(seqnames(chic.tile.diameter_60)),
        SIMPLIFY=FALSE
      ) %>%
        do.call(c, .)
    )
  ),

  # ChIC loaded BAM files to KDE and the bandwidth-equivalent rectangular windows.
  tar_map(
    cross_join(
      bowtie.refs,
      chic.samples
    ) %>%
      cross_join(
        tribble(
          ~bp_suffix, ~pos_fixer_callable,
          "CN", rlang::sym("paired_end_pos_to_midpoint"),
          "FE", rlang::sym("paired_end_pos_to_5_prime")
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
        chic.tile.diameter_60 = rlang::syms(str_glue("chic.tile.diameter_60_{name}"))
      ),
    names = suffix,

    tar_target(
      chic.granges.diameter_60,
      do.call(
        rbind,
        append(bulk_reads_split, list(bulk_reads_misc))
      ) %>%
        pos_fixer_callable %>%
        filter(
          between(length, 150, 200),
          between(mapq, 20, 254)
        ) %>%
        with(
          count_overlaps_bases(chic.tile.diameter_60, .) %>%
          `metadata<-`(value = list(est_library_size = length(pos)))
        )
    ),

    tar_target(
      chic.sparseVector,
      append(
        sapply(
          bulk_reads_split,
          \(df) bam_cover_paired_end_fragments_bp(
            df, min_mapq = 20, min_fl = 150, max_fl = 200, markdup = TRUE
          )[[1]],
          simplify=FALSE
        ),
        as.list(
          (
            bulk_reads_misc %>%
              paired_end_reads_to_fragment_lengths %>%
              # mutate(rname = rname %>% factor(setdiff(levels(.), names(if (name == "masked") masked.lengths else chr.lengths)))) %>%
              filter(between(mapq, 20, 254), between(length, 150, 200)) %>%
              split(.$rname) %>%
              sapply(\(df) nrow(df))
          )[setdiff(levels(bulk_reads_misc$rname), names(if (name == "masked") masked.lengths else chr.lengths))] %>%
            setNames(
              setdiff(levels(bulk_reads_misc$rname), names(if (name == "masked") masked.lengths else chr.lengths))
            )
        )
      )
    ),
    tar_target(
      chic.chr.fpkm,
      GRanges(
        seqnames = names(chic.sparseVector),
        IRanges(
          start = 1,
          width = pull(bulk_reads_idxstats, "rlength", "rname")[names(chic.sparseVector)]
        ),
        score = sapply(chic.sparseVector, sum) %>%
          `/`(sum(.)) %>%
          `*`(1000^3 / pull(bulk_reads_idxstats, "rlength", "rname")[names(chic.sparseVector)]),
        seqlengths = pull(bulk_reads_idxstats, "rlength", "rname")[names(chic.sparseVector)]
      ) %>%
        `metadata<-`(value = list(est_library_size = sum(sapply(chic.sparseVector, sum)))) %>%
        split(names(chic.sparseVector) %>% factor(., .))
    ),
    tar_target(
      chic.granges.fpkm,
      chic.chr.fpkm %>%
        replace(
          which(sapply(chic.sparseVector, length) > 1),
          chic.sparseVector %>%
            subset(sapply(chic.sparseVector, length) > 1) %>%
            mapply(
              \(rname, v) density_est_granges_approx(
                chic.chr.fpkm[[rname]],
                obs_vector = v,
                bw = 25,
                sample_rate = 20
              ),
              names(.),
              .
            ) %>%
            GRangesList
        )
    )
  ),

  tar_map(
    chic.experiments %>%
      cross_join(dplyr::rename(bowtie.refs, reference="name")) %>%
      cross_join(tibble(bp_suffix = c("CN", "FE"))) %>%
      mutate(driver = driver %>% replace(. == "Nos", "nos")) %>%
      rowwise %>%
      mutate(
        chic.tile.diameter_60_score = rlang::syms(str_glue("chic.tile.diameter_60_score_{reference}")),
        chic.tile.diameter_60_score_lookup = rlang::syms(str_glue("chic.tile.diameter_60_score_lookup_{reference}")),
        chic_sample_names = subset(
          chic.samples$sample,
          (chic.samples$driver %>% replace(. == "Nos", "nos")) == driver &
          # (chic.samples$molecule == "H3" | chic.samples$group == mark)
          chic.samples$group == mark
        ) %>%
          list,
        chic.granges.fpkm = rlang::syms(str_glue("chic.granges.fpkm_{chic_sample_names}_{bp_suffix}_{reference}")) %>%
          list,
        chic.granges.diameter_60 = rlang::syms(str_glue("chic.granges.diameter_60_{chic_sample_names}_{bp_suffix}_{reference}")) %>%
          list,
        experimental_group = (
          chic.samples$group[match(chic_sample_names, chic.samples$sample)]
        ) %>%
          list,
        molecule = (
          chic.samples$molecule[match(chic_sample_names, chic.samples$sample)]
        ) %>%
          list,
        rep = (
          chic.samples$rep[match(chic_sample_names, chic.samples$sample)]
        ) %>%
          list
      ),
    names = mark | name | bp_suffix | reference,
    tar_target(
      chic.experiment.sample.names, chic_sample_names
    ),
    tar_target(
      chic.experiment.granges.fpkm,
      unlist(chic.granges.fpkm[[1]]) %>%
        attributes %>%
        with(
          GRanges(
            seqnames,
            ranges,
            seqlengths = seqlengths(chic.granges.fpkm[[1]]),
            score = sapply(
              chic.granges.fpkm,
              \(lst) unlist(lst)$score
            )
          ) %>%
            `metadata<-`(
              value = list(
                est_library_size = sapply(
                  chic.granges.fpkm,
                  \(lst) lst[[1]]@metadata$est_library_size
                )
              )
            )
        )
    ),
    tar_target(
      chic.experiment.granges.diameter_60,
      chic.granges.diameter_60[[1]] %>%
        attributes %>%
        append(list(molecule=molecule, rep=rep)) %>%
        with(
          GRanges(
            seqnames,
            ranges,
            seqlengths = seqlengths(chic.granges.diameter_60[[1]]),
            score = sapply(chic.granges.diameter_60, \(gr) gr$score) %>%
              matrix(
                nrow = nrow(.),
                dimnames = list(NULL, str_glue("{molecule}_Rep{rep}"))
              )
          ) %>%
            `metadata<-`(
              value = list(
                est_library_size = sapply(
                  chic.granges.diameter_60,
                  \(gr) gr@metadata$est_library_size
                ) %>%
                  setNames(str_glue("{molecule}_Rep{rep}"))
              )
            )
        )
    ),
    tar_target(
      chic.experiment.granges.offset_diameter_60,
      with(
        list(ones = matrix(1, nrow = nrow(chic.experiment.granges.diameter_60@elementMetadata), ncol = ncol(chic.experiment.granges.diameter_50@elementMetadata))),
        GRanges(
          seqnames(chic.experiment.granges.diameter_60),
          ranges(chic.experiment.granges.diameter_60),
          offset = as.matrix(
            Diagonal(x = log(ranges(chic.experiment.granges.diameter_60)@width))
            %*% ones
            + ones %*%
            Diagonal(x = log(chic.experiment.granges.diameter_60@metadata$est_library_size))
          ) - 3 * log(1000)
        )
      )
    ),
    tar_target(
      chic.experiment.quantify,
      if (length(unique(rep)) > 1)
        glm_gp(
          as.matrix(chic.experiment.granges.diameter_60@elementMetadata),
          # We put "molecule" and "rep" in our tar_map, so create new names.
          ~ 0 + mol + R,
          tibble(
            mol = molecule,
            R = rep %>%
              factor %>%
              `contrasts<-`(value = contr.helmert(length(levels(.))))
          ),
          size_factors = 1,
          offset = as.matrix(
            elementMetadata(chic.experiment.granges.offset_diameter_60)
          ),
          overdispersion = "global",
          overdispersion_shrinkage = FALSE,
          verbose = TRUE
        ) %>%
          with(
            GRanges(
              seqnames(chic.experiment.granges.diameter_60),
              ranges(chic.experiment.granges.diameter_60),
              score = as.data.frame(exp(Beta)) %>%
                replace(. < 1e-8, 0),
              seqlengths = seqlengths(chic.experiment.granges.diameter_60)
            ) %>%
              `metadata<-`(
                value = list(overdispersions=overdispersions)
              )
          )
      else GRanges(
        seqnames(chic.experiment.granges.diameter_60),
        ranges(chic.experiment.granges.diameter_60),
        score = (
          as.matrix(elementMetadata(chic.experiment.granges.diameter_60))
          / exp(as.matrix(elementMetadata(chic.experiment.granges.offset_diameter_60)))
        ) %>%
          replace(. < 1e-8, 0) %>%
          as.data.frame
      )
    ),

    tar_map(
      tibble(
        bw = c(25, 100000),
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
      tribble(
        ~filename, ~score, ~score_smooth,
        "Rough_Input",
        list(elementMetadata(chic.experiment.quantify)[, 1]),
        list(elementMetadata(chic.experiment.quantify.smooth_bw100000)[, 1]),
        "Rough_Mark",
        list(elementMetadata(chic.experiment.quantify)[, 2]),
        list(elementMetadata(chic.experiment.quantify.smooth_bw100000)[, 2]),
        "FSeq_Input",
        list(elementMetadata(chic.experiment.quantify.smooth_bw25)[, 1]),
        list(elementMetadata(chic.experiment.quantify.smooth_bw100000)[, 1]),
        "FSeq_Mark",
        list(elementMetadata(chic.experiment.quantify.smooth_bw25)[, 2]),
        list(elementMetadata(chic.experiment.quantify.smooth_bw100000)[, 2]),
        "FSeq_Enrich",
        list(elementMetadata(chic.experiment.quantify.smooth_bw25)[, 2]
        / elementMetadata(chic.experiment.quantify.smooth_bw25)[, 1]),
        list(elementMetadata(chic.experiment.quantify.smooth_bw100000)[, 2]
        / elementMetadata(chic.experiment.quantify.smooth_bw100000)[, 1])
      ) %>%
        rowwise %>%
        summarise(
          filename = with(
            list(mark=mark, name=name, bp_suffix=bp_suffix, reference=reference),
            str_glue("chic/{reference}/{name}_{mark}_{bp_suffix}_{filename}.bw")
          ),
          dir.create = dir.create(dirname(filename), rec=T, showW=F),
          gr = GRanges(
            seqnames(chic.tile.diameter_50_score),
            ranges(chic.tile.diameter_50_score),
            seqlengths = seqlengths(chic.tile.diameter_50_score),
            score = score[[1]] %>%
              replace(
                elementMetadata(chic.experiment.quantify.smooth_bw50)[, 1] < 1,
                score_smooth[[1]][
                  elementMetadata(chic.experiment.quantify.smooth_bw50)[, 1] < 1
                ]
              ) %>%
              `[`(chic.tile.diameter_50_score_lookup) %>%
              replace_na(0)
          ) %>%
            list,
          filename = export(gr, BigWigFile(filename)) %>% as.character
        ) %>%
        pull(filename),
      packages = tar_option_get("packages") %>% c("tidyr")
    ),

    tar_target(
      chic.experiment.model.data,
      tibble(
        molecule = factor(molecule),
        experimental_group = factor(experimental_group) %>%
          fct_relevel(mark),
        rep = factor(as.character(rep)) %>%
          `contrasts<-`(value = contr.helmert(length(levels(.)))),
        # Let's use one window in the exon of "tj" as the example locus for fitting an lm.
        score = chic.experiment.granges.fpkm@elementMetadata[
          which.max(chic.experiment.granges.fpkm@seqnames == "2L" & chic.experiment.granges.fpkm@ranges@start == 19465001),
        ] %>%
          unlist,
        est_library_size = chic.experiment.granges.fpkm@metadata$est_library_size
      ) %>%
        as.list %>%
        with(
          append(
            .,
            {
              experimental_group_rep <- interaction(experimental_group, rep)
              contrasts(experimental_group_rep) <- contr.helmert(length(levels(experimental_group_rep)))
              list(experimental_group_rep = experimental_group_rep)
            }
          )
        )
    ),
    tar_target(
      chic.experiment.lm.sample,
      with(
        chic.experiment.model.data,
        lm(
          "score ~ 0 + molecule + experimental_group_rep",
          list(score, molecule, experimental_group_rep),
          weights = 1 / sqrt(chic.experiment.model.data$est_library_size),
          subset = chic.experiment.model.data$experimental_group == marks
        )
      )
    ),
    tar_target(
      chic.experiment.enrich,
      with(
        manova(
          score ~ 0 + molecule + experimental_group_rep,
          chic.experiment.model.data %>%
            replace(
              names(.) == "score",
              list(t(as.matrix(chic.experiment.granges.fpkm@elementMetadata)))
            ),
          weights = 1 / sqrt(chic.experiment.model.data$est_library_size)
        ),
        {
          score <- replace(
            coefficients[2, ] / coefficients[1, ],
            coefficients[1, ] < 1,
            NA
          )
          GRanges(
            chic.experiment.granges.fpkm@seqnames,
            chic.experiment.granges.fpkm@ranges,
            seqlengths = seqlengths(chic.experiment.granges.fpkm),
            score = score
          )
        }
      )
    )
  ),

  apply(
    chic.experiments %>%
      mutate(driver = driver %>% replace(. == "Nos", "nos")),
    1,
    \(v) tar_target_raw(
      str_glue("chic.regression.result_{v['mark']}_{v['name']}"),
      call(
        "chic_perform_regression",
        target_chic_load_raw_expr(group = v["mark"], driver = v["driver"])
      )
    )
  ),
  tar_map(
    chic.experiments %>%
      rowwise %>%
      mutate(
        driver = driver %>% replace(. == "Nos", "nos"),
        chic.regression.result = rlang::syms(str_glue("chic.regression.result_{mark}_{name}")),
        chic_samples_df = target_chic_load_raw_expr(group = mark, driver = driver) %>%
          list
      ),
    names = mark | name,
    tar_file(
      chic.bw.2,
      tibble(
        data = list(
          chic_samples_df %>%
            rowwise %>%
            mutate(
              bw = if (molecule == "H3") 250 else 25,
              density = smooth_sparse_vector_to_density(pileup, feature.rle, bw=bw) %>% list
            )
        ),
        mymanova = chic_perform_regression(data[[1]]) %>% list,
        myfoldchange = chic_regression_fold_change(mymanova[[1]], data[[1]]$density[[1]]) %>% list,
        mywindow = feature_lengths_to_sliding_windows(feature.lengths, flybase.lengths, 10L, 10L) %>% list,
        mytrack = track_kde_to_sliding(myfoldchange[[1]], mywindow[[1]]) %>% list,
        myt2 = coverage_interp_obs_granges(
          mytrack[[1]],
          replace_na(mytrack[[1]]$coverage >= 1, FALSE)
          &
          is.finite(mytrack[[1]]$fold_change)
        ) %>%
          list,
        output_path = paste0("chic/", driver, "_", mark, ".new.FE.bw"),
        do_export = rtracklayer::export(
          {
            seqlengths(myt2[[1]]) <- flybase.lengths
            myt2[[1]]$fold_change[is.na(myt2[[1]]$fold_change) & !as.logical(myt2[[1]]@seqnames %in% names(chr.lengths))] <- 0
            myt2[[1]]$fold_change <- pmax(-0.1, myt2[[1]]$fold_change)
            GRanges(myt2[[1]], score = myt2[[1]]$fold_change)
          },
          output_path,
          "bigwig"
        ) %>%
          list
      ) %>%
        pull(output_path),
      packages = tar_option_get("packages") %>% c("tidyr")
    )
  ),

  # Profiles of ChIC faceted by quantification
  tar_map(
    experiment.driver %>% mutate(quartile.factor = rlang::syms(str_glue("quartile.factor_{celltype}"))),
    names = celltype,
    tar_target(
      sc_chr_factor,
      tibble(
        rowname = names(quartile.factor),
        chr = read.csv(assay.data.sc)$chr %>%
          fct_recode(`2`="2L", `2`="2R", `3`="3L", `3`="3R") %>%
          factor(c("2", "3", "4", "X", "Y")),
        quant = quartile.factor %>%
          fct_recode(off="Q1", low="Q2", med="Q3", high="Q4"),
        activity = quartile.factor %>%
          fct_recode(off="Q1", active="Q2", active="Q3", active="Q4"),
        fct = interaction(chr, quant),
        fct_activity = interaction(chr, activity)
      )
    ),
    tar_target(
      tss_sc_chr_quartile_data,
      chic_average_profiles(
        pull(sc_chr_factor, fct, rowname),
        dirname(
          c(
            chic.bw_H3K4_Germline,
            chic.bw_H3K27_Germline,
            chic.bw_H3K9_Germline
          )[1]
        ),
        assay.data.sc,
        driver,
        # "CPM Quartile",
        "",
        # fake colors
        sprintf("#%06d", 10 * seq_along(levels(sc_chr_factor$fct)))
      )$data,
      format = "parquet"
    ),
    tar_target(
      tss_sc_chr_data,
      chic_average_profiles(
        pull(sc_chr_factor, chr, rowname),
        dirname(
          c(
            chic.bw_H3K4_Germline,
            chic.bw_H3K27_Germline,
            chic.bw_H3K9_Germline
          )[1]
        ),
        assay.data.sc,
        driver,
        # "CPM Quartile",
        "",
        # fake colors
        sprintf("#%06d", 10 * seq_along(levels(sc_chr_factor$chr)))
      )$data,
      format = "parquet"
    ),
    tar_target(
      tss_sc_chr_active_data,
      chic_average_profiles(
        pull(sc_chr_factor, fct_activity, rowname),
        dirname(
          c(
            chic.bw_H3K4_Germline,
            chic.bw_H3K27_Germline,
            chic.bw_H3K9_Germline
          )[1]
        ),
        assay.data.sc,
        driver,
        # "CPM Quartile",
        "",
        # fake colors
        sprintf("#%06d", 10 * seq_along(levels(sc_chr_factor$chr)))
      )$data,
      format = "parquet"
    )
  ),

  # Temp targets for loading chic.bw.2 and producing repli graphics
  tar_target(
    chic.track_chr,
    tibble(
      filename = chic.bw.2_H3K4_Germline,
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
      filename = chic.bw.2_H3K4_Germline,
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
  tar_target(chic.results_nos_chr, list(H3K4=chic.bw.2_H3K4_Germline, H3K27=chic.bw.2_H3K27_Germline, H3K9=chic.bw.2_H3K9_Germline)),
  tar_target(chic.results_tj_chr, list(H3K4=chic.bw.2_H3K4_Somatic, H3K27=chic.bw.2_H3K27_Somatic, H3K9=chic.bw.2_H3K9_Somatic)),
  tar_target(chic.results_nos_masked, chic.results_nos_chr),
  tar_target(chic.results_tj_masked, chic.results_tj_chr)
)