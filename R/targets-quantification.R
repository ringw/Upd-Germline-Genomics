targets.quantification <- list(
  tar_target(
    transcript.ids,
    read_transcript_ids(assay.data.sc, flybase.gtf)
  ),
  tar_file(
    transcript.isoforms.txt,
    tibble(
      filename = "scRNA-seq/quantify_transcript_ids.txt",
      write.table = write.table(
        transcript.ids %>% subset(select=tx),
        filename,
        row.names = F,
        col.names = F,
        quote = F
      )
    ) %>%
      pull(filename)
  ),
  tar_file(
    gene.ids.flybase.txt,
    tibble(
      filename = "scRNA-seq/quantify_gene_ids.txt",
      write.table = write.table(
        transcript.ids %>% subset(select=flybase) %>% subset(!duplicated(.)),
        filename,
        row.names = F,
        col.names = F,
        quote = F
      )
    ) %>%
      pull(filename)
  ),
  tar_map(
    sce.data,
    names = batch,
    tar_target(
      count.exons,
      Read10XBamExons(
        bam_path,
        transcript.ids,
        cells_to_use = read.csv(metadata, row.names=1) %>%
          rename(batch = "batch.") %>%
          subset(batch. == batch) %>%
          rownames %>%
          str_replace(paste0(batch, "_"), "")
      ),
      packages = tar_option_get("packages") %>% c("Rsamtools")
    ),
    tar_target(
      count.exons.transcript.pct,
      rowSums(count.exons != 0)
    ),
    tar_target(
      count.exons.decont.result,
      count.exons %>%
        `colnames<-`(paste0(batch, "_", colnames(.))) %>%
        decontX,
      packages = tar_option_get("packages") %>% c("decontX")
    ),
    tar_file(
      cell.barcodes.txt,
      tibble(
        filename = str_glue("scRNA-seq/cells_picked_{batch}.txt"),
        write.table = write.table(
          data.frame(
            cell = read.csv(metadata, row.names=1) %>%
              rename(batch = "batch.") %>%
              subset(batch. == batch) %>%
              rownames %>%
              str_replace(paste0(batch, "_"), "")
          ),
          filename,
          row.names = F,
          col.names = F,
          quote = F
        )
      ) %>%
        pull(filename)
    )
  ),
  tar_target(
    transcript.ids.pct.cells,
    with(
      read.csv(metadata, row.names=1),
      {
        nos.1 <- count.exons.transcript.pct_nos.1 / sum(batch == "nos.1")
        nos.2 <- count.exons.transcript.pct_nos.2 / sum(batch == "nos.2")
        tj.1 <- count.exons.transcript.pct_tj.1 / sum(batch == "tj.1")
        tj.2 <- count.exons.transcript.pct_tj.2 / sum(batch == "tj.2")
        (nos.1 + nos.2 + tj.1 + tj.2) / 4
      }
    )
  ),
  tar_target(
    transcript.ids.present,
    transcript.ids.pct.cells >= 0.001
  ),
  tar_target(
    Upd_glm_transcripts.offset.matrix,
    assemble_decont_to_glm_offset(
      list(count.exons_nos.1, count.exons_nos.2, count.exons_tj.1, count.exons_tj.2),
      list(
        count.exons.decont.result_nos.1, count.exons.decont.result_nos.2,
        count.exons.decont.result_tj.1, count.exons.decont.result_tj.2
      ),
      transcript.ids.present,
      metadata
    ) %>%
      as.matrix
  ),
  tar_target(
    Upd_glm_transcripts,
    fit_glm_with_offset(
      combine_counts_mats(
        list(nos.1=count.exons_nos.1, nos.2=count.exons_nos.2, tj.1=count.exons_tj.1, tj.2=count.exons_tj.2)
      )[transcript.ids.present, rownames(read.csv(metadata, row.names=1))],
      Upd_glm_transcripts.offset.matrix,
      Upd_model_matrix,
      metadata
    )
  ),
  tar_target(
    Upd_count_transcripts,
    right_join(
      transcript.ids,
      Upd_glm_transcripts$Beta %>%
        tcrossprod(Upd_celltype_contrasts_with_batch_effect) %>%
        exp %>%
        as.data.frame %>%
        rownames_to_column("tx"),
      "tx"
    )
  )
)