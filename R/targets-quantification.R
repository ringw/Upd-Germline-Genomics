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
  tar_file(
    cell_ranger_bam_feature_matrix_coords,
    "scripts/cell_ranger_bam_feature_matrix_coords.pl"
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
      count.exons.genes,
      tibble(
        coords_file = tempfile(),
        do_run_coords = run(
          "sh",
          c(
            "-c",
            paste0(
              "rclone cat sharepoint:Documents/Upd-Germline/",
              bam_path,
              " | ",
              "samtools view -@ 4 -D CB:",
              cell.barcodes.txt,
              " | perl ",
              cell_ranger_bam_feature_matrix_coords,
              " ",
              cell.barcodes.txt,
              " ",
              gene.ids.flybase.txt,
              " GX > ",
              coords_file
            )
          )
        ) %>%
          list,
        my_dimnames = list(
          list(
            read.table(gene.ids.flybase.txt, header=F) %>% pull(1),
            read.table(cell.barcodes.txt, header=F) %>% pull(1)
          )
        ),
        mm_file = coords_to_matrix_market(
          coords_file,
          my_dimnames[[1]],
          nrow(read.table(coords_file, header=F)),
          tempfile(fileext = ".mtx")
        ),
        mymatrix = list(readMM(mm_file))
      ) %>%
        with(
          {
            dimnames(mymatrix[[1]]) <- my_dimnames[[1]]
            mymatrix[[1]]
          }
        )
    ),
    tar_target(
      count.exons.decont.result,
      count.exons %>%
        `colnames<-`(paste0(batch, "_", colnames(.))) %>%
        decontX,
      packages = tar_option_get("packages") %>% c("decontX")
    ),
    tar_target(
      count.exons.genes.decont.result,
      count.exons.genes %>%
        `colnames<-`(paste0(batch, "_", colnames(.))) %>%
        decontX,
      packages = tar_option_get("packages") %>% c("decontX")
    ),
    tar_file(
      cell.barcodes.txt,
      tibble(
        filename = paste0("scRNA-seq/cells_picked_", batch, ".txt"),
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
    Upd_decontX_exons,
    create_decontx_assay_fbgn_objects(
      list(count.exons.genes.decont.result_nos.1, count.exons.genes.decont.result_nos.2,
      count.exons.genes.decont.result_tj.1, count.exons.genes.decont.result_tj.2),
      metadata,
      assay.data.sc
    )
  ),
  tar_target(
    Upd_decontX,
    Upd_sc %>% seurat_decontx_using_batches,
    packages = tar_option_get("packages") %>% c("decontX")
  ),
  tar_target(
    Upd_sc_pseudobulk_transcripts,
    cbind_pseudobulk_decontaminated(
      list(nos.1=count.exons_nos.1, nos.2=count.exons_nos.2, tj.1=count.exons_tj.1, tj.2=count.exons_tj.2),
      list(nos.1=count.exons.decont.result_nos.1[["decontXcounts"]], nos.2=count.exons.decont.result_nos.2[["decontXcounts"]], tj.1=count.exons.decont.result_tj.1[["decontXcounts"]], tj.2=count.exons.decont.result_tj.2[["decontXcounts"]]),
      metadata
    )
  ),
  tar_target(
    Upd_glm_transcripts,
    fit_glm_decontX(Upd_sc_pseudobulk_transcripts, ~ 0 + ident)
  ),
  tar_target(
    Upd_count_transcripts,
    right_join(
      transcript.ids,
      Upd_glm_transcripts$Beta[
        ,
        paste0("ident", sce.clusters$cluster)
      ] %>%
        exp %>%
        as.data.frame %>%
        rownames_to_column("tx"),
      "tx"
    ) %>%
      dplyr::rename(
        gene_id=flybase,
        germline=identgermline,
        somatic=identsomatic,
        spermatocyte=identspermatocyte,
        somaticprecursor=identsomaticprecursor,
        muscle=identmuscle
      ) %>%
      column_to_rownames("tx") %>%
      replace(is.na(.), 0)
  )
)