sc_cells_dendrogram <- function(seurat) {
  # do not subsample cells
  dendro_cells <- Cells(seurat)
  seurat <- seurat %>% NormalizeData
  dendro_pca <- GetAssayData(seurat)[
    sapply(seurat[["RNA"]]@meta.data$pass.min.cells, isTRUE),
  ] %>%
    as.matrix %>%
    RunPCA(verb = F, assay = "RNA", npcs = 25)
  dendro_obj <- dendro_pca@cell.embeddings[dendro_cells, ] %>%
    dist %>%
    hclust(method = "average") %>%
    as.dendrogram
}

sc_sort_cells_dendrogram <- function(seurat, dendro_obj) {
  # wts <- as.numeric(Idents(seurat))
  # dendro_obj <- dendro_obj %>% reorder(wts, agglo.FUN=mean)
  order_cells <- FetchData(seurat, c("ident", "nCount_RNA")) %>%
    rownames_to_column %>%
    group_by(ident) %>%
    reframe(cell = rowname[order(nCount_RNA, decreasing=T)]) %>%
    pull(cell)
  dendro_obj <- dendro_obj %>% rotate(order_cells)
}

sc_genes_dendrogram <- function(seurat, gene_logical = sapply(seurat[["RNA"]]@meta.data$pass.min.cells, isTRUE)) {
  seurat <- seurat %>% NormalizeData
  dendro_scale_data <- GetAssayData(seurat)[
    gene_logical,
  ] %>%
    ScaleData(verb=F)
  dendro_dist <- (dendro_scale_data %*% t(dendro_scale_data)) %>%
    `-`(max(.)) %>%
    `*`(-1 / sqrt(ncol(seurat))) %>%
    as.dist
  dendro_obj <- dendro_dist %>% hclust %>% as.dendrogram
}

sc_sort_genes_dendrogram <- function(Upd_cpm, dendro_obj) {
  order_genes <- apply(
    Upd_cpm[labels(dendro_obj), ],
    1,
    \(v) data.frame(index = which.max(v), value = max(v))
  ) %>%
    bind_rows(.id = "gene") %>%
    group_by(index) %>%
    reframe(gene = gene[order(value, decreasing=T)]) %>%
    pull(gene)
  dendro_obj <- dendro_obj %>% rotate(order_genes)
}

plot_single_cell_heatmap <- function(seurat, genes, cells, cells_cut = 10, genes_cut = 130) {
  seurat <- seurat %>% NormalizeData
  genes_plot <- labels(genes)
  cells_plot <- labels(cells)

  # Hang dendrogram - bring leaf tree nodes up. The leaves (0) may initially be
  # very far from the parent node (which tracks the avg distance).
  cells <- (cells %>% cut(h=cells_cut))$upper %>% hang.dendrogram
  genes <- (genes %>% cut(h=genes_cut))$upper %>% hang.dendrogram

  # Indicator bar for cell idents. How tall should it be as a percentage of the
  # height of the plot? The height is decided by the genes.
  scalbar <- 0.1 * length(genes_plot)
  cells_data <- dendro_data(cells)$segments
  scal <- -0.3 * length(genes_plot) / diff(range(c(cells_data$y, cells_data$yend)))
  cells_data$y <- scal * (cells_data$y - min(cells_data$yend))
  cells_data$yend <- scal * (cells_data$yend - min(cells_data$yend))
  cells_data <- with(cells_data, data.frame(x=x, y=y, xend=xend, yend=yend)) %>%
    mutate(y = y - scalbar, yend = yend - scalbar)

  genes_data <- dendro_data(genes)$segments
  scal <- -0.2 * length(cells_plot) / diff(range(c(genes_data$y, genes_data$yend)))
  genes_data$y <- scal * (genes_data$y - min(genes_data$yend))
  genes_data$yend <- scal * (genes_data$yend - min(genes_data$yend))
  genes_data <- with(genes_data, data.frame(x=y, y=x, xend=yend, yend=xend))
  all_genes <- GetAssayData(seurat)[genes_plot, cells_plot] %>%
    as.matrix %>%
    resizeImage(2000, floor(length(genes_plot) / 2), "bilinear")
  dimnames(all_genes) <- list(
    # genes
    y=NULL,
    # cells
    x=NULL
  )
  all_genes <- melt(all_genes, value.name="LogNormalize") %>%
    mutate(
      y = seq(1, length(genes_plot), length.out=nrow(all_genes))[y],
      x = seq(1, length(cells_plot), length.out=ncol(all_genes))[x]
    )
  gg <- (
    ggplot(
      all_genes,
      aes(x, y)
    )
    + geom_raster(aes(fill=LogNormalize))
    + scale_fill_viridis_c(limits = c(0, 3), oob=squish)
    + new_scale_fill()
    + geom_tile(
      aes(y = -scalbar/2, height = scalbar, fill = ident),
      mutate(FetchData(seurat, "ident")[cells_plot,, drop=F], x = seq_along(ident))
    )
    + scale_fill_manual(values = cluster_colors, guide=guide_none())
    
    # Show each dendro:
    + geom_segment(
      aes(xend=xend, yend=yend),
      rbind(cells_data, genes_data),
      linewidth = 0.2
    )
    + scale_y_continuous(pos = "right")
    + coord_cartesian(
      c(0.5, length(cells_plot) + 0.5),
      c(length(genes_plot) + 0.5, 0.5),
      expand = F, clip = "off"
    )
    + labs(
      x = "Cells",
      y = paste0("Genes")
    )
    + theme(
      plot.margin = unit(c(4.5, 0.2, 0.2, 3), "cm")
    )
  )
}