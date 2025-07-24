## Germline Analysis - Complete Bioinformatics and Genomic Characterization Project

The project is a complete and transparent (commented) workspace, for a
reduced-complexity *in vivo* stem cell system: single-cell gene expression
in a few main clusters, low-input ChIC-ChIP, and low-input Repliseq. The
[build system](https://books.ropensci.org/targets/) goes from a sample sheet
(of SRA accessions, after publication) to Bowtie2, Samtools view to data frame,
regression and analysis in R, to the published figures.

### Table of Contents

#### Targets

[targets-flybase.R](R/targets-flybase.R) - Ingesting FASTA genome/feature files and building Bowtie2 index

[targets-zpileup.R](R/targets-zpileup.R), [pileup.R](R/pileup.R) - ChIC OR ChIP OR Repliseq aligned BAM files -> SAM columns plus tags Apache Beam -> new column of fragment midpoints

[targets-chic](R/targets-chic.R) - Defines the sliding windows and all ChIC-seq tracks and figures

[targets-quantification](R/targets-quantification.R) - Defines the objects that we created from the Seurat object for single-cell figures, applying DecontX, 10X Genomics TX tag, and Limma

[targets-repli](R/targets-repli.R) - Defines all Repliseq analysis, tracks, posterior distribution, testing, and figures

[targets-sce](R/targets-sce.R) - Defines the single-cell figures

#### Helper Functions

[bulk-excel](R/bulk-excel.R) - Write some large supplemental tables of genes.

[chic-enrichment](R/chic-enrichment.R) - Wrapper for GLM for the ChIC-ChIP experiment.

[chic-heatmap](R/chic-heatmap.R) - Read into the sliding window GLM coefficient and write to a heatmap for genes.

[chic-lineplot](R/chic-lineplot.R) - Read into the heatmap matrix for genes and write arithmetic mean line plot.

[chic-tracks](R/chic-tracks.R) - Impute the GLM Coefficient with a 2 kb kernel smoothed version (to overcome lack of input in some regions). This writes the BW files to use for a genome browser.

[granges-regression](R/granges-regression.R) - Our GRanges (with score) to SummarizedExperiment converter for Repliseq

[repli-beta](R/repli-beta.R) - Inference for Repliseq with Dirichlet-Multinomial Regression, predictor has Logistic prior and Logistic link function, and responses parameterized by the regularized incomplete Beta function.

[repli-chic-heatmap](R/repli-chic-heatmap.R) - Render colorbars by taking thousands of bins of the genome for each possible RT value and the chromatin landscape.

[repli-logistic](R/repli-logistic.R) - Tanh-like logistic function for creating RT value scores.

[sce-analysis](R/sce-analysis.R) - Identifying DEGs visually (feature plot), total variance (PCA of GSC-like and CySC-like cells), and transcriptome landscape (volcano plot).

[sce-deg](R/sce-deg.R) - Fit our apeglm model for L2FC coefficient.

[sce-excel](R/sce-excel.R) - Write our Limma-fitting and apeglm-fitting details from the single-cell experiment into an Excel report.

[sce-quantification](R/sce-quantification.R) - Recapitulate Cell Ranger BAM tag (this time TX) to a feature matrix.

[sce](R/sce.R) - Identify the clusters in the single-cell experiment.

[smooth-tile-track](R/smooth-tile-track.R) - F-Seq implementation (applied to our ChIC-ChIP regression predictions).

### Requirements

Ubuntu 20.04.6 LTS

R 4.3.3

anndata 0.7.5.6

apeglm 1.24.0

arrow 15.0.1

decontX 1.0.0

dplyr 1.1.4

extraDistr 1.10.0

glmGamPoi 1.14.3

limma 3.56.2

Matrix 1.6-3

mvtnorm 1.2-4

openxlsx 4.2.5.2

pracma 2.4.4

rtracklayer 1.62.0

scran 1.30.2

scuttle 1.12.0

Seurat 5.0.3

SeuratObject 5.0.1

tibble 3.2.1
