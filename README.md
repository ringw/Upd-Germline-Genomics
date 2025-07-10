## Germline Analysis - Complete Bioinformatics and Genomic Characterization Project

The project is a complete and transparent (commented) workspace, for a
reduced-complexity *in vivo* stem cell system: single-cell gene expression
in a few main clusters, low-input ChIC-ChIP, and low-input Repliseq. The
[build system](https://books.ropensci.org/targets/) goes from a sample sheet
(of SRA accessions, after publication) to Bowtie2, Samtools view to data frame,
regression and analysis in R, to the published figures.

### Table of Contents

#### Bioinformatics

[targets-flybase.R](R/targets-flybase.R) - Ingesting FASTA genome/feature files and building Bowtie2 index

[targets-zpileup.R](R/targets-zpileup.R), [pileup.R](R/pileup.R) - ChIC OR ChIP OR Repliseq aligned BAM files -> Apache Parquet with the same bulk reads columns -> 

targets-chic
targets-quantification
targets-repli
targets-sce

### Requirements

Ubuntu 20.04.6 LTS

R 4.3.3

Seurat 5.0.3 from [CRAN tar.gz](https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_5.0.3.tar.gz)

SeuratObject 5.0.1 from cran2deb4ubuntu 5.0.1-1cran1.2004.0
