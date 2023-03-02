# nextflow-scRNAseq
Nextflow for single-cell RNAseq

**************************

# 1. Introduction

### Overview:

This pipeline manages a scRNA-Seq workflow starting from raw fastq files and converting
them to standard file formats for use by downstream tools. The steps involved are:

* Cell Ranger Count takes FASTQ files from cellranger mkfastq and performs alignment, filtering, barcode counting, and UMI counting. It uses the Chromium   cellular barcodes to generate feature-barcode matrices, determine clusters, and perform gene expression analysis. The count pipeline can take input     
  from multiple sequencing runs on the same GEM well. cellranger count also processes Feature Barcode data alongside Gene Expression reads.
* Scrublet Process: Single-Cell Remover of Doublets, Python code for identifying doublets in single-cell RNA-seq data.
* Add meta, will add the metadata from CellRanger Count and user provided samplesheet to .h5 file.
* Create QC Report of all Samples provided.
* Cell Clustering and Cell Type Annotation.
* Perform Trajectory Analysis.
<a id="dependencies"></a>

## Dependencies    

This repository requires that `pandoc` and `libhdf5-devel` libraries are installed as dependencies of the `H5MANIPULATOR` functions:
```
sudo apt-get install pandoc libhdf5-dev
```

CRAN packages can be installed in R using:
```
install.packages("devtools")
install.packages('viridis',repo="https://cloud.r-project.org")
install.packages("scCustomize",repo="https://cloud.r-project.org")
# Dot plot is depedent on GitHub Report (https://github.com/Simon-Leonard/FlexDotPlot)
devtools::install_github("Simon-Leonard/FlexDotPlot")
```
Some Packages are Dependent on BiocManager
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("scry")
BiocManager::install("ComplexHeatmap")
BiocManager::install("rhdf5")

```
To download dependencies that were developed internally by BWH Bioinformatics and Genomics Hub
```
Sys.setenv(GITHUB_PAT = "[your_PAT_here]")
devtools::install_github("bwh-bioinformatics-hub/H5MANIPULATOR")
devtools::install_github("bwh-bioinformatics-hub/qcreporter")

git clone https://github.com/bwh-bioinformatics-hub/rna_seq_pipeline_bwh.git
```
To create conda environment with dependencies install
```
conda env create -f environment.yml 
```
