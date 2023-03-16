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
This repository uses CellRanger Counts to generate the CellRangers outs directory that is used downstream: You can download it here:
CellRanger: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
```
You can Download CellRanger Software with this command: 
wget -O cellranger-7.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.1.0.tar.gz?Expires=1677814041&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2Nzc4MTQwNDF9fX1dfQ__&Signature=CeQnToHDIRIiiKInKBCYNLXM3TnZOI346o1XXSjTiPCaAO-B4r1kVheGJF3ZkWWZG1ea1DJN9P4kQ9BHzrP--PHPYhatI5gIB8pxD74WoNAxa4iZxiQAYUKRq7I4z58L2jVwgTbddeOWfSSi2atH2zUVVOOAepnmGkH554a-zdBw2wT4CX1SABsFJ9uODCKkYn5wjg~KxOAK2ULT6FAWcf6IJWLM4OKm9Lz~ill576WFYIfl3OMng~tp-MPC-i641I9mb3qB4O1rsYSvI-IrvULbJ~aCUqNqW9pkdBZxfx5RHzv-SNpSau7WYLN613UMyKshMm07W-GFtwZ~XYYU5w__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
```

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
Python Packages can be installed in Python using:
```
pip install scrublet
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
Github pages of dependencies developed internally: \
scRNA_seq Pipeline: https://github.com/bwh-bioinformatics-hub/rna_seq_pipeline_bwh \
H5MANIPULATOR: https://github.com/bwh-bioinformatics-hub/H5MANIPULATOR \
qcreporter: https://github.com/bwh-bioinformatics-hub/qcreporter \

To create conda environment with dependencies install
```
conda env create -f environment.yml 
```
Setup:

Nextflow will parse a .csv file for the sample_ids and the path to the directory the fastq files are in for your project. Please provide at the minimum a sample_id column to the .csv file. \
If working with a .xls/.xlsx file please create a .csv file called samples.csv with a column labeled == 'sample_id' \
example: \
https://github.com/bwh-bioinformatics-hub/nextflow-scRNAseq/blob/main/samplesheet/samplesheet.csv
```
sample_id,fastq_path
P1708_SP093_105,/mnt/data0/projects/biohub/susanne2022/raw_data/sussane_2022_raw_files/fastq_files/
```


