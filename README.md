# nextflow-scRNAseq
Nextflow for single-cell RNAseq

### [Version History](#version_history)

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
```
```
To download dependencies that were developed internally by BWH Bioinformatics and Genomics Hub
```
Sys.setenv(GITHUB_PAT = "[your_PAT_here]")
devtools::install_github("bwh-bioinformatics-hub/H5MANIPULATOR")
devtools::install_github("bwh-bioinformatics-hub/H5MANIPULATOR")

git clone https://github.com/bwh-bioinformatics-hub/rna_seq_pipeline_bwh.git
