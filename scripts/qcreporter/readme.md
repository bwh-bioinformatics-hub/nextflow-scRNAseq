# RNA/ATAC Seq QC Reporter

R functions and scripts for summarizing and performing sample-level QC on sequencing pipeline data.  

<a id="contents"></a>

## Contents

#### [Dependencies](#dependencies)

#### [Installation](#installation)

#### [Available Reports](#available_report)

#### [Running Reports](#batch_report)


#### [NGS Batch Report](#ngs_batch_report)
- [Parameters](#ngs_report_param)
- [Sample Sheet Guidelines](#ngs_sample_sheet)
- [Outputs](#ngs_report_out)
- [Tests](#ngs_report_test)  

#### [Version History](#version_history)

<a id="dependencies"></a>

## Dependencies    

This repository requires that `pandoc` and `libhdf5-devel` libraries are installed as dependencies of the `H5MANIPULATOR` functions:
```
sudo apt-get install pandoc libhdf5-dev
```

CRAN packages can be installed in R using:
```
# Example
install.packages("jsonlite")
install.packages("rmarkdown")
install.packages("optparse")
```

`H5MANIPULATOR`, Because it is a private repository, you may need to provided a [Github Personal Access Token](https://github.com/settings/tokens) for installation:
```
Sys.setenv(GITHUB_PAT = "[your_PAT_here]")
devtools::install_github("acicalo2/H5MANIPULATOR")
```

[Return to Contents](#contents)

<a id="installation"></a>

## Installation
`qcreporter` is an R package with associated executable scripts. First install dependencies. Then install the R package from Github repository:

```
Sys.setenv(GITHUB_PAT = "[your_PAT_here]")
devtools::install_github("acicalo2/qcreporter")
```  
To run scripts, clone the GitHub repository and run the desired wrapper script within the local clone.

[Return to Contents](#contents)  

<a id="available_report"></a>

## Available Reports and Scripts 
Available batch reports are as follows:
- [scRNA cell report](#scrna_batch_report): Batch report for stand alone scRNA 

[Return to Contents](#contents)  

<a id="scrna_report_param"></a>

#### Input Parameters

There are 9 parameters for this script:  

* `-b or --experiment`:  The experiment name, ie B001
* `-i or --in_dir`: The input directory containing the files to process. Should include the following subdirectories:  
  * `h5`: Contains all sample .h5 files for the batch generated by [`tenx_rna_metadata_update`](https://github.com/acicalo2/rna_seq_pipeline_bwh)  
* `-k or --in_key`: A 6-column .csv Sample Sheet (see format [below](#scrna_sample_sheet)) of identifiers for all samples in the batch.
* `-n or --n_cores`: An integer value of number of cores to use for multithreaded processes, used by Seurat functions  
* `-m or --mc_mb_limit`: An integer value of number of maximum size in Mb allowed for exporting globals to each worker in multicore processing (for futures R package). Defaults to 50000, suggest >= 20000.  
* `-d or --out_dir`: A directory path to use to output the batch report.
* `-o or --out_html`: A filename to use to output the HTML summary report file. For example "B001_scRNA_batch_report.html"

[Return to Contents](#contents)  


[Return to Contents](#contents)

<a id="ngs_batch_report"></a>

## RNA-Seq Sample Report

The `qc_batch_summary.r` wrapper script renders `qc_report_rna_seq.rmd to create a RNA-Seq Sample QC report (html). 

[Return to Contents](#contents)  

<a id="rna_seq_report_param"></a>

#### Input Parameters

There are 8 parameters for this script:  

* `-b or --batch_id`:  The experiment name,
* `-m or --in_method`:  A ";"-delimited string of the data streams being 
processed enclosed in quotes, for example "scrna" 
  * `scrna` (single cell RNA)
* `-i or --in_dir`: The input directory containing all results files to process. Must 
include one subdirectory per modality named in `in_method` above.
* `-k or --in_key`: 
* `-d or --out_dir`: A directory path to use to output the batch report.
* `-o or --out_html`: A filename to use to output the HTML summary report file. 
For example "experiment_id_scRNA_report.html"

An example run:
```
git clone https://github.com/acicalo2/rna_atac_seq_qc_reporter.git

Rscript --vanilla \
    /home/jupyter/batchreporter/qc_batch_summary.r \
    -e experiment_id  \
    -m 'scrna;scatac' \
    -i /home/acicalo/output/scrna \
    -k /home/acicalo/samplesheet/samplesheet.csv   \
    -d /home/acicalo/output/ \
    -o  experiment_id_rnaseq_sample_report.html 
    ```

[Return to Contents](#contents)

<a id="rnaseq_report_out"></a>

#### Output Files

`qc_batch_summary.r` will generate the HTML reporting file with name as defined by input parameter -o. 
```
experiment_id_rnaseq_sample_report.html 

```

[Return to Contents](#contents)


#### Version History  

|Version|Date|Description|
|--------|---------|------------------------------------------------|
|1.1.0|01/12/2023|Addition of sample report|
|1.0.0 |          |Initial release, scRNA |  


[Return to Contents](#contents)