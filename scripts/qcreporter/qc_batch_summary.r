library(optparse)

option_list <- list(
  make_option(opt_str = c("-e","--experiment_id"),
              type = "character",
              default = NULL,
              help = "Sample identifier",
              metavar = "character"),
  make_option(opt_str = c("-m","--in_method"),
              type = "character",
              default = NULL,
              help = "Input batch pipeline modality string",
              metavar = "character"),
  make_option(opt_str = c("-i","--in_dir"),
              type = "character",
              default = NULL,
              help = "Input directory containing h5 and json files",
              metavar = "character"),
  make_option(opt_str = c("-z","--cellrangers_dir"),
              type = "character",
              default = NULL,
              help = "cellrangers out directory",
              metavar = "character"),
  make_option(opt_str = c("-u","--scrubletdir"),
              type = "character",
              default = NULL,
              help = "scrublet directory",
              metavar = "character"),
  make_option(opt_str = c("-f","--refdir"),
              type = "character",
              default = NULL,
              help = "reference directory",
              metavar = "character"),
  make_option(opt_str = c("-k","--in_key"),
              type = "character",
              default = NULL,
              help = "Input sample sheet",
              metavar = "character"),
  make_option(opt_str = c("-d","--out_dir"),
              type = "character",
              default = NULL,
              help = "Output Directory",
              metavar = "character"),
  make_option(opt_str = c("-o","--out_html"),
              type = "character",
              default = NULL,
              help = "Output HTML run summary file",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)

args <- parse_args(opt_parser)

if(is.null(args$experiment_id)) {
  print_help(opt_parser)
  stop("No parameters supplied.")
}

if(!dir.exists(args$out_dir)) {
  dir.create(args$out_dir)
}

rmd_path <- file.path(args$out_dir,
                      paste0(args$experiment_id,
                             "_ngs_sample_qc_report.rmd"))

file.copy(system.file("rmarkdown/ngs_sample_qc_report.rmd", package = "qcreporter"),
          rmd_path,
          overwrite = TRUE)

rmarkdown::render(
  input = rmd_path,
  params = list(  experiment_id    = args$experiment_id,
                  cellrangers_dir  = args$cellrangers_dir,
                  scrubletdir      = args$scrubletdir,
                  outdir           = args$out_dir,
                  refdir           = args$refdir,
                  projectName      = args$experiment_id,
                  in_method_string = args$in_method,
                  in_dir           = args$in_dir,  
                  in_key           = args$in_key,  
                  out_dir          = args$out_dir),
  output_file = args$out_html,
  quiet = TRUE
)

file.remove(rmd_path)