#!/usr/bin/env Rscript

library(knitr)
library(rmarkdown)
library(argparse)
library(quarto)


parser <- ArgumentParser()

parser$add_argument("-t", "--template", required=T, help="Path to Rmarkdown template file.")
parser$add_argument("-p", "--patient-id", required=T, help="Patient ID")
parser$add_argument("-r", "--results", required=T, help="Path to CSV file containing mutation and genotype results for the patient.")
parser$add_argument("-rc", "--cna-results-dir", required=T, help="Path to directory containing CNA results for the patient.")
parser$add_argument("-tt", "--tumor-type", required=T, help="Tumor type")
parser$add_argument("-m", "--metadata", required=T, help="Path to file containing meta data for each sample. Should contain a \'cmo_sample_id_plasma\', \'sex\', and \'collection_date\' columns. Can also optionally include a \'timepoint\' column (e.g. for treatment information).")
parser$add_argument("-d", "--dmp-id", help="DMP patient ID (optional).")
parser$add_argument("-ds", "--dmp-sample-id", help="DMP sample ID (optional).")
parser$add_argument("-dm", "--dmp-maf", help="Path to DMP MAF file (optional).")
parser$add_argument("-o", "--output", help="Output file with .html extension")
parser$add_argument(
  "-md", "--keep-rmarkdown", help="Dont make tmp file for markdown, keep it in the same directory", action="store_true")
parser$add_argument(
  "-ca", "--combine-access", help="Don't split VAF plots by clonality.", action="store_true")
parser$add_argument(
  "-pi", "--plot-impact", help="Also plot VAFs from IMPACT samples.", action="store_true")

args <- parser$parse_args()

dmp_maf <- NULL
if (!is.null(args$dmp_maf)) {
  dmp_maf <- normalizePath(args$dmp_maf)
}


input_text <- knitr::knit_expand(
  normalizePath(args$template),
  PATIENT_ID=args$patient_id,
  RESULTS=normalizePath(args$results),
  CNA_RESULTS_DIR=normalizePath(args$cna_results_dir),
  TUMOR_TYPE=args$tumor_type,
  DMP_ID=args$dmp_id,
  DMP_SAMPLE_ID=args$dmp_sample_id,
  DMP_MAF_PATH=dmp_maf,
  METADATA=normalizePath(args$metadata),
  COMBINE_ACCESS=args$combine_access,
  PLOT_IMPACT=args$plot_impact
)

tmp <- tempfile(fileext = ".Rmd")
cat(input_text, file = tmp)

if (args$keep_rmarkdown){
  rmd_name <- gsub(".html",".Rmd", args$output)
  output_cwd <- normalizePath(dirname(args$output)) 
  output_rmd_path <- paste(output_cwd,"/",rmd_name, sep='')
  file.copy(tmp,output_rmd_path)
}
#rmarkdown::render(
#  tmp,
#  output_format = "html_document",
#  output_dir = normalizePath(dirname(args$output)),
#  output_file=args$output)

quarto::quarto_render(tmp,output_file=args$output)
