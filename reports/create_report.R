#!/usr/bin/env Rscript

library(knitr)
library(rmarkdown)
library(argparse)


parser <- ArgumentParser()

parser$add_argument("-t", "--template", required=T, help="Path to Rmarkdown template file.")
parser$add_argument("-p", "--patient-id", required=T, help="Patient ID")
parser$add_argument("-r", "--results", required=T, help="Path to CSV file containing mutation and genotype results for the patient.")
parser$add_argument("-rc", "--cna-results-dir", required=T, help="Path to directory containing CNA results for the patient.")
parser$add_argument("-tt", "--tumor-type", required=T, help="Tumor type")
parser$add_argument("-m", "--metadata", required=T, help="Path to file containing meta data.")
parser$add_argument("-d", "--dmp-id", help="DMP patient ID.")
parser$add_argument("-ds", "--dmp-sample-id", help="DMP sample ID")
parser$add_argument("-dm", "--dmp-maf", help="Path to DMP MAF file")
parser$add_argument("-o", "--output", help="Output file")

args <- parser$parse_args()


treatment_file <- NULL
if (!is.null(args$treatment)) {
  treatment_file <- normalizePath(args$treatment)
}

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
  METADATA=normalizePath(args$metadata)
)

tmp <- tempfile(fileext = ".Rmd")
cat(input_text, file = tmp)

rmarkdown::render(
  tmp,
  output_format = "html_document",
  output_dir = normalizePath(dirname(args$output)),
  output_file=args$output)
