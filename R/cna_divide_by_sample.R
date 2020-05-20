# Executable -----------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(stringr)
  library(dplyr)
})

if (!interactive()) {
  
  parser=ArgumentParser()
  parser$add_argument('-i', '--input', type='character', help='file path to CNA pipeline output')
  parser$add_argument('-o', '--output', type='character', help='output directory to write all sample CNA results')
  args=parser$parse_args()
  
  input_cna_file = args$input
  outputdir = args$output
  
  # input_cna_file = '/juno/work/bergerm1/MSK-ACCESS/ACCESS-Projects/test_access/access_data_analysis/cnas/cna_results.txt'
  # outputdir = '/juno/work/bergerm1/MSK-ACCESS/ACCESS-Projects/test_access/access_data_analysis/cnas/'
  cna_results = fread(input_cna_file)[,sample := gsub('\\.','-',gsub('_mean_cvg','',sample))]
  cna_results
  lapply(unique(cna_results$sample),function(x){
    write.table(cna_results[sample == x],
                paste0(outputdir,'/',x,'_cna_results.txt'),
                sep = '\t',quote = F,row.names = F)
  })
  
}
