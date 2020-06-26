# library(data.table)
# library(stringr)
# library(tidyr)
# library(dplyr)

#' @export
CNA_processing = function(
  master.ref,results.dir,
  dmp.dir = '/juno/work/access/production/resources/cbioportal/current/mskimpact'
){
  # # test input section -----------------------------------------------------------
  # master.ref = fread('/juno/work/bergerm1/bergerlab/zhengy1/access_data_analysis/data/example_master_file.csv')
  # results.dir = paste0('/juno/work/bergerm1/MSK-ACCESS/ACCESS-Projects/test_access/access_data_analysis/output_042020/')
  # dmp.dir = '/ifs/work/bergerm1/zhengy1/dmp-2020/mskimpact/'
  # 
  # DMP CNA calls --------------------------------------------------------
  DMP.cna <- fread(paste0(dmp.dir,'data_CNA.txt'))
  
  DMP.cna.edited = DMP.cna[!Hugo_Symbol %in% c('CDKN2Ap14ARF','CDKN2Ap16INK4A'),c(
    'Hugo_Symbol',grep(paste0(master.ref$dmp_patient_id,collapse = '|'),colnames(DMP.cna),value = T)
  ),with = F] %>% melt.data.table(id.vars = 'Hugo_Symbol',variable.name = 'Tumor_Sample_Barcode',value.name = 'fc') %>%
    filter(fc != 0) %>% mutate(CNA_tumor = case_when(fc < 0 ~ 'HOMDEL',fc > 0 ~ 'AMP'),
                               dmp_patient_id = gsub('-T..-IM.','',Tumor_Sample_Barcode)) %>% 
    merge(unique(master.ref[!is.na(dmp_patient_id),.(cmo_patient_id,dmp_patient_id)]),by = 'dmp_patient_id',all.x = T) %>%
    transmute(Hugo_Symbol,CNA_tumor,cmo_patient_id,dmp_patient_id) %>% unique() %>% data.table()
  
  
  # access CNA calls --------------------------------------------------------
  access.cna = do.call(rbind,lapply(master.ref$cna_path,function(x){
    fread(x) %>% transmute(Tumor_Sample_Barcode = gsub('\\.','-',gsub('_mean_cvg','',sample)),Hugo_Symbol = region, p.adj,fc) %>% 
      filter(p.adj <= 0.05) %>% merge(unique(master.ref[,.(cmo_patient_id,cmo_sample_id_plasma)]),
                                      by.x = 'Tumor_Sample_Barcode',by.y = 'cmo_sample_id_plasma',all.x = T) %>% data.table()
  }))
  
  
  access.cna.genes = c('AR','BRCA1','BRCA2','CDK4','EGFR','ERBB2','MET','MDM2','MLH1','MSH2','MSH6','MYC')
  merge(access.cna,DMP.cna.edited,by = c('Hugo_Symbol','cmo_patient_id'),all.x = T) %>% mutate(CNA = case_when(
    (fc > 1.5 & Hugo_Symbol %in% access.cna.genes) ~ 'AMP',
    (fc < -1.5 & Hugo_Symbol %in% access.cna.genes) ~ 'HOMDEL',
    # & Hugo_Symbol %in% access.cna.genes
    (fc > 1.2 & !is.na(CNA_tumor)) ~ 'AMP',
    (fc > 1.2 & !is.na(CNA_tumor)) ~ 'HOMDEL',
    TRUE ~ ''
  )) %>% filter(CNA != '') %>% data.table() -> access.call.set
  
  if(any(access.final.call.set$CNA_tumor != access.final.call.set$CNA))  warning('Tumor and plasma have conflicting CNA detected!')
  
  dir.create(paste0(results.dir,'/CNA_final_call_set'))
  lapply(unique(master.ref$cmo_sample_id_plasma),function(x){
    write.table(access.call.set[Tumor_Sample_Barcode == x],
                paste0(results.dir,'/CNA_final_call_set/',x,'_cna_final_call_set.txt'),
                sep = '\t',quote = F,row.names = F)
  })
  
}

# Executable -----------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(stringr)
  library(dplyr)
})

if (!interactive()) {
  
  parser=ArgumentParser()
  parser$add_argument('-m', '--masterref', type='character', help='File path to master reference file')
  parser$add_argument('-o', '--resultsdir', type='character', help='Output directory')
  parser$add_argument('-dmp', '--dmpdir', type='character', default = '/juno/work/access/production/resources/cbioportal/current/mskimpact',
                      help='Directory of clinical DMP IMPACT repository [default]')
  args=parser$parse_args()
  
  master.ref = args$masterref
  results.dir = args$resultsdir
  dmpdir = args$dmpdir
  
  CNA_processing(master.ref,results.dir,dmp.dir)  
}
