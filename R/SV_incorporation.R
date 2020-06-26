# library(data.table)
# library(stringr)
# library(tidyr)
# library(dplyr)

#' @export
SV_incorporation = function(
  master.ref,results.dir,
  dmp.dir = '/juno/work/access/production/resources/cbioportal/current/mskimpact',
  criteria = 'stringent'
){
  # # test input section -----------------------------------------------------------
  # master.ref = fread('/juno/work/bergerm1/bergerlab/zhengy1/access_data_analysis/data/example_master_file.csv')
  # results.dir = paste0('/juno/work/bergerm1/MSK-ACCESS/ACCESS-Projects/test_access/access_data_analysis/output_042020/')
  # dmp.dir = '/ifs/work/bergerm1/zhengy1/dmp-2020/mskimpact/'
  # # criteria <- 'permissive'
  # criteria <- 'stringent'
  # 
  # DMP fusion calls --------------------------------------------------------
  DMP.fusion <- fread(paste0(dmp.dir,'/data_SV.txt')) %>%
    transmute(DMP_SAMPLE_ID = SampleId,EventType = Sv_Class_Name,Gene1 = Site1_Gene,Gene2 = Site2_Gene,
              Chr1 = Site1_Chrom,Chr2 = Site2_Chrom,Pos1 = Site1_Pos,Pos2 = Site2_Pos,PairedReadCount = Paired_End_Read_Support,
              SplitReadCount = Split_Read_Support,TumorReadCount = Tumor_Read_Count,EventInfo = Annotation) %>% data.table()
  
  # execution ---------------------------------------------------------------
  
  dir.create(paste0(results.dir,'/results_',criteria,'_combined'))
  gene.list <- c('ALK','BRAF','ERBB2','EGFR','FGFR1','FGFR2','FGFR3','KIT','MET','NTRK1','NTRK2','NTRK3','PDGFRA','RET','ROS1')
  # x <- unique(master.ref$cmo_patient_id)[1]
  lapply(unique(master.ref$cmo_patient_id),function(x){
    # get sample sheet --------------------------------------------------------
    sample.sheet <- fread(paste0(results.dir,'/',x,'/',x,'_sample_sheet.tsv'))
    # get plasma SV calls -----------------------------------------------------
    total.sv <- do.call(rbind,lapply(sample.sheet[Sample_Type == 'duplex']$Sample_Barcode,function(y){
      SV.filename <- master.ref[cmo_sample_id_plasma == y]$sv_path
      if(!file.exists(SV.filename)){stop(paste0('SV file: ',SV.filename,' ----- does not exist'))}
      tmp.SV <- fread(SV.filename) %>% filter(Significance == 'KeyGene')
      return(tmp.SV)
    })) %>% transmute(TumorId = paste0(TumorId,'___total'),SV_Type,Gene1,Gene2,Chr1,Chr2,Pos1,Pos2,PairedReadCount = PairEndReadSupport,SplitReadCount = SplitReadSupport,TumorReadCount,Fusion) 
    
    # get DMP SV calls --------------------------------------------------------
    DMP.sv <- do.call(rbind,lapply(sample.sheet[Sample_Type == 'DMP_Tumor']$Sample_Barcode,function(y){
      DMP.fusion[DMP_SAMPLE_ID == y]
    })) 
    if(!is.null(DMP.sv)){
      DMP.sv <- DMP.sv %>% dplyr::transmute(TumorId = paste0(DMP_SAMPLE_ID,'___DMP_Tumor'),SV_Type = case_when(EventType == 'INVERSION' ~ 'INV',EventType == 'DELETION' ~ 'DEL',
                                                                                                               EventType == 'INSERTION' ~ 'INS',EventType == 'DUPLICATION' ~ 'DUP',
                                                                                                               EventType == 'TRANSLOCATION' ~ 'BND',TRUE ~ 'UNKNOWN'),
                                            Gene1,Gene2,Chr1,Chr2,Pos1,Pos2,PairedReadCount,SplitReadCount,TumorReadCount = NA,Fusion = EventInfo)  
    }else{
      # dummy df if there is no DMP fusion found
      DMP.sv <- data.frame(matrix(nrow = 0,ncol = ncol(DMP.fusion)))
      colnames(DMP.sv) <- colnames(DMP.fusion)
    }
    
    # event desc. reconciliating possible DMP vs manta ------------------------
    rbind(total.sv,DMP.sv) %>% 
      # consolidate information for dcasting data frame
      unite(sample_info,PairedReadCount,SplitReadCount,TumorReadCount,sep = '-') %>% unite(gene_pair,Gene1,Gene2,sep = '__') %>% unite(chr_pair,Chr1,Chr2,sep = '__') %>% 
      # selecting all unique events with their description (manta outputs will be the same, within DMP will be the same)
      unite(event_info,gene_pair,chr_pair,SV_Type,Pos1,Pos2,sep = '---') %>% select(event_info,Fusion) %>% unique() %>% data.table() -> event.desc 
    event.desc <- event.desc[,.(HGVSp_Short = paste0(unique(Fusion),collapse = ' | ')),.(event_info)]
    
    # make into a row ---------------------------------------------------------
    rbind(total.sv,DMP.sv) %>% select(-Fusion) %>%
      # consolidate information for dcasting data frame
      unite(sample_info,PairedReadCount,SplitReadCount,TumorReadCount,sep = '-') %>% unite(gene_pair,Gene1,Gene2,sep = '__') %>% unite(chr_pair,Chr1,Chr2,sep = '__') %>% 
      unite(event_info,gene_pair,chr_pair,SV_Type,Pos1,Pos2,sep = '---') %>% spread(TumorId,sample_info) %>% 
      merge(event.desc,by = 'event_info',all.x = T) %>%
      # parse out information
      separate(event_info,into = c('gene_pair','chr_pair','SV_Type','Pos1','Pos2'),sep = '---') %>%
      # process information for row binding into SNV table
      setnames(c('gene_pair','chr_pair','SV_Type','Pos1','Pos2'),c('Hugo_Symbol','Chromosome','Variant_Classification','Start_Position','End_Position')) %>%
      mutate(Reference_Allele = '',Tumor_Seq_Allele2 = '',ExAC_AF = '',Hotspot = ifelse(grepl(paste0(gene.list,collapse = '|'),Hugo_Symbol),'Significance gene',''),
             DMP = '',duplex_support_num = '',call_confidence = '') %>%
      select(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,HGVSp_Short,Reference_Allele,Tumor_Seq_Allele2,ExAC_AF,Hotspot,DMP,duplex_support_num,call_confidence,everything()) %>%
      data.table() -> SV.table
    
    # adding some annotation columns specific to snv table --------------------
    apply(sample.sheet[Sample_Type != 'plasma_simplex'],1,function(y){
      current.colname <- paste0(y[which(colnames(sample.sheet) == 'Sample_Barcode')],'___',
                                case_when(y[which(colnames(sample.sheet) == 'Sample_Type')] == 'duplex' ~'total',
                                          y[which(colnames(sample.sheet) == 'Sample_Type')] == 'DMP_Tumor' ~'DMP_Tumor',
                                          y[which(colnames(sample.sheet) == 'Sample_Type')] == 'DMP_Normal' ~'DMP_Normal',
                                          y[which(colnames(sample.sheet) == 'Sample_Type')] == 'unfilterednormal' ~'unfilterednormal'))
      # If any of the sample does not already have a column
      if(!current.colname %in% colnames(SV.table)){
        SV.table[,eval(current.colname) := NA]
      }
      # only for plasma sample -- '.called' column
      if(y[which(colnames(sample.sheet) == 'Sample_Type')] == 'duplex'){
        SV.table[,paste0(y[which(colnames(sample.sheet) == 'Sample_Barcode')],'___duplex.called') := case_when(
          is.na(get(paste0(y[which(colnames(sample.sheet) == 'Sample_Barcode')],'___total'))) ~ 'Not Called',
          is.na(all(!str_split(Hugo_Symbol,'__') %in% access.gene.list)) ~ 'Not Covered',
          T ~ 'Called'
        )]
      }
      # edit DMP signout column
      if(y[which(colnames(sample.sheet) == 'Sample_Type')] == 'DMP_Tumor'){
        SV.table[!is.na(get(current.colname)),DMP := 'Signed out'] 
      }
    })
    # adding CH filler column
    SV.table$CH = ''
    
    # append and write --------------------------------------------------------
    tmp.snv.table.path = paste0(results.dir,'/results_',criteria,'/',x,'_SNV_table.csv')
    if(!file.exists(tmp.snv.table.path)){
      stop(paste0(tmp.snv.table.path,' does not exist. Check if filter_calls was run correcly'))
    }
    SNV.table <- fread(tmp.snv.table.path)
    snv.sv.table.directory <- paste0(results.dir,'/results_',criteria,'_combined/',x,'_table.csv') 
    write.csv(rbind(SNV.table,SV.table),snv.sv.table.directory,quote = F,row.names = F)
    
  })
  
}

# Executable -----------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(tidyr)
  library(stringr)
  library(dplyr)
  library(argparse)
})

if (!interactive()) {
  
  parser=ArgumentParser()
  parser$add_argument('-m', '--masterref', type='character', help='File path to master reference file')
  parser$add_argument('-o', '--resultsdir', type='character', help='Output directory')
  parser$add_argument('-dmp', '--dmpdir', type='character', default = '/juno/work/access/production/resources/cbioportal/current/mskimpact',
                      help='Directory of clinical DMP IMPACT repository [default]')
  parser$add_argument('-c', '--criteria', type='character', default = 'stringent',
                      help='Calling criteria [default]')
  args=parser$parse_args()
  
  master.ref = args$masterref
  results.dir = args$resultsdir
  dmpdir = args$dmpdir
  criteria = args$criteria
  
  cat(paste0(paste0(c(paste0(rep('-',15),collapse = ''),'Arguments input: ',master.ref,results.dir,dmpdir,criteria,
                      paste0(rep('-',15),collapse = '')),collapse = "\n"),'\n'))
  
  if(!criteria %in% c('stringent','permissive')){
    stop('Criteria argument should be either stringent or permissive')
  }
  
  suppressWarnings(SV_incorporation(fread(master.ref),results.dir,dmpdir,criteria))
  
}

