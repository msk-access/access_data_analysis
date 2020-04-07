library(data.table)
library(stringr)
library(tidyr)
library(dplyr)


master.ref <- fread('/juno/work/bergerm1/bergerlab/zhengy1/TDM1/Sample_mapping/master_ref_031920.csv')
results.dir <- paste0('/juno/work/bergerm1/bergerlab/zhengy1/TDM1/Analysis_files/run_022520/')
manta.dir <- '/juno/work/bergerm1/bergerlab/zhengy1/TDM1/Analysis_files/manta_032220/'
criteria <- 'stringent'
# criteria <- 'permissive'
dir.create(paste0(results.dir,'/results_',criteria,'_combined'))


if(length(list.files(paste0(results.dir,'/results_',criteria,'_reviewed/'))) != length(unique(master.ref$`CMO patient ID`))){
  stop('Not all patients had SNV reported')
}

# DMP fusion calls --------------------------------------------------------
DMP.fusion <- fread('/juno/work/bergerm1/bergerlab/zhengy1/TDM1/Original_files/SVs_Peter.txt')
colnames(DMP.fusion) <- c('Run','DMP_SAMPLE_ID','alys2sample_id','EventType','Chr1','Pos1','Chr2','Pos2','Gene1','Gene2','Gene1desc','Gene2desc',
                          'EventInfo','Breakpoint_Type','PairedReadCount','SplitReadCount','TumorType','Annotation','Comments')

# execution ---------------------------------------------------------------

gene.list <- c('ALK','BRAF','ERBB2','EGFR','FGFR1','FGFR2','FGFR3','KIT','MET','NTRK1','NTRK2','NTRK3','PDGFRA','RET','ROS1')
x <- list.files(paste0(results.dir,'/results_',criteria,'_reviewed/'),full.names = T)[1]
lapply(list.files(paste0(results.dir,'/results_',criteria,'_reviewed/'),full.names = T),function(x){
# get sample sheet --------------------------------------------------------
  TDM1.ID <- gsub(paste0('.*./results_',criteria,'_reviewed/+','|_SNV_table.csv'),'',x)
  print(TDM1.ID)
  CMO.ID <- unique(master.ref[GRAIL_PT_ID == TDM1.ID]$`CMO patient ID`) 
  if(length(CMO.ID) != 1){
    stop('Something wrong with non-unique EDD ID to CMO ID mapping, try to debug')
  }
  sample.sheet <- fread(paste0(results.dir,'/',CMO.ID,'/',CMO.ID,'_sample_sheet.tsv'))
# get plasma SV calls -----------------------------------------------------
  total.sv <- do.call(rbind,lapply(sample.sheet[Sample_Type == 'plasma']$Sample_Barcode,function(y){
    SV.filename <- paste0(manta.dir,'/final_results/',gsub('_IGO.*.','',y),'_Annotated_Evidence-annotated.txt')
    if(!file.exists(SV.filename)){stop(paste0('SV file: ',SV.filename,' ----- does not exist'))}
    tmp.SV <- fread(SV.filename) %>% filter(Significance == 'KeyGene')
    return(tmp.SV)
  })) %>% transmute(TumorId = paste0(TumorId,'___total'),SV_Type,Gene1,Gene2,Chr1,Chr2,Pos1,Pos2,PairedReadCount = PairEndReadSupport,SplitReadCount = SplitReadSupport,TumorReadCount,Fusion) 
  
# get DMP SV calls --------------------------------------------------------
  DMP.sv <- do.call(rbind,lapply(sample.sheet[Sample_Type == 'Tumor']$Sample_Barcode,function(y){
    DMP.fusion[DMP_SAMPLE_ID == y]
  })) 
  if(!is.null(DMP.sv)){
    DMP.sv <- DMP.sv %>% dplyr::transmute(TumorId = paste0(DMP_SAMPLE_ID,'___Tumor'),SV_Type = case_when(EventType == 'INVERSION' ~ 'INV',EventType == 'DELETION' ~ 'DEL',
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
                              case_when(y[which(colnames(sample.sheet) == 'Sample_Type')] == 'plasma' ~'total',
                                        y[which(colnames(sample.sheet) == 'Sample_Type')] == 'Tumor' ~'Tumor',
                                        y[which(colnames(sample.sheet) == 'Sample_Type')] == 'normal_DMP' ~'normal_DMP',
                                        y[which(colnames(sample.sheet) == 'Sample_Type')] == 'normal' ~'unfilterednormal'))
    # If any of the sample does not already have a column
    if(!current.colname %in% colnames(SV.table)){
      SV.table[,eval(current.colname) := NA]
    }
    # only for plasma sample -- '.called' column
    if(y[which(colnames(sample.sheet) == 'Sample_Type')] == 'plasma'){
      SV.table[,paste0(y[which(colnames(sample.sheet) == 'Sample_Barcode')],'___duplex.called') := NA]
    }
    # edit DMP signout column
    if(y[which(colnames(sample.sheet) == 'Sample_Type')] == 'Tumor'){
      SV.table[!is.na(get(current.colname)),DMP := 'Signed out'] 
    }
  })
  # adding CH filler column
  SV.table$CH = ''
  
# append and write --------------------------------------------------------
  SNV.table <- fread(x)
  snv.sv.table.directory <- paste0(results.dir,'/results_',criteria,'_combined/',TDM1.ID,'_table.csv') 
  write.csv(rbind(SNV.table,SV.table),snv.sv.table.directory,quote = F,row.names = F)
  
})
