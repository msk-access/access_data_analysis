library(data.table)
library(tidyr)
library(stringr)
library(dplyr)
source('/ifs/work/bergerm1/zhengy1/RET_all/Code/genotype_source_code.R')

master.ref <- fread('/juno/work/bergerm1/bergerlab/zhengy1/TDM1/Sample_mapping/master_ref_022520.csv')
# vardict.dir <- '/juno/work/bergerm1/bergerlab/zhengy1/TDM1/snv_pipeline/maf_081919/'
results.dir <- '/juno/work/bergerm1/bergerlab/zhengy1/TDM1/Analysis_files/run_022520/'
CH.calls = fread('/ifs/work/bergerm1/zhengy1/RET_all/Original_files/signedout_CH.txt')

# criteria <- 'permissive'
criteria <- 'stringent'

if(criteria == 'permissive'){
  hotspot.support <- 1
  non.hotspot.support <- 3
}else{
  hotspot.support <- 3
  non.hotspot.support <- 5
}

dir.create(paste0(results.dir,'/results_',criteria))

# DMP stuff ---------------------------------------------------------------
DMP.key <- fread('/ifs/dmprequest/12-245/key.txt')

# Pooled normals ----------------------------------------------------------
duplexsupport <- function(x) {
  print(x)
  length(which(x >= 2))
}
pooled.normal.mafs <-
  fread(paste0(results.dir,'/pooled/all_all_unique_fillout_maf2maf.maf')) %>% mutate(Tumor_Sample_Barcode = paste0(Tumor_Sample_Barcode,'___pooled')) %>%
  select(Hugo_Symbol,Tumor_Sample_Barcode,Chromosome,Start_Position,End_Position,Variant_Classification,HGVSp_Short,Reference_Allele,Tumor_Seq_Allele2,t_alt_count_fragment) %>% 
  group_by(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,Reference_Allele,Tumor_Seq_Allele2) %>% 
  summarise(duplex_support_num = length(which(t_alt_count_fragment >= 2))) %>% filter(duplex_support_num >= 2) %>% data.table()

# for each patient produce the correct results ----------------------------
x <- 'C-2NMC49'
x <- unique(master.ref$`CMO patient ID`)[3]
x <- 'C-JKUXLJ'
all.fillout.dim <- lapply(unique(master.ref$`CMO patient ID`),function(x){
  print(x)
# Inputs and sanity checks ------------------------------------------------
  fillouts.filenames <- list.files(paste0(results.dir,'/',x,'/fillout'),full.names = T)  
  fillouts.dt <- do.call(rbind,lapply(fillouts.filenames,function(y){
    maf.file <- fread(y) 
    # fragment counts replacing actual allele counts
    if(!grepl('normal_DMP_fillout_maf2maf|Tumor_fillout_maf2maf',y)){
      maf.file <- maf.file[,c('t_alt_count','t_ref_count','t_depth') := list(t_alt_count_fragment,t_ref_count_fragment,t_total_count_fragment)][,!c('t_alt_count_fragment','t_ref_count_fragment','t_total_count_fragment'),with = F] %>%
        mutate(Tumor_Sample_Barcode = gsub('_cl.*.FX-|_S.._..._cl.*.FX-|_cl.*.FX|_S.._..._cl.*.FX','___',Tumor_Sample_Barcode)) %>% 
        mutate(Tumor_Sample_Barcode = ifelse(grepl('-N[0-9][0-9][0-9]-d',Tumor_Sample_Barcode),paste0(Tumor_Sample_Barcode,'unfilterednormal'),Tumor_Sample_Barcode)) %>% data.table()
    }else{
      maf.file <- maf.file %>% merge(DMP.key[,.(V1,V2)],by.x = 'Tumor_Sample_Barcode',by.y = 'V2',all.x = T) %>%
        mutate(Tumor_Sample_Barcode = gsub(' ','',ifelse(grepl('-T',V1),paste0(V1,'___Tumor'),paste(V1,'___normal_DMP')))) %>% 
        select(-c(V1,t_alt_count_fragment,t_ref_count_fragment,t_total_count_fragment)) %>% data.table()
    }
    return(maf.file)
  }))
# # Testing (commented out because testing finished)  -----------------------
#   if(length(unique(table(fillouts.dt$Tumor_Sample_Barcode))) == 1){
#     print(paste0('All mafs for patient ',x,' are the same!'))
#     return(TRUE)
#   }else{
#     print(paste0('All mafs for patient ',x,' are not the same!'))
#     return(FALSE)
#   }
# merging and melting -----------------------------------------------------
  sample.sheet <- fread(paste0(results.dir,'/',x,'/',x,'_sample_sheet.tsv')) %>% rowwise() %>%
    mutate(Sample_Type = ifelse(Sample_Type == 'plasma','duplex',
                                ifelse(Sample_Type == 'normal','unfilterednormal',
                                       ifelse(Sample_Type == 'plasma_simplex','simplex',Sample_Type)))) %>%
    mutate(column.names = paste0(Sample_Barcode,'___',Sample_Type)) %>% data.table()
  hotspot.maf <- fread(paste0(results.dir,'/',x,'/',x,'_all_unique_calls_hotspots.maf')) %>% rowwise() %>%
    transmute(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,
              # HGVSp_Short,
              Reference_Allele,Tumor_Seq_Allele2,Hotspot = ifelse(hotspot_whitelist,'Hotspot',NA)) %>% data.table()
  dmp.maf <- fread(paste0(results.dir,'/',x,'/',x,'_impact_calls.maf')) %>%
    select(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,
           # HGVSp_Short,
           Reference_Allele,Tumor_Seq_Allele2) %>% mutate(DMP = 'Signed out') %>% unique() %>% data.table()
  
  fillouts.dt <-
    fillouts.dt %>% mutate(t_var_freq = paste0(t_alt_count,'/',t_depth,'(',round(t_alt_count/t_depth,4),')')) %>% 
      select(Hugo_Symbol,Tumor_Sample_Barcode,Chromosome,Start_Position,End_Position,Variant_Classification,HGVSp_Short,Reference_Allele,
             Tumor_Seq_Allele2,t_var_freq,ExAC_AF) %>% spread(Tumor_Sample_Barcode,t_var_freq) %>%
      # hotspot information 
      merge(hotspot.maf,by = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification',
                               # 'HGVSp_Short',
                               'Reference_Allele','Tumor_Seq_Allele2'),all.x = T) %>%
      # Identifying signed out calls
      merge(dmp.maf,by = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification','Reference_Allele','Tumor_Seq_Allele2'),all.x = T) %>%
      # pooled normal for systemic artifacts
      merge(pooled.normal.mafs,by = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification','Reference_Allele','Tumor_Seq_Allele2'),all.x = T) %>% 
      filter(is.na(duplex_support_num) | !is.na(DMP)) %>% data.table()
  # Interesting cases where DMP signed out calls are artifacets
  if(any(!is.na(fillouts.dt$DMP) & !is.na(fillouts.dt$duplex_support_num))){
    print(paste0('Look at ',x,' for DMP signed out plasma artifacts...'))
  }
# germline filtering for matched and unmatched ----------------------------
  plasma.samples <- sample.sheet[Sample_Type %in% c('duplex')]$column.names
  normal.samples <- sample.sheet[Sample_Type %in% c('unfilterednormal','normal_DMP')]$column.names
  fillouts.dt[,c(
    paste0(plasma.samples,'.called')
    # paste0(gsub('duplex','simplex',plasma.samples),'.called')
    
  ) := 'Not Called']
  # preliminary calling
  # tmp.col.name <- plasma.samples[1]
  lapply(plasma.samples,function(tmp.col.name){
    # genotyping (signed out stuff)
    fillouts.dt[(as.numeric(gsub("/.*.$",'',get(tmp.col.name))) >= 1 | as.numeric(gsub("/.*.$",'',get(paste0(gsub('duplex','simplex',tmp.col.name))))) > 1) & DMP == 'Signed out',
                eval(paste0(tmp.col.name,'.called')) := 'Genotyped']
                # c(eval(paste0(tmp.col.name,'.called')),eval(paste0(gsub('duplex','simplex',tmp.col.name),'.called'))) := list('Called','Called')]
    # hotspot reads
    fillouts.dt[as.numeric(gsub("/.*.$",'',get(tmp.col.name))) >= hotspot.support & Hotspot == 'Hotspot',
                eval(paste0(tmp.col.name,'.called')) := 'Called']
                # c(eval(paste0(tmp.col.name,'.called')),eval(paste0(gsub('duplex','simplex',tmp.col.name),'.called'))) := list('Called','Called')]
    # non hotspot reads
    fillouts.dt[as.numeric(gsub("/.*.$",'',get(tmp.col.name))) >= non.hotspot.support & is.na(Hotspot),
                eval(paste0(tmp.col.name,'.called')) := 'Called']
                # c(eval(paste0(tmp.col.name,'.called')),eval(paste0(gsub('duplex','simplex',tmp.col.name),'.called'))) := list('Called','Called')]
    print(table(fillouts.dt[,get(paste0(tmp.col.name,'.called'))]))
  })
  if(all(!c('unfilterednormal','normal_DMP') %in% sample.sheet$Sample_Type)){
    tmp.col.name <- plasma.samples[1]
    lapply(plasma.samples,function(tmp.col.name){
      fillouts.dt[as.numeric(gsub("\\(|\\)",'',str_extract(get(tmp.col.name),"\\(.*.\\)"))) >= 0.3 | ExAC_AF >= 0.0001,eval(paste0(tmp.col.name,'.called')) := 'Not Called']
      fillouts.dt[get(tmp.col.name)  == '0/0(NaN)',eval(paste0(tmp.col.name,'.called')) := 'Not Covered']
    })
  }else{
    lapply(plasma.samples,function(tmp.col.name){
      lapply(normal.samples,function(tmp.col.name.normal){
                    # duplex tvar/nvar > 5
        fillouts.dt[(as.numeric(gsub("\\(|\\)",'',str_extract(get(tmp.col.name),"\\(.*.\\)")))/as.numeric(gsub("\\(|\\)",'',str_extract(get(tmp.col.name.normal),"\\(.*.\\)"))) < 5) |
                      # if duplex have no reads, use simplex tvar
                      (as.numeric(gsub("\\(|\\)",'',str_extract(get(gsub('duplex','simplex',tmp.col.name)),"\\(.*.\\)")))/as.numeric(gsub("\\(|\\)",'',str_extract(get(tmp.col.name.normal),"\\(.*.\\)"))) < 5 &
                         as.numeric(gsub("/.*.$",'',get(tmp.col.name)))  == 0),
                    eval(paste0(tmp.col.name,'.called')) := 'Not Called']
        fillouts.dt[get(tmp.col.name)  == '0/0(NaN)',eval(paste0(tmp.col.name,'.called')) := 'Not Covered']
      })
    })
  }

# final processing --------------------------------------------------------
  # Save only the useful column
  fillouts.dt <-  fillouts.dt[DMP == 'Signed out' | fillouts.dt[,apply(.SD,1,function(x){any(x == 'Called')})]]  
  # combining duplex and simplex counts
  lapply(plasma.samples,function(tmp.col.name){
    # hotspot reads
    fillouts.dt[,eval(gsub('duplex','total',tmp.col.name)) := paste0(
      as.numeric(gsub("/.*.$",'',get(tmp.col.name)))+as.numeric(gsub("/.*.$",'',get(gsub('duplex','simplex',tmp.col.name)))),'/',
      as.numeric(gsub("^.*./|\\(.*.$",'',get(tmp.col.name)))+as.numeric(gsub("^.*./|\\(.*.$",'',get(gsub('duplex','simplex',tmp.col.name)))),'(',
      round((as.numeric(gsub("/.*.$",'',get(tmp.col.name)))+as.numeric(gsub("/.*.$",'',get(gsub('duplex','simplex',tmp.col.name)))))/
              (as.numeric(gsub("^.*./|\\(.*.$",'',get(tmp.col.name)))+as.numeric(gsub("^.*./|\\(.*.$",'',get(gsub('duplex','simplex',tmp.col.name))))),4),')'
    )]
    fillouts.dt[,c(eval(gsub('duplex','simplex',tmp.col.name)),eval(tmp.col.name)):= list(NULL,NULL)]
  })
  fillouts.dt <- fillouts.dt[,order(colnames(fillouts.dt)),with = F] %>% 
    # filter for artifacts
    mutate(call_confidence = case_when(
      (Hugo_Symbol == 'TERT' & is.na(Hotspot)) | (Hugo_Symbol == 'ERBB2' & grepl('[A-Z]90[0-9][A-Z]',HGVSp_Short)) | 
        (Hugo_Symbol == 'BRAF' & grepl('711',HGVSp_Short)) | (Hugo_Symbol == 'NF1' & grepl('[A-Z]106[0-9][A-Z]',HGVSp_Short)) ~ 'Low',
      DMP == 'Signed out' ~ 'High',
      TRUE ~ ''
    )) %>%
    merge(CH.calls[,.(Hugo_Symbol = Gene,Chromosome = Chrom,Start_Position = Start,Reference_Allele = Ref,Tumor_Seq_Allele2 = Alt,HGVSp_Short = AAchange,Variant_Classification = VariantClass,CH = 'Yes')],
          by = c('Hugo_Symbol','Chromosome','Start_Position','Variant_Classification','HGVSp_Short','Reference_Allele','Tumor_Seq_Allele2'),
          all.x = T) %>%
    mutate(CH = ifelse(is.na(CH),'No','Yes')) %>%
    select(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,HGVSp_Short,Reference_Allele,Tumor_Seq_Allele2,
           ExAC_AF,Hotspot,DMP,CH,duplex_support_num,call_confidence,sort(everything())) 
  
  write.csv(fillouts.dt,paste0(results.dir,'/results_',criteria,'/',
                               # gsub('-','_',master.ref[`CMO patient ID` == x]$GRAIL_PT_ID[1]),
                               master.ref[`CMO patient ID` == x]$GRAIL_PT_ID[1],'_SNV_table.csv'),row.names = F)
})





if(all(unlist(all.fillout.dim))){
  print('All dimension of  fillout mafs for each patient looks correct')
}
