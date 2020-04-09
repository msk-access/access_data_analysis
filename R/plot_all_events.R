library(data.table)
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
source('/ifs/work/bergerm1/zhengy1/RET_all/Code/table_to_maf.R')
# convert naming to timepoint, get rid of uncovered impact and access calls
source('/ifs/work/bergerm1/zhengy1/RET_all/Code/process_maf_for_graphs.R')
# for plotting consistency
status_id = c('Called' = 19, 'Not Called' =  4, 'Signed out' = 15,
              'Not Signed out' = 13, 'Not Covered' = 8, 'Genotyped' = 17)

master.ref <- fread('/juno/work/bergerm1/bergerlab/zhengy1/TDM1/Sample_mapping/master_ref_031920.csv')
run.dir <- paste0('/juno/work/bergerm1/bergerlab/zhengy1/TDM1/Analysis_files/run_022520/')
cna.dir <- '/juno/work/bergerm1/bergerlab/zhengy1/TDM1/Analysis_files/cna_032120/'
criteria <- 'stringent'
# criteria <- 'permissive'
results.dir = paste0(run.dir,'/results_',criteria,'_combined_reviewed')
output.dir = paste0('/juno/work/bergerm1/bergerlab/zhengy1/TDM1/For_Bob/plot_',format(Sys.time(),'%m%d%y'))
dir.create(output.dir)
# process cna calls
source('/juno/work/bergerm1/bergerlab/zhengy1/TDM1/Code/CNA_processing.R')
# returning `total.cna.call.set` variable


x = list.files(results.dir,full.names = T)[2]
# x = '/juno/work/bergerm1/bergerlab/zhengy1/TDM1/Analysis_files/run_022520//results_stringent_combined_reviewed/TDM_PT49_table.csv'
x = '/juno/work/bergerm1/bergerlab/zhengy1/TDM1/Analysis_files/run_022520//results_stringent_combined_reviewed/MSK-ML-0039_table.csv'
lapply(list.files(results.dir,full.names = T),function(x){
  TDM1.ID <- gsub(paste0('.*./results_',criteria,'_combined_reviewed/+','|_table.csv'),'',x)
  CMO.ID <- unique(master.ref[GRAIL_PT_ID == TDM1.ID]$`CMO patient ID`) 
  print(TDM1.ID)
  print(CMO.ID)
  # THIS PLOTS PLASMA SAMPLES ONLY
  
  # SNV
  tmp.table = fread(x)[call_confidence == 'High' | grepl('Protein Fusion: in frame',HGVSp_Short)]
  tmp.sample.sheets = fread(list.files(paste0(run.dir,'/',CMO.ID),pattern = 'sample_sheet',full.names = T)) %>% rowwise() %>%
    mutate(Sample_Type = ifelse(Sample_Type == 'plasma','duplex',
                                ifelse(Sample_Type == 'normal','duplexnormal',
                                       ifelse(Sample_Type == 'plasma_simplex','simplex',Sample_Type)))) %>%
    mutate(column.names = paste0(Sample_Barcode,'___',Sample_Type)) %>% data.table()
  tmp.table = table_to_maf(tmp.table,tmp.sample.sheets)
  tmp.table = data.table(process_maf_for_graph(tmp.table))
  # preserving the original process_maf_for_graphs so there are not two versions going around
  # converting TSB back to original form and merge it back to the collection timeline info from Bob
  tmp.table$Tumor_Sample_Barcode = paste0(CMO.ID,'-',tmp.table$Tumor_Sample_Barcode,'-d')
  
  # CNA
  tmp.cna = total.cna.call.set[Tumor_Sample_Barcode %in% tmp.sample.sheets$Sample_Barcode] %>%
    filter(grepl('L...-d',Tumor_Sample_Barcode)) %>% data.table()
      
  # transform sample IDs into times
  transform.vector = structure(as.character(master.ref[GRAIL_PT_ID == TDM1.ID]$collection_date),
                               names = master.ref[GRAIL_PT_ID == TDM1.ID]$`CMO sample ID.plasma`)
  print(transform.vector)
  tmp.table$Tumor_Sample_Barcode = transform.vector[tmp.table$Tumor_Sample_Barcode]
  
  if(nrow(tmp.table) == 0 | all(tmp.table$t_alt_count == 0)){
    print('skiping to the next')
    if(nrow(tmp.cna)) stop(paste0('Need to make CNA only file for: ',x))
    return()
  }
  
  colourCount = nrow(unique(tmp.table[,.(Hugo_Symbol,HGVSp_Short)]))
  getPalette = colorRampPalette(brewer.pal(8, "Set2"))
  SNV.SV.plot = ggplot(tmp.table) + 
    geom_line(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                  color = paste0(Hugo_Symbol,' ',ifelse(grepl('^p\\.',HGVSp_Short),HGVSp_Short,'')),group = paste0(Hugo_Symbol,'_',HGVSp_Short))) + 
    geom_point(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                   color = paste0(Hugo_Symbol,' ',ifelse(grepl('^p\\.',HGVSp_Short),HGVSp_Short,'')),shape = call_confidence),size = 1.5) + 
    labs(title=TDM1.ID,x='Time Point', y='VAF') + #scale_x_date(date_labels = "%Y %b %d") +
    scale_shape_manual(values=status_id,name = 'Call Status') + scale_color_manual(values = getPalette(colourCount),name = 'Alteration') + 
    theme_minimal() + scale_y_log10() +
    theme(panel.grid.major = element_blank(),legend.position="top",legend.box = "vertical",
          axis.text.x = element_text(angle=45, hjust=1, face = 'bold'))
  print(SNV.SV.plot)
  
  if(nrow(tmp.cna) > 0){
    tmp.cna = tmp.cna %>% mutate(Tumor_Sample_Barcode = factor(Tumor_Sample_Barcode,unique(tmp.sample.sheets[Sample_Type == 'duplex']$Sample_Barcode))) %>%
      # expand table on all empty samples without any calls
      data.table() %>% dcast.data.table(Hugo_Symbol + CNA + Chromosome + Cyt ~ Tumor_Sample_Barcode,drop = c(TRUE, FALSE),fill = 0,value.var = 'fc') %>% 
      melt.data.table(id.vars = c('Hugo_Symbol','CNA','Chromosome','Cyt'),variable.name = 'Tumor_Sample_Barcode',value.name = 'fc') %>% data.table()
    tmp.cna$Tumor_Sample_Barcode = transform.vector[tmp.cna$Tumor_Sample_Barcode]

    colourCount = nrow(unique(tmp.cna[,.(Hugo_Symbol,CNA)]))
    getPalette = colorRampPalette(brewer.pal(8, "Set2"))
    CNA.plot = ggplot(tmp.cna) + 
      geom_bar(aes(x = Tumor_Sample_Barcode,y = abs(fc),fill = paste0(Hugo_Symbol,'_',CNA)),position="dodge", stat="identity") +
      labs(x='Time Point', y='Absolute fc') + #scale_x_date(date_labels = "%Y %b %d") +
      scale_fill_manual(values = getPalette(colourCount),name = 'Alteration') +
      theme_minimal() + theme(panel.grid.major = element_blank(),legend.position="bottom",axis.text.x = element_text(angle=45, hjust=1,face = 'bold'))
    print(CNA.plot)
    
    pdf(paste0(output.dir,'/',TDM1.ID,'_all_events.pdf'),width = 10,height = 7)
    print(ggarrange(SNV.SV.plot,CNA.plot,ncol = 1,heights = c(2,1)))
    dev.off()
    
  }else{
    pdf(paste0(output.dir,'/',TDM1.ID,'_all_events.pdf'),width = 10,height = 7)
    print(SNV.SV.plot)
    dev.off()
    
  }
  
})
  