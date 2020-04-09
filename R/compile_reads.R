library(data.table)
library(tidyr)
library(stringr)
library(dplyr)
source('/ifs/work/bergerm1/zhengy1/RET_all/Code/genotype_source_code.R')

master.ref <- fread('/juno/work/bergerm1/bergerlab/zhengy1/TDM1/Sample_mapping/master_ref_022520.csv')
# vardict.dir <- '/ifs/work/bergerm1/zhengy1/RET_all/vardict_output/vardict_run_053119/'
vardict.dir <- '/juno/work/bergerm1/bergerlab/zhengy1/TDM1/snv_pipeline/maf_101019/'
results.dir <- paste0('/juno/work/bergerm1/bergerlab/zhengy1/TDM1/Analysis_files/run_',format(Sys.time(),'%m%d%y'))
# results.dir <- paste0('/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/run_temp')
dir.create(results.dir)
# make tmp directory in output directory
dir.create(paste0(results.dir,'/tmp'))

# data from DMP -----------------------------------------------------------
DMP.maf <- fread('/ifs/work/bergerm1/zhengy1/dmp/mskimpact/data_mutations_extended.txt') %>% filter(Mutation_Status != 'GERMLINE') %>% data.table()
dim(DMP.maf)
DMP.key <- fread('/ifs/dmprequest/12-245/key.txt')
dim(DMP.key)

DMP.RET.maf <- DMP.maf[which(grepl(paste0(unique(master.ref[!is.na(DMP_ID)]$DMP_ID),collapse = '|'),DMP.maf$Tumor_Sample_Barcode)),]

# data from plasma --------------------------------------------------------
# total.maf <- do.call(rbind,lapply(list.files(paste0(vardict.dir,'frag/'),pattern = 'anno.maf',full.names = T),function(x){
#                                                fread(x) %>% filter(as.numeric(t_alt_count_fragment) > 0) %>% data.table()
#                                              })) %>% mutate(Tumor_Sample_Barcode = gsub('_S.._001.*|_cl.*','',Tumor_Sample_Barcode)) %>% data.table()
total.maf <- do.call(rbind,lapply(list.files(vardict.dir,pattern = 'combined-variants.vep_keptrmv_taggedHotspots.maf',recursive = T,full.names = T),function(x){
  fread(x) %>% filter(as.numeric(t_alt_count) > 0) %>% data.table()
})) %>% mutate(Tumor_Sample_Barcode = gsub('_S.._001.*|_cl.*','',Tumor_Sample_Barcode)) %>% data.table()
# unique(total.maf$Tumor_Sample_Barcode)

# Pooled normal samples ---------------------------------------------------
pooled.bams <- list.files('/ifs/work/bergerm1/ACCESS-Projects/novaseq_curated_duplex_v2/',pattern = '.bam',full.names = T)

# For each patient --------------------------------------------------------
# x <- unique(master.ref$`CMO patient ID`)[4]
# x <- unique(master.ref$`CMO patient ID`)[16]
# x = 'C-YW82CY'
all.fillout.id <- lapply(unique(master.ref$`CMO patient ID`),function(x){
  print(x)
  dir.create(paste0(results.dir,'/',x))
  # sample sheet with colummns -- TSB, sample type, bam path, treatm --------
  # need to get DMP tumor, DMP normal, plasma, plasma normal (if there is any), pooled normal
  # DMP sample sheet
  if(
    !any(is.na(master.ref[`CMO patient ID` == x]$DMP_ID)) & 
    x != 'C-YW82CY'
  ){
    dmp.sample.sheet <- data.frame(Sample_Barcode = DMP.key[grepl(unique(master.ref[`CMO patient ID` == x]$DMP_ID),V1)]$V1,
                                   BAM_path = paste0('/ifs/dmpshare/share/irb12_245/',DMP.key[grepl(unique(master.ref[`CMO patient ID` == x]$DMP_ID),V1)]$V2,'.bam')) %>%
      mutate(`CMO patient ID` = x,Sample_Type = ifelse(grepl('-T',Sample_Barcode),'Tumor','normal_DMP'),DMP_ID = unique(master.ref[`CMO patient ID` == x]$DMP_ID))
  }else{dmp.sample.sheet = NULL}
  # pooled.normal.sample.sheet <-  data.frame(Sample_Barcode = gsub('^.*.duplex_bams//|.bam','',pooled.bams),BAM_path = pooled.bams) %>%
  #   mutate(`CMO patient ID` = x,Sample_Type = 'normal_pooled',DMP_ID = unique(master.ref[`CMO patient ID` == x]$DMP_ID))
  # total sample sheet
  sample.sheet <- master.ref[`CMO patient ID` == x,.(Tumor_Sample_Barcode.plasma,Tumor_Sample_Barcode.normal,BAM_path.plasma,BAM_path.normal,`CMO patient ID`,DMP_ID)] %>%
    # CMO sample sheet
    unite(plasma, matches('plasma'),sep = '___') %>% unite(normal,matches('normal'),sep = '___') %>% gather(Sample_Type,Sample_Info,plasma,normal,na.rm = T) %>%
    separate(Sample_Info,c('Sample_Barcode','BAM_path'),sep = '___') %>% filter(!is.na(Sample_Barcode) & Sample_Barcode != 'NA') %>%
    rbind(dmp.sample.sheet) %>%
    # rbind(pooled.normal.sample.sheet) %>%
    data.table()
  # adding simplex
  sample.sheet <- rbind(sample.sheet,data.table(t(apply(sample.sheet[Sample_Type == 'plasma'],1,
                                                        function(x){gsub('duplex','simplex',x)}))) %>%
                          mutate(Sample_Type = 'plasma_simplex')) %>% unique()
  write.table(sample.sheet,paste0(results.dir,'/',x,'/',x,'_sample_sheet.tsv'),sep = '\t',quote = F,row.names = F)
  # piece together all unique calls -----------------------------------------
  # get duplex calls
  duplex.calls <- total.maf[Tumor_Sample_Barcode %in% sample.sheet$Sample_Barcode]
  # get impact calls
  impact.calls <- DMP.RET.maf[Tumor_Sample_Barcode %in% sample.sheet$Sample_Barcode]
  write.table(impact.calls[,.(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,HGVSp_Short,Reference_Allele,Tumor_Seq_Allele2)],
              paste0(results.dir,'/',x,'/',x,'_impact_calls.maf'),sep = '\t',quote = F,row.names = F)
  all.calls <- rbind(duplex.calls[,intersect(colnames(total.maf),colnames(DMP.RET.maf)),with = F],impact.calls[,intersect(colnames(total.maf),colnames(DMP.RET.maf)),with = F])
  all.calls <- all.calls[which(!duplicated(all.calls[,.(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,HGVSp_Short,Reference_Allele,Tumor_Seq_Allele2)])),] %>%
    filter(Variant_Classification != 'Silent' & !grepl('RP11-',Hugo_Symbol) & !grepl('Intron',Variant_Classification))
  write.table(all.calls,paste0(results.dir,'/',x,'/',x,'_all_unique_calls.maf'),sep = '\t',quote = F,row.names = F)
  # tagging hotspots
  system(paste0(
    'bsub  -R "rusage[mem=4]" -cwd ',results.dir,'/',x,'/ -oo %J_hotspot.o -eo %J_hotspot.e -W 00:59 ',
    ' python /ifs/work/bergerm1/zhengy1/ACCESS-Pipeline-DEV/cwl_tools/hotspots/tag_hotspots.py ',
    ' -m ',results.dir,'/',x,'/',x,'_all_unique_calls.maf',
    ' -itxt /ifs/work/bergerm1/Innovation/Resources/Hotspots/hotspot-list-union-v1-v2_with_TERT.txt ',
    ' -o ',results.dir,'/',x,'/',x,'_all_unique_calls_hotspots.maf',
    ' -outdir ',results.dir,'/',x,'/',x
  ))

  # genotype all bams in this patient directory -----------------------------
  # genotyping plasma samples -- plasma duplex&simplex, plasma normal, pooled plasma normal
  system(paste0(
    'python /ifs/work/bergerm1/Innovation/sandbox/Maysun/projects/PanCancer/PR_BRCA/tier1+2/fillout/patient_mafs/ConvertMafForFillout.py ',
    results.dir,'/',x,'/',x,'_all_unique_calls.maf'
  ))
  dir.create(paste0(results.dir,'/',x,'/fillout'))
  dir.create(paste0(results.dir,'/',x,'/portal'))
  job.ids <- apply(sample.sheet[Sample_Type %in% c('plasma','normal','normal_pooled','plasma_simplex')],1,function(y){
    genotype_plasma_bams(bam_path = y[which(colnames(sample.sheet) == 'BAM_path')],
                         sample_barcode = y[which(colnames(sample.sheet) == 'Sample_Barcode')],
                         CMO_ID = y[which(colnames(sample.sheet) == 'CMO patient ID')], maf_dir = results.dir,sample_type = y[which(colnames(sample.sheet) == 'Sample_Type')],
                         tmp_dir = paste0(results.dir,'/tmp'))
  })
  # genotyping impact samples -- DMP tumor, DMP normal
  if(nrow(sample.sheet[!Sample_Type %in% c('plasma','normal','normal_pooled','plasma_simplex')]) > 0){
    job.ids <- c(job.ids,apply(sample.sheet[!Sample_Type %in% c('plasma','normal','normal_pooled','plasma_simplex')],1,function(y){
      genotype_impact_bams(bam_path = y[which(colnames(sample.sheet) == 'BAM_path')],
                           sample_barcode = y[which(colnames(sample.sheet) == 'Sample_Barcode')],
                           CMO_ID = y[which(colnames(sample.sheet) == 'CMO patient ID')], maf_dir = results.dir,sample_type = y[which(colnames(sample.sheet) == 'Sample_Type')],
                           tmp_dir = paste0(results.dir,'/tmp'))
    }))
  }
})


# Get base count multi sample in pooled normal ----------------------------
# all all unique calls in entire cohort
dir.create(paste0(results.dir,'/pooled'))
all.all.unique.temp.mafs <- do.call(rbind,lapply(unique(master.ref$`CMO patient ID`),function(x){
  fread(list.files(paste0(results.dir,'/',x),pattern = 'temp.maf',full.names = T))
}))
all.all.unique.temp.mafs <- all.all.unique.temp.mafs[-which(duplicated(all.all.unique.temp.mafs[,.(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,HGVSp_Short,Reference_Allele,Tumor_Seq_Allele2)]))]
write.table(all.all.unique.temp.mafs,paste0(results.dir,'/pooled/all_all_unique-temp.maf'),sep = '\t',quote = F,row.names = F)

pooled.sample.job.id <- system(paste0(
  'bsub -W 12:00  -R "rusage[mem=30]" ',
  '-w ',' \"',paste0(paste0('done(',unlist(all.fillout.id),')'),collapse = '&&'),'\" ',
  ' -oo ',results.dir,'/pooled/GetBaseCountsMultiSample.o ',' -eo ',results.dir,'/pooled/GetBaseCountsMultiSample.e ',
  ' /home/hasanm/Innovation/software/maysun/GetBaseCountsMultiSample/GetBaseCountsMultiSample ',
  '--omaf --filter_duplicate 0 --thread 10 --maq 20 --fasta /ifs/depot/assemblies/H.sapiens/b37/b37.fasta ',
  '--maf ',results.dir,'/pooled/all_all_unique-temp.maf --fragment_count 1 --output ',
  results.dir,'/pooled/all_all_unique_fillout.maf ',
  paste0(paste0(' --bam ',gsub('^.*./|.bam','',pooled.bams),':',pooled.bams),collapse = '')
),intern = T)
pooled.sample.job.id <- as.numeric(gsub('Job <|> is.*.$','',pooled.sample.job.id))
print(pooled.sample.job.id)
system(paste0(
  # 'bsub -W 02:00 -R "rusage[mem=8]" -cwd ',results.dir,'/pooled',
  # ' -oo %J_fillout.o -eo %J_fillout.e -w \"done(',pooled.sample.job.id,')\"',' cmo_maf2maf --input-maf ',
  # results.dir,'/pooled/all_all_unique_fillout.maf ',
  # ' --output-maf ',results.dir,'/pooled/all_all_unique_fillout_maf2maf.maf ',
  # ' --retain-cols Center,Verification_Status,Validation_Status,Mutation_Status,',
  # 'Sequencing_Phase,Sequence_Source,Validation_Method,Score,BAM_file,Sequencer,',
  # 'Tumor_Sample_UUID,Matched_Norm_Sample_UUID,Caller,t_total_count_fragment,',
  # 't_ref_count_fragment,t_alt_count_fragment && rm -rf /ifs/work/bergerm1/zhengy1/RET_all/Code/*.maf'
  # 
  'bsub -W 12:00  -R "rusage[mem=8]" -cwd ',results.dir,'/pooled',
  ' -oo %J_fillout.o -eo %J_fillout.e ',
  ' -w \"done(',pooled.sample.job.id,')\"',
  ' /opt/common/CentOS_7-dev/bin/perl /opt/common/CentOS_6-dev/vcf2maf/v1.6.17/maf2maf.pl ',
  ' --input-maf ',results.dir,'/pooled/all_all_unique_fillout.maf ',
  ' --ncbi-build GRCh37 --tum-vad-col t_alt_count --custom-enst /opt/common/CentOS_6-dev/vcf2maf/v1.6.17/data/isoform_overrides_at_mskcc ',
  ' --tmp-dir ',tempfile(pattern = "file", tmpdir = paste0(results.dir,'/tmp'), fileext = ""),
  ' --vep-path /opt/common/CentOS_6-dev/vep/v95 --vep-forks 4 --nrm-rad-col n_ref_count --buffer-size 5000 ',
  ' --filter-vcf /opt/common/CentOS_6-dev/vep/cache/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz --nrm-vad-col n_alt_count --max-filter-ac 10 ',
  ' --vep-data /opt/common/CentOS_6-dev/vep/cache/ --tum-rad-col t_ref_count --nrm-depth-col n_depth ',
  ' --ref-fasta /ifs/depot/pi/resources/genomes/GRCh37/fasta/b37.fasta ',
  ' --output-maf ',results.dir,'/pooled/all_all_unique_fillout_maf2maf.maf ',
  ' --retain-cols Center,Verification_Status,Validation_Status,Mutation_Status,Sequencing_Phase,Sequence_Source,Validation_Method,Score,BAM_file,Sequencer,',
  'Tumor_Sample_UUID,Matched_Norm_Sample_UUID,Caller,t_total_count_fragment,t_ref_count_fragment,t_alt_count_fragment --species homo_sapiens --tum-depth-col t_depth ',
  '&& rm -rf /ifs/work/bergerm1/zhengy1/RET_all/Code/*.maf'
))



# # checking if all jobs ran correctly --------------------------------------
# out.log.filenames <- list.files(results.dir,pattern = '\\.o',recursive = T,full.names = T)
# out.log <- unlist(lapply(out.log.filenames,function(x){
#   file.content <- readLines(x)
#   job.status <- gsub('^.*.cluster <solar. ','',file.content[which(grepl('Subject: ',file.content))])
# }))
# if(all(out.log == 'Done')){
#   print('Everything ran correctly')
# }
