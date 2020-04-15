library(data.table)

repo.dir = '/juno/work/bergerm1/bergerlab/zhengy1/access_data_analysis/'
test.data.dir = '/juno/work/bergerm1/MSK-ACCESS/ACCESS-Projects/test_access/access_data_analysis/' 
test.bam.dir = paste0(test.data.dir,'/bams')
test.maf.dir = paste0(test.data.dir,'/mafs')
test.cna.dir = paste0(test.data.dir,'/cnas')
test.sv.dir = paste0(test.data.dir,'/svs')
master.ref = fread(paste0(repo.dir,'/data/example_master_file_raw.csv')) %>% rowwise %>%
  mutate(BAM_path.plasma.duplex = list.files(test.bam.dir,paste0(`CMO sample ID.plasma`,'.*.duplex.bam'),full.names = T),
         BAM_path.plasma.simplex = list.files(test.bam.dir,paste0(`CMO sample ID.plasma`,'.*.simplex.bam'),full.names = T),
         BAM_path.normal = list.files(test.bam.dir,paste0(`CMO sample ID.normal`,'.*.FX.bam'),full.names = T),
         maf_path = list.files(test.maf.dir,paste0(`CMO sample ID.plasma`,'.*.maf'),full.names = T),
         cna_path = list.files(test.cna.dir,paste0(`CMO sample ID.plasma`,'.*.txt'),full.names = T),
         sv_path = list.files(test.sv.dir,paste0(`CMO sample ID.plasma`,'.*.txt'),full.names = T)) %>% select(-BAM_path.plasma) %>%
  setnames(c('CMO patient ID','CMO sample ID.plasma','BAM_path.plasma.duplex','BAM_path.plasma.simplex','CMO sample ID.normal',
             'BAM_path.normal','Paired','Sex','collection_date','DMP_ID','maf_path','cna_path','sv_path'),
           c('cmo_patient_id','cmo_sample_id_plasma','bam_path_plasma_duplex','bam_path_plasma_simplex','cmo_sample_id_normal',
             'bam_path_normal','paired','sex','collection_date','dmp_sample_id','maf_path','cna_path','sv_path'))
write.csv(master.ref,paste0(repo.dir,'/data/example_master_file.csv'),quote = F,row.names = F)