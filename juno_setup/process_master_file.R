library(data.table)

repo.dir = '/juno/work/bergerm1/bergerlab/zhengy1/access_data_analysis/'
test.bam.dir = '/juno/work/bergerm1/MSK-ACCESS/ACCESS-Projects/test_access/access_data_analysis/bams'
master.ref = fread(paste0(repo.dir,'/data/example_master_file_raw.csv')) %>% rowwise %>%
  mutate(BAM_path.plasma = list.files(test.bam.dir,paste0(`CMO sample ID.plasma`,'.*.duplex.bam'),full.names = T),
         BAM_path.normal = list.files(test.bam.dir,paste0(`CMO sample ID.normal`,'.*.duplex.bam'),full.names = T))
write.csv(master.ref,paste0(repo.dir,'/data/example_master_file.csv'),quote = F,row.names = F)