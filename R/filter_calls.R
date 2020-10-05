# library(data.table)
# library(tidyr)
# library(stringr)
# library(dplyr)

#' @export
filter_calls = function(
  master.ref,results.dir,
  CH.path = '/juno/work/access/production/resources/dmp_signedout_CH/current/signedout_CH.txt',
  criteria = 'stringent'
){
  # # test input section -----------------------------------------------------------
  # master.ref = fread('/juno/work/bergerm1/bergerlab/zhengy1/access_data_analysis/data/example_master_file.csv')
  # results.dir = paste0('/juno/work/bergerm1/MSK-ACCESS/ACCESS-Projects/test_access/access_data_analysis/output_042020/')
  # dmp.key.path = '/ifs/dmprequest/12-245/key.txt'
  # CH.calls = fread('/ifs/work/bergerm1/zhengy1/RET_all/Original_files/signedout_CH.txt')
  # # criteria <- 'permissive'
  # criteria <- 'stringent'
  #
  # criteria definition -----------------------------------------------------
  if(criteria == 'permissive'){
    hotspot.support <- 1
    non.hotspot.support <- 3
  }else{
    hotspot.support <- 3
    non.hotspot.support <- 5
  }

  dir.create(paste0(results.dir,'/results_',criteria))

  # inputs ---------------------------------------------------------------
  # DMP.key <- fread(dmp.key.path)
  CH.calls = fread(CH.path)
  pooled.normal.mafs <-
    fread(paste0(results.dir,'/pooled/all_all_unique.maf')) %>%
    mutate(Tumor_Sample_Barcode = paste0(Tumor_Sample_Barcode,'___pooled')) %>%
    select(Hugo_Symbol,Tumor_Sample_Barcode,Chromosome,Start_Position,End_Position,Variant_Classification,HGVSp_Short,Reference_Allele,Tumor_Seq_Allele2,t_alt_count) %>%
    group_by(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,Reference_Allele,Tumor_Seq_Allele2) %>%
    summarise(duplex_support_num = length(which(t_alt_count >= 2))) %>%
    filter(duplex_support_num > 0,.preserve = T) %>%
    data.table()

  # for each patient produce the correct results ----------------------------
  # x <- unique(master.ref$cmo_patient_id)[1]
  all.fillout.dim <- lapply(unique(master.ref$cmo_patient_id),function(x){
    print(paste0('Processing patient ',x))

    # Inputs and sanity checks ------------------------------------------------
    fillouts.filenames <- list.files(paste0(results.dir,'/',x,'/'),'ORG-STD_genotyped.maf|ORG-SIMPLEX-DUPLEX_genotyped.maf',full.names = T)

    # compiling a sample sheet with duplex, simplex, normal, DMP tumor and DMP normal
    sample.sheet <- fread(paste0(results.dir,'/',x,'/',x,'_sample_sheet.tsv'))[,.(Sample_Barcode,cmo_patient_id,Sample_Type)]
    simplex.sample.sheet = sample.sheet[Sample_Type == 'duplex',.(Sample_Barcode,cmo_patient_id,Sample_Type = 'simplex')]
    sample.sheet = rbind(sample.sheet,simplex.sample.sheet) %>%
      mutate(column.names = paste0(Sample_Barcode,'___',Sample_Type)) %>%
      data.table()

    # compiling different genotype files from step 1
    fillouts.dt <- do.call(rbind,lapply(fillouts.filenames,function(y){
      sample.name = gsub('.*./|-ORG.*.','',y)
      sample.type = sample.sheet[Sample_Barcode == sample.name]$Sample_Type

      # t_alt_count,t_ref_count,t_depth these columns are useless, have to use duplex/simplex/standard columms
      maf.file <- fread(y) %>%
        select(-c(t_alt_count,t_ref_count)) %>%
        data.table()

      if (nrow(maf.file) == 0) {

        columns <- c(
          "Hugo_Symbol", "Tumor_Sample_Barcode", "Chromosome", "Start_Position",
          "End_Position", "Variant_Classification", "HGVSp_Short",
          "Reference_Allele", "Tumor_Seq_Allele2", "t_var_freq", "ExAC_AF")
        df <- data.frame(matrix(ncol = length(columns), nrow = 0))
        colnames(df) <- columns

        return(df)
      }

      # fragment counts replacing actual allele counts
      if(grepl('SIMPLEX-DUPLEX_genotyped',y)){
        melt.id.vars = colnames(maf.file)[!grepl('fragment',colnames(maf.file))]

        # get rid of simplex duplex aggregate columns
        maf.file %>%
          select(-c(contains('simplex_duplex'))) %>%
          # melting and dcasting columns back but separating duplex and simplex columns
          # t_alt_duplex, t_depth_duplex, t_alt_simplex, t_depth_simplex --> t_alt, t_depth
          melt.data.table(id.vars = melt.id.vars,variable.name = 'variable',value.name = 'value') %>%
          mutate(variable = gsub('fragment','_',variable)) %>%
          separate(variable,c('variable','Sample_Type'),sep = '___') %>%
          mutate(Tumor_Sample_Barcode = paste0(sample.name,'___',Sample_Type)) %>%
          select(-Sample_Type) %>%
          data.table() %>%
          unique() %>%
          dcast.data.table(as.formula(paste0(paste0(melt.id.vars,collapse = ' + '),' ~ variable')),value.var = 'value') -> maf.file
      }else{
        maf.file <- maf.file %>%
          mutate(Tumor_Sample_Barcode = paste0(sample.name,'___',sample.type)) %>%
          # swaping the t_alt_count(etc)_standard for t_alt_count(etc)
          mutate(t_alt_count = t_alt_count_standard,t_total_count = t_total_count_standard)
      }

      maf.file = maf.file %>%
        mutate(t_var_freq = paste0(t_alt_count,'/',t_total_count,'(',round(t_alt_count/t_total_count,4),')')) %>%
        transmute(Hugo_Symbol,Tumor_Sample_Barcode,Chromosome = as.character(Chromosome),Start_Position,End_Position,Variant_Classification,
                  HGVSp_Short,Reference_Allele,Tumor_Seq_Allele2,t_var_freq,ExAC_AF) %>%
        data.table()

      return(maf.file)
    })) %>%
      unique() %>%
      data.table()

    # merging and melting -----------------------------------------------------
    hotspot.maf <- fread(paste0(results.dir,'/',x,'/',x,'_all_unique_calls_hotspots.maf')) %>%
      rowwise() %>%
      transmute(Hugo_Symbol,Chromosome = as.character(Chromosome),Start_Position,End_Position,Variant_Classification,
                # HGVSp_Short,
                Reference_Allele,Tumor_Seq_Allele2,Hotspot = ifelse(hotspot_whitelist,'Hotspot',NA)) %>%
      data.table()

    dmp.maf <- fread(paste0(results.dir,'/',x,'/',x,'_impact_calls.maf')) %>%
      transmute(Hugo_Symbol,Chromosome = as.character(Chromosome),Start_Position,End_Position,Variant_Classification,
                # HGVSp_Short,
                Reference_Allele,Tumor_Seq_Allele2) %>%
      mutate(DMP = 'Signed out') %>%
      unique() %>%
      data.table()

    if((nrow(dmp.maf) > 0) && (nrow(fillouts.dt) > 0)){
        fillouts.dt <- fillouts.dt %>%
          dcast.data.table(Hugo_Symbol + Chromosome + Start_Position + End_Position + Variant_Classification +
                           HGVSp_Short + Reference_Allele + Tumor_Seq_Allele2 + ExAC_AF ~ Tumor_Sample_Barcode,
                           value.var = 't_var_freq') %>%
          # hotspot information
          merge(
            hotspot.maf,
            by = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification','Reference_Allele','Tumor_Seq_Allele2'),
            all.x = T) %>%
          # Identifying signed out calls
          merge(
            dmp.maf,
            by = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification','Reference_Allele','Tumor_Seq_Allele2'),
            all.x = T) %>%
          # pooled normal for systemic artifacts
          merge(
            pooled.normal.mafs,
            by = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification','Reference_Allele','Tumor_Seq_Allele2'),
            all.x = T) %>%
          data.table()
    } else if (nrow(fillouts.dt) > 0){
      fillouts.dt <- fillouts.dt %>%
        dcast.data.table(
          Hugo_Symbol + Chromosome + Start_Position + End_Position + Variant_Classification +
          HGVSp_Short + Reference_Allele + Tumor_Seq_Allele2 + ExAC_AF ~ Tumor_Sample_Barcode,
          value.var = 't_var_freq') %>%
        # hotspot information
        merge(
          hotspot.maf,
          by = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification','Reference_Allele','Tumor_Seq_Allele2'),
          all.x = T) %>%
        # pooled normal for systemic artifacts
        merge(
          pooled.normal.mafs,
          by = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification','Reference_Allele','Tumor_Seq_Allele2'),
          all.x = T) %>%
        mutate(DMP = NA) %>%
        data.table()
    } else {
      # if fillouts.dt has no data, then add the needed columns with no data
      fillouts.dt[,c("DMP", "Hotspot", "duplex_support_num") := NA]
    }

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

    # check if the table is emtpy, if so then add the remaining columns and exit

    if (nrow(fillouts.dt) == 0) {
      fillouts.dt[,c("call_confidence", "CH") := NA]

      fillouts.dt <- fillouts.dt %>% select(
        Hugo_Symbol,Chromosome,Start_Position,End_Position,
        Variant_Classification,HGVSp_Short,Reference_Allele,Tumor_Seq_Allele2,
        ExAC_AF,Hotspot,DMP,CH,duplex_support_num,call_confidence,sort(everything()))

      write.csv(
        fillouts.dt,
        paste0(results.dir,'/results_',criteria,'/',x,'_SNV_table.csv'),
        row.names = F)
      return
    }

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
      # print(table(fillouts.dt[,get(paste0(tmp.col.name,'.called'))]))
    })

    if(all(!c('unfilterednormal','normal_DMP') %in% sample.sheet$Sample_Type)){
      tmp.col.name <- plasma.samples[1]
      lapply(plasma.samples,function(tmp.col.name){
        #fillouts.dt[as.numeric(gsub("\\(|\\)",'',str_extract(get(tmp.col.name),"\\(.*.\\)"))) >= 0.3 | ExAC_AF >= 0.0001,eval(paste0(tmp.col.name,'.called')) := 'Not Called']
        fillouts.dt[get(tmp.col.name)  == '0/0(NaN)',eval(paste0(tmp.col.name,'.called')) := 'Not Covered']
      })
    }else{
      lapply(plasma.samples,function(tmp.col.name){
        lapply(normal.samples,function(tmp.col.name.normal){
          # duplex tvar/nvar > 5
          fillouts.dt[(as.numeric(gsub("\\(|\\)",'',str_extract(get(tmp.col.name),"\\(.*.\\)")))/as.numeric(gsub("\\(|\\)",'',str_extract(get(tmp.col.name.normal),"\\(.*.\\)"))) < 2) |
                        # if duplex have no reads, use simplex tvar
                        (as.numeric(gsub("\\(|\\)",'',str_extract(get(gsub('duplex','simplex',tmp.col.name)),"\\(.*.\\)")))/as.numeric(gsub("\\(|\\)",'',str_extract(get(tmp.col.name.normal),"\\(.*.\\)"))) < 2 &
                           as.numeric(gsub("/.*.$",'',get(tmp.col.name)))  == 0),
                      eval(paste0(tmp.col.name,'.called')) := 'Not Called']
          fillouts.dt[get(tmp.col.name)  == '0/0(NaN)',eval(paste0(tmp.col.name,'.called')) := 'Not Covered']
        })
      })
    }

    # final processing --------------------------------------------------------
    # Save only the useful column
    #print(fillouts.dt)
    #print("#######")
    fillouts.dt <-  fillouts.dt[DMP == 'Signed out' | fillouts.dt[,apply(.SD,1,function(x){any(x == 'Called')})]]
    #print(fillouts.dt)
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

    write.csv(fillouts.dt,paste0(results.dir,'/results_',criteria,'/',x,'_SNV_table.csv'),row.names = F)
  })

  if(all(unlist(all.fillout.dim))){
    print('All dimension of  fillout mafs for each patient looks correct')
  }

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
  parser$add_argument('-ch', '--chlist', type='character', default = '/juno/work/access/production/resources/dmp_signedout_CH/current/signedout_CH.txt',
                      help='List of signed out CH calls [default]')
  parser$add_argument('-c', '--criteria', type='character', default = 'stringent',
                      help='Calling criteria [default]')
  args=parser$parse_args()

  master.ref = args$masterref
  results.dir = args$resultsdir
  chlist = args$chlist
  criteria = args$criteria

  cat(paste0(paste0(c(paste0(rep('-',15),collapse = ''),'Arguments input: ',master.ref,results.dir,chlist,criteria,
                      paste0(rep('-',15),collapse = '')),collapse = "\n"),'\n'))

  if(!criteria %in% c('stringent','permissive')){
    stop('Criteria argument should be either stringent or permissive')
  }

  suppressWarnings(filter_calls(fread(master.ref),results.dir,chlist,criteria))

}
