library(data.table)
library(tidyr)
library(stringr)
library(dplyr)


#' @export
compile_reads <- function(
                          master.ref, results.dir, project.ID, pooled.bam.dir = "/juno/work/access/production/resources/msk-access/current/novaseq_curated_duplex_bams_dmp/current/",
                          fasta.path = "/juno/work/access/production/resources/reference/current/Homo_sapiens_assembly19.fasta",
                          genotyper.path = "/ifs/work/bergerm1/Innovation/software/maysun/GetBaseCountsMultiSample/GetBaseCountsMultiSample",
                          dmp.dir = "/juno/work/access/production/resources/cbioportal/current/mskimpact", mirror.bam.dir = "/juno/res/dmpcollab/dmpshare/share/irb12_245",
                          dmp.key.path = "/juno/res/dmpcollab/dmprequest/12-245/key.txt") {
  # # test input section -----------------------------------------------------------
  # master.ref = fread('/juno/work/bergerm1/bergerlab/zhengy1/access_data_analysis/data/example_master_file.csv')
  # results.dir = paste0('/juno/work/bergerm1/MSK-ACCESS/ACCESS-Projects/test_access/access_data_analysis/output_',format(Sys.time(),'%m%d%y'))
  # pooled.bam.dir = '/ifs/work/bergerm1/ACCESS-Projects/novaseq_curated_duplex_v2/'
  # fasta.path = '/work/access/production/resources/reference/current/Homo_sapiens_assembly19.fasta'
  # genotyper.path = '/ifs/work/bergerm1/Innovation/software/maysun/GetBaseCountsMultiSample/GetBaseCountsMultiSample'
  # dmp.dir = '/ifs/work/bergerm1/zhengy1/dmp/mskimpact/'
  # mirror.bam.dir = '/ifs/dmpshare/share/irb12_245/'
  # dmp.key.path = '/ifs/dmprequest/12-245/key.txt'
  # setting up directory ----------------------------------------------------
  dir.create(results.dir)
  # make tmp directory in output directory
  dir.create(paste0(results.dir, "/tmp"))
  # checking virtualenv -----------------------------------------------------
  geno.bash <- system("which genotype_variants", intern = T)
  if (length(geno.bash) == 0) {
    # print(pyclone.path)
    stop("needs to run \nsource /home/accessbot/miniconda3/bin/activate && conda activate genotype-variants-0.3.0")
  }

  # data from DMP -----------------------------------------------------------
  DMP.key <- fread(dmp.key.path)
  if (any(!master.ref[grepl("^P-", dmp_patient_id)]$dmp_patient_id %in% gsub("-T..-IM.", "", DMP.key[grepl("IM", V1)]$V1))) {
    stop(paste0(
      "These DMP IDs are not found in DMP key file: ",
      paste0(master.ref[grepl("^P-", dmp_patient_id)]$dmp_patient_id[which(!master.ref[grepl("^P-", dmp_patient_id)]$dmp_patient_id %in%
        gsub("-T..-IM.", "", DMP.key[grepl("IM", V1)]$V1))], collapse = " ,")
    ))
  }
  DMP.maf <- fread(paste0(dmp.dir, "/data_mutations_extended.txt")) %>%
    filter(Mutation_Status != "GERMLINE") %>%
    data.table()
  DMP.RET.maf <- DMP.maf[grepl(paste0(unique(master.ref[grepl("^P-", dmp_patient_id)]$dmp_patient_id), collapse = "|"), Tumor_Sample_Barcode), ]

  # Pooled normal samples ---------------------------------------------------
  pooled.bams <- list.files(pooled.bam.dir, pattern = ".bam", full.names = T)

  # For each patient --------------------------------------------------------
  x <- unique(master.ref$cmo_patient_id)[1]
  # x = unique(master.ref$cmo_sample_id_plasma)[16]
  # x = 'C-YW82CY'
  print("Compiling reads per patient")
  all.fillout.id <- lapply(unique(master.ref$cmo_patient_id), function(x) {
    print(x)
    dir.create(paste0(results.dir, "/", x))
    dmp_id <- unique(master.ref[cmo_patient_id == x]$dmp_patient_id)
    # sample sheet with colummns -- TSB, sample type, bam path, treatm --------
    # need to get DMP tumor, DMP normal, plasma, plasma normal (if there is any), pooled normal
    # DMP sample sheet
    if (!is.na(dmp_id)) {
      all.dmp.ids <- DMP.key[grepl(paste0(dmp_id, "-T..-IM."), V1)]$V1
      all.dmp.bam.ids <- DMP.key[grepl(paste0(dmp_id, "-T..-IM."), V1)]$V2
      bam.sub.dir <- unlist(lapply(strsplit(substr(all.dmp.bam.ids, 1, 2), ""), function(x) {
        paste0(x, collapse = "/")
      }))
      dmp.sample.sheet <- data.frame(
        Sample_Barcode = all.dmp.ids,
        standard_bam = paste0(mirror.bam.dir, "/", bam.sub.dir, "/", all.dmp.bam.ids, ".bam") %>%
          mutate(cmo_patient_id = x, Sample_Type = ifelse(grepl("-T", Sample_Barcode), "DMP_Tumor", "DMP_Normal"), dmp_patient_id = dmp_id)
      )
    } else {
      dmp.sample.sheet <- NULL
    }
    # total sample sheet
    sample.sheet <- master.ref[
      cmo_patient_id == x,
      # plasma bams -- duplex and simplex bam
      .(
        Sample_Barcode = cmo_sample_id_plasma, duplex_bam = bam_path_plasma_duplex,
        simplex_bam = bam_path_plasma_simplex, cmo_patient_id, Sample_Type = "duplex", dmp_patient_id
      )
    ] %>%
      merge(rbind(
        unique(master.ref[
          cmo_patient_id == x,
          # buffy coat + DMP bams -- standard bam only
          .(
            Sample_Barcode = cmo_sample_id_normal, standard_bam = bam_path_normal,
            cmo_patient_id, Sample_Type = "unfilterednormal", dmp_patient_id
          )
        ]),
        dmp.sample.sheet
      ), all = T)
    # catch '' or NA for empty cells for some cmo_sample_id_normal
    sample.sheet <- sample.sheet[-which(is.na(Sample_Barcode) | Sample_Barcode == "")]
    write.table(sample.sheet, paste0(results.dir, "/", x, "/", x, "_sample_sheet.tsv"), sep = "\t", quote = F, row.names = F)
    # piece together all unique calls -----------------------------------------
    # get duplex calls
    duplex.calls <- do.call(rbind, lapply(master.ref[cmo_patient_id == x]$maf_path, function(x) {
      # fread(x) %>% filter(as.numeric(D_t_alt_count_fragment) > 0) %>% data.table()
      fread(x) %>%
        filter(as.numeric(t_alt_count) > 0) %>%
        data.table()
    }))
    # get impact calls
    impact.calls <- DMP.RET.maf[Tumor_Sample_Barcode %in% sample.sheet$Sample_Barcode]
    write.table(impact.calls[, .(Hugo_Symbol, Chromosome, Start_Position, End_Position, Variant_Classification, HGVSp_Short, Reference_Allele, Tumor_Seq_Allele2)],
      paste0(results.dir, "/", x, "/", x, "_impact_calls.maf"),
      sep = "\t", quote = F, row.names = F
    )
    # combining plasma and impact calls
    all.calls <- rbind(
      duplex.calls[, intersect(colnames(duplex.calls), colnames(DMP.RET.maf)), with = F],
      impact.calls[, intersect(colnames(duplex.calls), colnames(DMP.RET.maf)), with = F]
    )
    # getting rid of duplicate calls and take the first occurence of all events
    all.calls <- all.calls[which(!duplicated(all.calls[, .(Hugo_Symbol, Chromosome, Start_Position, End_Position, Variant_Classification, HGVSp_Short, Reference_Allele, Tumor_Seq_Allele2)])), ] %>%
      filter(Variant_Classification != "Silent" & !grepl("RP11-", Hugo_Symbol) & !grepl("Intron", Variant_Classification))
    write.table(all.calls, paste0(results.dir, "/", x, "/", x, "_all_unique_calls.maf"), sep = "\t", quote = F, row.names = F)
    # tagging hotspots
    system(paste0(
      'bsub  -R "rusage[mem=4]" -cwd ', results.dir, "/", x, "/ -oo hotspot.o -eo hotspot.e -W 00:59 ",
      " -P ", project.ID, " -J ", x, "_tag_hotspot ",
      " python /ifs/work/bergerm1/zhengy1/ACCESS-Pipeline-DEV/cwl_tools/hotspots/tag_hotspots.py ",
      " -m ", results.dir, "/", x, "/", x, "_all_unique_calls.maf",
      " -itxt /ifs/work/bergerm1/Innovation/Resources/Hotspots/hotspot-list-union-v1-v2_with_TERT.txt ",
      " -o ", results.dir, "/", x, "/", x, "_all_unique_calls_hotspots.maf",
      " -outdir ", results.dir, "/", x, "/", x
    ))

    # genotype all bams in this patient directory -----------------------------
    # genotyping plasma samples -- plasma duplex&simplex, plasma normal, pooled plasma normal
    write.table(sample.sheet[, .(
      sample_id = Sample_Barcode, maf = paste0(results.dir, "/", x, "/", x, "_all_unique_calls.maf"),
      standard_bam, duplex_bam, simplex_bam
    )],
    paste0(results.dir, "/", x, "/", x, "_genotype_metadata.tsv"),
    sep = "\t", quote = F, row.names = F
    )
    job.ids <- system(paste0(
      "bsub -cwd ", results.dir, "/", x, ' -W 12:00  -R "rusage[mem=8]"  -oo genotyping.o -eo genotyping.e ',
      " -P ", project.ID, " -J ", x, "_genotype_variants ",
      " genotype_variants small_variants multiple-samples -i ", results.dir, "/", x, "/", x, "_genotype_metadata.tsv",
      " -r ", fasta.path, " -g ", genotyper.path, " -v DEBUG "
    ), intern = T)
    job.ids <- as.numeric(gsub("Job <|> is.*.$", "", job.ids))
  })


  # Get base count multi sample in pooled normal ----------------------------
  # all all unique calls in entire cohort
  print("Compiling reads in pooled samples")
  dir.create(paste0(results.dir, "/pooled"))
  all.all.unique.mafs <- do.call(rbind, lapply(unique(master.ref$cmo_patient_id), function(x) {
    fread(list.files(paste0(results.dir, "/", x), pattern = "unique_calls.maf$", full.names = T))
  }))
  all.all.unique.mafs <- all.all.unique.mafs[-which(duplicated(all.all.unique.mafs[, .(Hugo_Symbol, Chromosome, Start_Position, End_Position, Variant_Classification, HGVSp_Short, Reference_Allele, Tumor_Seq_Allele2)]))]
  write.table(all.all.unique.mafs, paste0(results.dir, "/pooled/all_all_unique.maf"), sep = "\t", quote = F, row.names = F)

  write.table(data.frame(
    sample_id = gsub("^.*./|.bam", "", pooled.bams), maf = paste0(results.dir, "/pooled/all_all_unique.maf"),
    standard_bam = pooled.bams, duplex_bam = "", simplex_bam = ""
  ),
  paste0(results.dir, "/pooled/pooled_metadata.tsv"),
  sep = "\t", quote = F, row.names = F
  )

  pooled.sample.job.id <- system(paste0(
    "bsub -cwd ", results.dir, '/pooled -W 12:00  -R "rusage[mem=8]" -oo genotyping.o -eo genotyping.e ',
    " -w ", ' \"', paste0(paste0("done(", unlist(all.fillout.id), ")"), collapse = "&&"), '\" ',
    " -P ", project.ID, " -J pooled_genotype_variants ",
    " genotype_variants small_variants multiple-samples -i ", results.dir, "/pooled/pooled_metadata.tsv",
    " -r ", fasta.path, " -g ", genotyper.path, " -v DEBUG "
  ), intern = T)
  pooled.sample.job.id <- as.numeric(gsub("Job <|> is.*.$", "", pooled.sample.job.id))
  while (
    !any(grepl("Done successfully", system(paste0("bjobs -l ", pooled.sample.job.id), intern = T)))
  ) {
    Sys.sleep(120)
  }
  print("Compile reads done!")
}

# Executable -----------------------------------------------------------------------------------------------------------
# Minimal columns for input mafs
#
# Hugo_Symbol,Chromosome,Start_Position,End_Position,Tumor_Sample_Barcode,Variant_Classification,HGVSp_Short,Reference_Allele,Tumor_Seq_Allele2,D_t_alt_count_fragment

suppressPackageStartupMessages({
  library(data.table)
  library(tidyr)
  library(stringr)
  library(dplyr)
  library(argparse)
})

if (!interactive()) {
  parser <- ArgumentParser()
  parser$add_argument("-m", "--masterref", type = "character", help = "File path to master reference file")
  parser$add_argument("-o", "--resultsdir", type = "character", help = "Output directory")
  parser$add_argument("-pid", "--projectid",
    type = "character", default = "",
    help = "Project ID for submitted jobs involved in this run"
  )
  parser$add_argument("-pb", "--pooledbamdir",
    type = "character", default = "/juno/work/access/production/resources/msk-access/current/novaseq_curated_duplex_bams_dmp/current/",
    help = "Directory for all pooled bams [default]"
  )
  parser$add_argument("-fa", "--fastapath",
    type = "character", default = "/juno/work/access/production/resources/reference/current/Homo_sapiens_assembly19.fasta",
    help = "Reference fasta path [default]"
  )
  parser$add_argument("-gt", "--genotyperpath",
    type = "character", default = "/ifs/work/bergerm1/Innovation/software/maysun/GetBaseCountsMultiSample/GetBaseCountsMultiSample",
    help = "Genotyper executable path [default]"
  )
  parser$add_argument("-dmp", "--dmpdir",
    type = "character", default = "/juno/work/access/production/resources/cbioportal/current/mskimpact",
    help = "Directory of clinical DMP IMPACT repository [default]"
  )
  parser$add_argument("-mb", "--mirrorbamdir",
    type = "character", default = "/juno/res/dmpcollab/dmpshare/share/irb12_245",
    help = "Mirror BAM file directory [default]"
  )
  parser$add_argument("-dmpk", "--dmpkeypath",
    type = "character", default = "/juno/res/dmpcollab/dmprequest/12-245/key.txt",
    help = "DMP mirror BAM key file [default]"
  )
  args <- parser$parse_args()

  master.ref <- args$masterref
  results.dir <- args$resultsdir
  project.ID <- args$projectid
  pooled.bam.dir <- args$pooledbamdir
  fasta.path <- args$fastapath
  genotyper.path <- args$genotyperpath
  dmp.dir <- args$dmpdir
  mirror.bam.dir <- args$mirrorbamdir
  dmp.key.path <- args$dmpkeypath


  if (project.ID == "") {
    project.ID <- paste0(sample(c(0:9), size = 10, replace = T), collapse = "")
  }

  print(paste0("Input parameters for run ", project.ID))
  print(master.ref)
  print(results.dir)
  print(pooled.bam.dir)
  print(fasta.path)
  print(genotyper.path)
  print(dmp.dir)
  print(mirror.bam.dir)
  print(dmp.key.path)
  suppressWarnings(compile_reads(fread(master.ref), results.dir, project.ID, pooled.bam.dir, fasta.path, genotyper.path, dmp.dir, mirror.bam.dir, dmp.key.path))
  print("compile reads function finished")
}
