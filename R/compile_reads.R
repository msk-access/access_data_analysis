#library(data.table)
#library(tidyr)
#library(stringr)
#library(dplyr)


#' @export
compile_reads <- function(master.ref,
                          results.dir,
                          project.ID,
                          pooled.bam.dir = "/juno/work/access/production/resources/msk-access/current/novaseq_curated_duplex_bams_dmp/current/",
                          fasta.path = "/juno/work/access/production/resources/reference/current/Homo_sapiens_assembly19.fasta",
                          genotyper.path = "/work/access/production/resources/tools/GetBaseCountsMultiSample/current/GetBaseCountsMultiSample",
                          dmp.dir = "/juno/work/access/production/resources/cbioportal/current/msk_solid_heme",
                          mirror.bam.dir = "/juno/res/dmpcollab/dmpshare/share/irb12_245",
                          mirror.access.bam.dir = "/juno/res/dmpcollab/dmpshare/share/access_12_245/",
                          dmp.key.path = "/juno/res/dmpcollab/dmprequest/12-245/key.txt",
                          access.key.path = "/juno/res/dmpcollab/dmprequest/ACCESS-12-245/key.txt") {
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
    stop(
      "needs to run \nsource /home/accessbot/miniconda3/bin/activate && conda activate genotype-variants-0.3.0"
    )
  }

  # data from DMP -----------------------------------------------------------
  DMP.key <- fread(dmp.key.path)
  if (any(
    !master.ref[grepl("^P-", dmp_patient_id)]$dmp_patient_id %in% gsub("-T..-IH.|-T..-IM.|-T..-XS", "", DMP.key[grepl("IH|IM|XS", V1)]$V1)
  )) {
    warning(paste0(
      "These DMP IDs are not found in DMP key file: ",
      paste0(master.ref[grepl("^P-", dmp_patient_id)]$dmp_patient_id[which(
        !master.ref[grepl("^P-", dmp_patient_id)]$dmp_patient_id %in%
          gsub("-T..-IH.|-T..-IM.|-T..-XS", "", DMP.key[grepl("IH|IM|XS", V1)]$V1)
      )], collapse = " ,")
    ))
  }
  access.key <- fread(access.key.path)
  if (any(
    !master.ref[grepl("^P-", dmp_patient_id)]$dmp_patient_id %in% gsub("-T..-IH.|-T..-IM.|-T..-XS", "", access.key[grepl("IH|IM|XS", V1)]$V1)
  )) {
    warning(paste0(
      "These DMP IDs are not found in DMP key file: ",
      paste0(master.ref[grepl("^P-", dmp_patient_id)]$dmp_patient_id[which(
        !master.ref[grepl("^P-", dmp_patient_id)]$dmp_patient_id %in%
          gsub("-T..-IH.|-T..-IM.|-T..-XS", "", access.key[grepl("IH|IM|XS", V1)]$V1)
      )], collapse = " ,")
    ))
  }
  DMP.maf <-
    fread(paste0(dmp.dir, "/data_mutations_extended.txt")) %>%
    filter(Mutation_Status != "GERMLINE") %>%
    data.table()
  DMP.RET.maf <-
    DMP.maf[grepl(paste0(unique(master.ref[grepl("^P-", dmp_patient_id)]$dmp_patient_id), collapse = "|"), Tumor_Sample_Barcode),]

  # Pooled normal samples ---------------------------------------------------
  pooled.bams <-
    list.files(pooled.bam.dir, pattern = ".bam", full.names = T)

  # For each patient --------------------------------------------------------
  x <- unique(master.ref$cmo_patient_id)[1]
  # x = unique(master.ref$cmo_sample_id_plasma)[16]
  # x = 'C-YW82CY'
  print("Compiling reads per patient")
  all.fillout.id <-
    lapply(unique(master.ref$cmo_patient_id), function(x) {
      print(x)
      dir.create(paste0(results.dir, "/", x))
      dmp_id <-
        unique(master.ref[cmo_patient_id == x]$dmp_patient_id)
      # sample sheet with colummns -- TSB, sample type, bam path, treatm --------
      # need to get DMP tumor, DMP normal, plasma, plasma normal (if there is any), pooled normal
      # DMP sample sheet
      if (is.na(dmp_id) | dmp_id == '') {
        dmp.sample.sheet <- NULL
      } else {
        all.dmp.ids.IM <-
          DMP.key[grepl(paste0(dmp_id, "-(T|N)..-IM."), V1)]$V1
        all.dmp.ids.IH <-
          DMP.key[grepl(paste0(dmp_id, "-(T|N)..-IH."), V1)]$V1
        all.dmp.ids.XS <-
          access.key[grepl(paste0(dmp_id, "-(T|N)..-XS."), V1)]$V1
        all.dmp.ids <- c(all.dmp.ids.IM, all.dmp.ids.IH)
        all.dmp.bam.ids.IM <-
          DMP.key[grepl(paste0(dmp_id, "-(T|N)..-IM."), V1)]$V2
        all.dmp.bam.ids.IH <-
          DMP.key[grepl(paste0(dmp_id, "-(T|N)..-IH."), V1)]$V2
        all.dmp.bam.ids.XS <-
          access.key[grepl(paste0(dmp_id, "-(T|N)..-XS."), V1)]$V2
        print(all.dmp.bam.ids.XS)
        all.dmp.bam.ids.XS <-
          gsub("-standard|-unfilter|-simplex|-duplex",
               "",
               access.key[grepl(paste0(dmp_id, "-(T|N)..-XS."), V1)]$V2)
        print(all.dmp.bam.ids.XS)
        all.dmp.bam.ids <-
          c(all.dmp.bam.ids.IM,
            all.dmp.bam.ids.IH)
        bam.sub.dir <-
          unlist(lapply(strsplit(substr(
            all.dmp.bam.ids, 1, 2
          ), ""), function(x) {
            paste0(x, collapse = "/")
          }))
        dmp.sample.sheet <- data.frame(
          Sample_Barcode = all.dmp.ids,
          standard_bam = paste0(
            mirror.bam.dir,
            "/",
            bam.sub.dir,
            "/",
            all.dmp.bam.ids,
            ".bam"
          )
        )
        access.bam.sub.dir <-
          unlist(lapply(strsplit(
            substr(all.dmp.bam.ids.XS, 1, 2), ""
          ), function(x) {
            paste0(x, collapse = "/")
          }))
        access.sample.sheet <- unique(
          data.frame(
            Sample_Barcode = all.dmp.ids.XS,
            duplex_bam = paste0(
              mirror.access.bam.dir,
              "/",
              access.bam.sub.dir,
              "/",
              all.dmp.bam.ids.XS,
              "-duplex.bam"
            ),
            simplex_bam = paste0(
              mirror.access.bam.dir,
              "/",
              access.bam.sub.dir,
              "/",
              all.dmp.bam.ids.XS,
              "-simplex.bam"
            )
          )
        )
        dmp.sample.sheet <-
          bind_row(dmp.sample.sheet, access.sample.sheet)  %>%
          mutate(
            cmo_patient_id = x,
            Sample_Type = ifelse(grepl("-T", Sample_Barcode), "DMP_Tumor", "DMP_Normal"),
            dmp_patient_id = dmp_id
          )
      }
      # total sample sheet
      sample.sheet <- master.ref[cmo_patient_id == x,
                                 # plasma bams -- duplex and simplex bam
                                 .(
                                   Sample_Barcode = as.character(cmo_sample_id_plasma),
                                   duplex_bam = bam_path_plasma_duplex,
                                   simplex_bam = bam_path_plasma_simplex,
                                   cmo_patient_id,
                                   Sample_Type = "duplex",
                                   dmp_patient_id
                                 )] %>%
        merge(rbind(unique(master.ref[cmo_patient_id == x &
                                        paired == 'Paired',
                                      # buffy coat + DMP bams -- standard bam only
                                      .(
                                        Sample_Barcode = as.character(cmo_sample_id_normal),
                                        standard_bam = bam_path_normal,
                                        cmo_patient_id,
                                        Sample_Type = "unfilterednormal",
                                        dmp_patient_id
                                      )]),
                    dmp.sample.sheet), all = T)
      # catch '' or NA for empty cells for some cmo_sample_id_normal
      sample.sheet <-
        sample.sheet[!is.na(Sample_Barcode) | Sample_Barcode != ""]
      write.table(
        sample.sheet,
        paste0(results.dir, "/", x, "/", x, "_sample_sheet.tsv"),
        sep = "\t",
        quote = F,
        row.names = F
      )
      # piece together all unique calls -----------------------------------------
      # get duplex calls
      duplex.calls <-
        do.call(rbind, lapply(master.ref[cmo_patient_id == x]$maf_path, function(x) {
          # fread(x) %>% filter(as.numeric(D_t_alt_count_fragment) > 0) %>% data.table()
          selectcolumns <-
            c(
              "Hugo_Symbol",
              "Entrez_Gene_Id",
              "Center",
              "NCBI_Build",
              "Chromosome",
              "Start_Position",
              "End_Position",
              "Strand",
              "Variant_Classification",
              "Variant_Type",
              "Reference_Allele",
              "Tumor_Seq_Allele1",
              "Tumor_Seq_Allele2",
              "dbSNP_RS",
              "dbSNP_Val_Status",
              "Tumor_Sample_Barcode",
              "caller_Norm_Sample_Barcode",
              "Match_Norm_Seq_Allele1",
              "Match_Norm_Seq_Allele2",
              "Tumor_Validation_Allele1",
              "Tumor_Validation_Allele2",
              "Match_Norm_Validation_Allele1",
              "Match_Norm_Validation_Allele2",
              "Verification_Status",
              "Validation_Status",
              "Mutation_Status",
              "Sequencing_Phase",
              "Sequence_Source",
              "Validation_Method",
              "Score",
              "BAM_File",
              "Sequencer",
              "Tumor_Sample_UUID",
              "Matched_Norm_Sample_UUID",
              "HGVSc",
              "HGVSp",
              "HGVSp_Short",
              "Transcript_ID",
              "Exon_Number",
              "caller_t_depth",
              "caller_t_ref_count",
              "caller_t_alt_count",
              "caller_n_depth",
              "caller_n_ref_count",
              "caller_n_alt_count",
              "all_effects",
              "Allele",
              "Gene",
              "Feature",
              "Feature_type",
              "Consequence",
              "cDNA_position",
              "CDS_position",
              "Protein_position",
              "Amino_acids",
              "Codons",
              "Existing_variation",
              "ALLELE_NUM",
              "DISTANCE",
              "STRAND_VEP",
              "SYMBOL",
              "SYMBOL_SOURCE",
              "HGNC_ID",
              "BIOTYPE",
              "CANONICAL",
              "CCDS",
              "ENSP",
              "SWISSPROT",
              "TREMBL",
              "UNIPARC",
              "RefSeq",
              "SIFT",
              "PolyPhen",
              "EXON",
              "INTRON",
              "DOMAINS",
              "AF",
              "AFR_AF",
              "AMR_AF",
              "ASN_AF",
              "EAS_AF",
              "EUR_AF",
              "SAS_AF",
              "AA_AF",
              "EA_AF",
              "CLIN_SIG",
              "SOMATIC",
              "PUBMED",
              "MOTIF_NAME",
              "MOTIF_POS",
              "HIGH_INF_POS",
              "MOTIF_SCORE_CHANGE",
              "IMPACT",
              "PICK",
              "VARIANT_CLASS",
              "TSL",
              "HGVS_OFFSET",
              "PHENO",
              "MINIMISED",
              "ExAC_AF",
              "ExAC_AF_AFR",
              "ExAC_AF_AMR",
              "ExAC_AF_EAS",
              "ExAC_AF_FIN",
              "ExAC_AF_NFE",
              "ExAC_AF_OTH",
              "ExAC_AF_SAS",
              "GENE_PHENO",
              "FILTER",
              "flanking_bps",
              "variant_id",
              "variant_qual",
              "ExAC_AF_Adj",
              "ExAC_AC_AN_Adj",
              "ExAC_AC_AN",
              "ExAC_AC_AN_AFR",
              "ExAC_AC_AN_AMR",
              "ExAC_AC_AN_EAS",
              "ExAC_AC_AN_FIN",
              "ExAC_AC_AN_NFE",
              "ExAC_AC_AN_OTH",
              "ExAC_AC_AN_SAS",
              "ExAC_FILTER",
              "gnomAD_AF",
              "gnomAD_AFR_AF",
              "gnomAD_AMR_AF",
              "gnomAD_ASJ_AF",
              "gnomAD_EAS_AF",
              "gnomAD_FIN_AF",
              "gnomAD_NFE_AF",
              "gnomAD_OTH_AF",
              "gnomAD_SAS_AF",
              "CallMethod",
              "VCF_POS",
              "VCF_REF",
              "VCF_ALT",
              "hotspot_whitelist",
              "Status",
              "D_t_alt_count_fragment",
              "D_t_ref_count_fragment",
              "D_t_vaf_fragment",
              "SD_t_alt_count_fragment",
              "SD_t_ref_count_fragment",
              "SD_t_vaf_fragment",
              "Matched_Norm_Sample_Barcode",
              "Matched_Norm_Bamfile",
              "n_alt_count_fragment",
              "n_ref_count_fragment",
              "n_vaf_fragment"
            )
          if ("Status" %in% names(fread(x))) {
            fread(x) %>% select(one_of(selectcolumns)) %>% subset((Status == "") |
                                                                    (is.na(Status)))
          } else {
            fread(x) %>% select(one_of(selectcolumns))
          }
          #      fread(x)
          # %>%
          # filter(as.numeric(t_alt_count) > 0) %>%
          # data.table()
        }))
      # get impact calls
      impact.calls <-
        DMP.RET.maf[Tumor_Sample_Barcode %in% sample.sheet$Sample_Barcode]
      write.table(
        impact.calls[, .(
          Hugo_Symbol,
          Chromosome,
          Start_Position,
          End_Position,
          Variant_Classification,
          HGVSp_Short,
          Reference_Allele,
          Tumor_Seq_Allele2
        )],
        paste0(results.dir, "/", x, "/", x, "_impact_calls.maf"),
        sep = "\t",
        quote = F,
        row.names = F
      )
      # combining plasma and impact calls
      all.calls <-
        rbind(duplex.calls[, intersect(colnames(duplex.calls), colnames(DMP.RET.maf)), with = F],
              impact.calls[, intersect(colnames(duplex.calls), colnames(DMP.RET.maf)), with = F])
      # getting rid of duplicate calls and take the first occurence of all events
      all.calls <-
        all.calls[which(!duplicated(all.calls[, .(
          Hugo_Symbol,
          Chromosome,
          Start_Position,
          End_Position,
          Variant_Classification,
          HGVSp_Short,
          Reference_Allele,
          Tumor_Seq_Allele2
        )])),] %>%
        mutate(
          t_ref_count = 0,
          t_alt_count = 0,
          n_ref_count = 0,
          n_alt_count = 0,
          Matched_Norm_Sample_Barcode = NA
        ) %>%
        filter(
          Variant_Classification != "Silent" &
            !grepl("RP11-", Hugo_Symbol) &
            !grepl("Intron", Variant_Classification)
        )
      write.table(
        all.calls,
        paste0(results.dir, "/", x, "/", x, "_all_unique_calls.maf"),
        sep = "\t",
        quote = F,
        row.names = F
      )
      # tagging hotspots
      system(
        paste0(
          'bsub  -R "rusage[mem=4]" -cwd ',
          results.dir,
          "/",
          x,
          "/ -oo hotspot.o -eo hotspot.e -W 00:59 ",
          " -P ",
          project.ID,
          " -J ",
          x,
          "_tag_hotspot ",
          " python /work/access/production/workflows/access_workflows/v1/pipeline_2.0.0/ACCESS-Pipeline/cwl_tools/hotspots/tag_hotspots.py ",
          " -m ",
          results.dir,
          "/",
          x,
          "/",
          x,
          "_all_unique_calls.maf",
          " -itxt /work/access/production/resources/msk-access/current/regions_of_interest/current/hotspot-list-union-v1-v2_with_TERT.txt ",
          " -o ",
          results.dir,
          "/",
          x,
          "/",
          x,
          "_all_unique_calls_hotspots.maf",
          " -outdir ",
          results.dir,
          "/",
          x,
          "/",
          x
        )
      )
      # genotype all bams in this patient directory -----------------------------
      # genotyping plasma samples -- plasma duplex&simplex, plasma normal, pooled plasma normal
      write.table(
        sample.sheet[, .(
          sample_id = Sample_Barcode,
          maf = paste0(results.dir, "/", x, "/", x, "_all_unique_calls.maf"),
          standard_bam,
          duplex_bam,
          simplex_bam
        )],
        paste0(results.dir, "/", x, "/", x, "_genotype_metadata.tsv"),
        sep = "\t",
        quote = F,
        row.names = F
      )
      job.ids <- system(
        paste0(
          "bsub -cwd ",
          results.dir,
          "/",
          x,
          ' -W 12:00  -R "rusage[mem=8]"  -oo genotyping.o -eo genotyping.e ',
          " -P ",
          project.ID,
          " -J ",
          x,
          "_genotype_variants ",
          " genotype_variants small_variants multiple-samples -i ",
          results.dir,
          "/",
          x,
          "/",
          x,
          "_genotype_metadata.tsv",
          " -r ",
          fasta.path,
          " -g ",
          genotyper.path,
          " -v DEBUG "
        ),
        intern = T
      )
      job.ids <- as.numeric(gsub("Job <|> is.*.$", "", job.ids))
    })


  # Get base count multi sample in pooled normal ----------------------------
  # all all unique calls in entire cohort
  print("Compiling reads in pooled samples")
  dir.create(paste0(results.dir, "/pooled"))
  all.all.unique.mafs <-
    do.call(rbind, lapply(unique(master.ref$cmo_patient_id), function(x) {
      fread(list.files(
        paste0(results.dir, "/", x),
        pattern = "unique_calls.maf$",
        full.names = T
      ))
    }))
  all.all.unique.mafs <-
    all.all.unique.mafs[!duplicated(all.all.unique.mafs[, .(
      Hugo_Symbol,
      Chromosome,
      Start_Position,
      End_Position,
      Variant_Classification,
      HGVSp_Short,
      Reference_Allele,
      Tumor_Seq_Allele2
    )]), ]
  write.table(
    all.all.unique.mafs,
    paste0(results.dir, "/pooled/all_all_unique.maf"),
    sep = "\t",
    quote = F,
    row.names = F
  )

  write.table(
    data.frame(
      sample_id = gsub("^.*./|.bam", "", pooled.bams),
      maf = paste0(results.dir, "/pooled/all_all_unique.maf"),
      standard_bam = pooled.bams,
      duplex_bam = "",
      simplex_bam = ""
    ),
    paste0(results.dir, "/pooled/pooled_metadata.tsv"),
    sep = "\t",
    quote = F,
    row.names = F
  )

  pooled.sample.job.id <- system(
    paste0(
      "bsub -cwd ",
      results.dir,
      '/pooled -W 12:00  -R "rusage[mem=8]" -oo genotyping.o -eo genotyping.e ',
      " -w ",
      ' \"',
      paste0(paste0("done(", unlist(all.fillout.id), ")"), collapse = "&&"),
      '\" ',
      " -P ",
      project.ID,
      " -J pooled_genotype_variants ",
      " genotype_variants small_variants multiple-samples -i ",
      results.dir,
      "/pooled/pooled_metadata.tsv",
      " -r ",
      fasta.path,
      " -g ",
      genotyper.path,
      " -v DEBUG "
    ),
    intern = T
  )
  pooled.sample.job.id <-
    as.numeric(gsub("Job <|> is.*.$", "", pooled.sample.job.id))
  while (!any(grepl("Done successfully", system(
    paste0("bjobs -l ", pooled.sample.job.id), intern = T
  )))) {
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
  parser$add_argument(
    "-pid",
    "--projectid",
    type = "character",
    default = "",
    help = "Project ID for submitted jobs involved in this run"
  )
  parser$add_argument(
    "-pb",
    "--pooledbamdir",
    type = "character",
    default = "/juno/work/access/production/resources/msk-access/current/novaseq_curated_duplex_bams_dmp/current/",
    help = "Directory for all pooled bams [default]"
  )
  parser$add_argument(
    "-fa",
    "--fastapath",
    type = "character",
    default = "/juno/work/access/production/resources/reference/current/Homo_sapiens_assembly19.fasta",
    help = "Reference fasta path [default]"
  )
  parser$add_argument(
    "-gt",
    "--genotyperpath",
    type = "character",
    default = "/work/access/production/resources/tools/GetBaseCountsMultiSample/current/GetBaseCountsMultiSample",
    help = "Genotyper executable path [default]"
  )
  parser$add_argument(
    "-dmp",
    "--dmpdir",
    type = "character",
    default = "/juno/work/access/production/resources/cbioportal/current/msk_solid_heme",
    help = "Directory of clinical DMP repository [default]"
  )
  parser$add_argument(
    "-mb",
    "--mirrorbamdir",
    type = "character",
    default = "/juno/res/dmpcollab/dmpshare/share/irb12_245",
    help = "Mirror BAM file directory [default]"
  )
  parser$add_argument(
    "-mab",
    "--mirroraccessbamdir",
    type = "character",
    default = "/juno/res/dmpcollab/dmpshare/share/access_12_245",
    help = "Mirror BAM file directory for MSK-ACCESS [default]"
  )
  parser$add_argument(
    "-dmpk",
    "--dmpkeypath",
    type = "character",
    default = "/juno/res/dmpcollab/dmprequest/12-245/key.txt",
    help = "DMP mirror BAM key file [default]"
  )
  parser$add_argument(
    "-dmpak",
    "--dmpaccesskeypath",
    type = "character",
    default = " /juno/res/dmpcollab/dmprequest/ACCESS-12-245/key.txt",
    help = "DMP mirror BAM key file for MSK-ACCESS [default]"
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
  mirror.access.bam.dir <- args$mirroraccessbamdir
  dmp.key.path <- args$dmpkeypath
  access.key.path <- args$dmpaccesskeypath


  if (project.ID == "") {
    project.ID <-
      paste0(sample(c(0:9), size = 10, replace = T), collapse = "")
  }

  print(paste0("Input parameters for run ", project.ID))
  print(master.ref)
  print(results.dir)
  print(pooled.bam.dir)
  print(fasta.path)
  print(genotyper.path)
  print(dmp.dir)
  print(mirror.bam.dir)
  print(mirror.access.bam.dir)
  print(dmp.key.path)
  print(access.key.path)
  suppressWarnings(
    compile_reads(
      fread(master.ref),
      results.dir,
      project.ID,
      pooled.bam.dir,
      fasta.path,
      genotyper.path,
      dmp.dir,
      mirror.bam.dir,
      dmp.key.path
    )
  )
  print("compile reads function finished")
}
