# library(data.table)
# library(stringr)
# library(tidyr)
# library(dplyr)

#' @export
CNA_processing = function(
  master.ref, results.dir, dmp.dir = '/work/access/production/resources/cbioportal/current/msk_solid_heme/',
  gene_set = "v1"
) {
  # Define gene sets
  access.cna.genes = c('AR','BRCA1','BRCA2','CDK4','EGFR','ERBB2','MET','MDM2','MLH1','MSH2','MSH6','MYC')
  access.v2.cna.genes = c(
    "ALK", "APC", "AR", "ARID1A", "ASXL1", "ATM", "BAP1", "BRCA1", "BRCA2", "CDH1", "CDK12", "CDK4", "CHEK2",
    "DNMT3A", "EGFR", "ERBB2", "ERBB3", "ERCC2", "FBXW7", "FGFR2", "FGFR3", "KDM6A", "KIT", "MAP2K1", "MAP2K2",
    "MDM2", "MET", "MLH1", "MSH2", "MSH6", "MTOR", "MYC", "MYCN", "NF1", "PALB2", "PDGFRA", "PIK3CA", "PMS2",
    "PTCH1", "RB1", "RET", "SMAD4", "TP53", "TSC1", "TSC2"
  )

  # Select the appropriate gene set
  selected_genes = if (gene_set == "v2") access.v2.cna.genes else access.cna.genes

  # DMP CNA calls --------------------------------------------------------
  DMP.cna <- fread(paste0(dmp.dir, '/data_CNA.txt'))

  dmp.ids <- master.ref$dmp_patient_id[(master.ref$dmp_patient_id != "") & !is.na(master.ref$dmp_patient_id)]
  dmp.pattern <- paste0(dmp.ids, collapse = '|')

  DMP.cna.edited = DMP.cna[
    !Hugo_Symbol %in% c('CDKN2Ap14ARF', 'CDKN2Ap16INK4A'),
    c('Hugo_Symbol', grep(dmp.pattern, colnames(DMP.cna), value = T)),
    with = F
  ] %>%
    melt.data.table(
      id.vars = 'Hugo_Symbol',
      variable.name = 'Tumor_Sample_Barcode',
      value.name = 'fc'
    ) %>%
    filter(fc != 0) %>%
    mutate(
      CNA_tumor = case_when(fc < 0 ~ 'HOMDEL', fc > 0 ~ 'AMP'),
      dmp_patient_id = gsub('-T..-IM.', '', Tumor_Sample_Barcode)
    ) %>%
    merge(
      unique(master.ref[!is.na(dmp_patient_id), .(cmo_patient_id, dmp_patient_id)]),
      by = 'dmp_patient_id',
      all.x = T
    ) %>%
    transmute(Hugo_Symbol, CNA_tumor, cmo_patient_id, dmp_patient_id) %>%
    unique() %>%
    data.table()

  # Access CNA calls --------------------------------------------------------
  access.cna = do.call(rbind, lapply(master.ref$cna_path, function(x) {
    fread(x) %>%
      transmute(
        Tumor_Sample_Barcode = gsub('\\.', '-', gsub('_mean_cvg', '', sample)),
        Hugo_Symbol = region,
        p.adj,
        fc
      ) %>%
      filter(p.adj <= 0.05) %>%
      merge(
        unique(master.ref[, .(cmo_patient_id, cmo_sample_id_plasma)]),
        by.x = 'Tumor_Sample_Barcode',
        by.y = 'cmo_sample_id_plasma',
        all.x = T
      ) %>%
      data.table()
  }))

  merge(
    access.cna,
    DMP.cna.edited,
    by = c('Hugo_Symbol', 'cmo_patient_id'),
    all.x = T
  ) %>%
    mutate(
      CNA = case_when(
        (fc >= 1.5 & Hugo_Symbol %in% selected_genes) ~ 'AMP',
        (fc <= -1.5 & Hugo_Symbol %in% selected_genes) ~ 'HOMDEL',
        (fc >= 1.2 & !is.na(CNA_tumor)) ~ 'AMP',
        (fc <= -1.2 & !is.na(CNA_tumor)) ~ 'HOMDEL',
        TRUE ~ ''
      )
    ) %>%
    filter(CNA != '') %>%
    data.table() -> access.call.set

  # Compare tumor and plasma CNAs
  condition <- any(access.call.set$CNA_tumor != access.call.set$CNA)
  if (condition || is.na(condition)) warning('Tumor and plasma have conflicting CNA detected!')

  # Save results
  dir.create(paste0(results.dir, '/CNA_final_call_set'), showWarnings = FALSE)
  lapply(unique(master.ref$cmo_sample_id_plasma), function(x) {
    write.table(access.call.set[Tumor_Sample_Barcode == x],
                paste0(results.dir, '/CNA_final_call_set/', x, '_cna_final_call_set.txt'),
                sep = '\t', quote = F, row.names = F)
  })
}
# Executable -----------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(stringr)
  library(dplyr)
})

if (!interactive()) {

  parser = ArgumentParser()
  parser$add_argument('-m', '--masterref', type = 'character', help = 'File path to master reference file')
  parser$add_argument('-o', '--resultsdir', type = 'character', help = 'Output directory')
  parser$add_argument('-dmp', '--dmpdir', type = 'character', default = '/work/access/production/resources/cbioportal/current/mskimpact/',
                      help = 'Directory of clinical DMP IMPACT repository [default]')
  parser$add_argument('-g', '--gene-set', type = 'character', default = 'v1',
                      help = 'Gene set to use: "v1" for access.cna.genes or "v2" for access.v2.cna.genes [default: v1]')
  args = parser$parse_args()

  master.ref = args$masterref
  results.dir = args$resultsdir
  dmp.dir = args$dmpdir
  gene_set = args$gene_set

  CNA_processing(fread(master.ref), results.dir, dmp.dir, gene_set)
}
