# library(data.table)
# library(stringr)
# library(tidyr)
# library(dplyr)
# library(ggplot2)
# library(ggpubr)
# library(RColorBrewer)


# helper methods ----------------------------------------------------------

# collapsing read counts from the same rearrangement events into one total count
collapse_AF <- function(x) {
  # x = c("32-117-190(0.7842)|53-0-959(0.0553)","NA|16-0-1035(0.0155)","63-0-954(0.066)|NA" )
  # x = c('3-5-1300')
  # x = c('NA|NA|NA|0-56-4183','NA|NA|NA|3-4-76')
  print(paste0(x, collapse = "','"))
  # number of samples
  samples.num <- attr(gregexpr("\\|", x[1])[[1]], "match.length")
  print(samples.num)
  # when not found, return -1
  if (samples.num != -1) {
    paste0(round(apply(separate(data.frame(AF = x), "AF", paste0("timpoint_", c(1:(length(samples.num) + 1))), sep = "\\|"), 2, function(one.tp.afs) {
      mean(unlist(lapply(one.tp.afs, function(tmp.af) {
        if (tmp.af == "NA") {
          return(0)
        }
        else {
          read.counts <- as.numeric(str_split(gsub("\\(.*", "", tmp.af), "-")[[1]])
          return((read.counts[1] + read.counts[2]) / read.counts[3])
        }
      })))
    }), digits = 3), collapse = "|")
  } else {
    print(x)
    as.character(round(mean(unlist(lapply(x, function(tmp.af) {
      if (tmp.af == "NA") {
        return(0)
      }
      else {
        read.counts <- as.numeric(str_split(gsub("\\(.*", "", tmp.af), "-")[[1]])
        return((read.counts[1] + read.counts[2]) / read.counts[3])
      }
    }))), digits = 3))
  }
}

# convert naming to timepoint, get rid of uncovered impact and access calls
process_maf_for_graph <- function(tmp.maf) {
  print("convert naming to timepoint, get rid of uncovered impact and access calls")
  # tmp.maf = ret.054.calls
  # tumor sample
  tumor.sample <- structure(gsub("-", "", str_extract(unique(tmp.maf$Tumor_Sample_Barcode[grep("IM[0-9]$", tmp.maf$Tumor_Sample_Barcode)]), "-T..-")),
    names = as.character(unique(tmp.maf$Tumor_Sample_Barcode[grep("IM[0-9]$", tmp.maf$Tumor_Sample_Barcode)]))
  )
  print(tumor.sample)
  # rest of the samples are plasma
  plasma.sample <- setdiff(tmp.maf$Tumor_Sample_Barcode, names(tumor.sample))
  # filter for plasma sample only
  tmp.maf <- tmp.maf[Tumor_Sample_Barcode %in% plasma.sample]
  # change samples into timepoint information
  plasma.sample <- structure(case_when(
    # some of the DA-ret sample need to be renamed
    grepl("-T0._", plasma.sample) ~ gsub("-|_", "", gsub("T", "L0", str_extract(plasma.sample, "-T0._"))),
    # otherwise extract L00 something
    TRUE ~ gsub("-", "", str_extract(plasma.sample, "-L...-"))
  ), names = plasma.sample)
  print(plasma.sample)
  sample.name.conversion <- c(tumor.sample, plasma.sample)
  print(sample.name.conversion)
  # get all not covered calls
  not.covered.df <- unique(tmp.maf[call_confidence == "Not Covered", .N, .(
    Hugo_Symbol, Chromosome, Start_Position, End_Position, Variant_Classification,
    HGVSp_Short, Reference_Allele, Tumor_Seq_Allele2
  )])[N > length(plasma.sample) / 2]
  only.covered.tmp.maf <- anti_join(tmp.maf, not.covered.df, by = c(
    "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Classification",
    "HGVSp_Short", "Reference_Allele", "Tumor_Seq_Allele2"
  )) %>% data.table()
  # only the converted timepoint names
  # only.covered.tmp.maf$Tumor_Sample_Barcode = sample.name.conversion[as.character(only.covered.tmp.maf$Tumor_Sample_Barcode)]
  if (any(grepl("__", only.covered.tmp.maf$Hugo_Symbol))) {
    fusion.only.covered.tmp.maf <- data.table(only.covered.tmp.maf)[grepl("__", Hugo_Symbol)]
    # process original hugo_symbol column (sort two genes by name)
    fusion.only.covered.tmp.maf$Hugo_Symbol <- unlist(lapply(fusion.only.covered.tmp.maf$Hugo_Symbol, function(x) {
      paste0(sort(str_split(x, "__")[[1]]), collapse = "-")
    }))
    fusion.only.covered.tmp.maf$Chromosome <- unlist(lapply(fusion.only.covered.tmp.maf$Chromosome, function(x) {
      paste0(sort(str_split(x, "__")[[1]]), collapse = "-")
    }))
    # collapsing AF for rows of the same events (i.e. reciprocal rearrangement) while perserving the sample level seaparation in AF
    fusion.only.covered.tmp.maf <- fusion.only.covered.tmp.maf[
      , .(
        Start_Position = Start_Position[1], End_Position = End_Position[1], HGVSp_Short = HGVSp_Short[1],
        Reference_Allele = Reference_Allele[1], Tumor_Seq_Allele2 = Tumor_Seq_Allele2[1],
        ExAC_AF = ExAC_AF[1], Hotspot = Hotspot[1], DMP = DMP[1], duplex_support_num = duplex_support_num[1],
        call_confidence = ifelse(any(call_confidence == "Called"), "Called", "Not Called"),
        call_info = paste0(call_info, collapse = " | "), CH = "No",
        t_alt_count = sum(t_alt_count, na.rm = T), t_total_count = sum(t_total_count, na.rm = T)
      ),
      .(Hugo_Symbol, Chromosome, Variant_Classification, Tumor_Sample_Barcode)
    ]
    only.covered.tmp.maf <- only.covered.tmp.maf[-grep("__", Hugo_Symbol)]
    only.covered.tmp.maf <- rbind(only.covered.tmp.maf, fusion.only.covered.tmp.maf)
  }
  only.covered.tmp.maf$t_alt_count <- ifelse(is.na(only.covered.tmp.maf$t_alt_count), 0, only.covered.tmp.maf$t_alt_count)
  only.covered.tmp.maf$t_total_count <- ifelse(is.na(only.covered.tmp.maf$t_total_count), 0, only.covered.tmp.maf$t_total_count)
  return(only.covered.tmp.maf)
}

# melting genotype tables into maf-like format
table_to_maf <- function(tmp.table, sample.table) {
  # tmp.table = fillouts.dt
  # sample.table = sample.sheet
  # tmp.table = ret.006.table
  # sample.table = ret.006.sample.sheet
  # extract information for plasma and tumor
  tmp.table <- data.table(tmp.table)
  lapply(sample.table[Sample_Type %in% c("duplex")]$Sample_Barcode, function(y) {
    sample.call.status.colname <- paste0(y, "___duplex.called")
    sample.af.colname <- paste0(y, "___total")
    tmp.table[, eval(y) := paste0(get(sample.call.status.colname), " | ", get(sample.af.colname))]
  })
  lapply(sample.table[Sample_Type %in% c("DMP_Tumor")]$Sample_Barcode, function(y) {
    tmp.table[, eval(y) := paste0(case_when(
      !is.na(get("DMP")) & get(paste0(sample.table[Sample_Type %in% c("duplex")]$Sample_Barcode[1], "___duplex.called")) != "Not Covered" ~ "Called",
      !is.na(get("DMP")) & get(paste0(sample.table[Sample_Type %in% c("duplex")]$Sample_Barcode[1], "___duplex.called")) == "Not Covered" ~ "Called (but not covered in ACCESS)",
      is.na(get("DMP")) & as.numeric(gsub("/.*", "", get(paste0(y, "___DMP_Tumor")))) > 3 ~ "Genotyped",
      TRUE ~ "Not Called"
    ), " | ", get(paste0(y, "___DMP_Tumor")))]
  })
  processed.tmp.table <- tmp.table[, !grep("___", colnames(tmp.table)), with = F] %>%
    # melting data frame by tumor samples
    melt(
      id.vars = c(
        "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Classification", "HGVSp_Short",
        "Reference_Allele", "Tumor_Seq_Allele2", "ExAC_AF", "Hotspot", "DMP", "duplex_support_num", "call_confidence", "CH"
      ),
      variable.name = "Tumor_Sample_Barcode", value.name = "call_info"
    ) %>%
    mutate(call_confidence = gsub(" \\| ", "", str_extract(call_info, ".*.\\| ")), call_info = gsub(".*.\\| ", "", call_info)) %>%
    rowwise() %>%
    mutate(
      t_alt_count = ifelse(grepl("-[0-9]+-", call_info),
        # SV parsing
        sum(as.numeric(str_split(call_info, "-|\\(")[[1]][1:2])),
        # SNV parsing
        as.numeric(gsub(" |\\/.*.", "", call_info))
      ),
      t_total_count = ifelse(grepl("-[0-9]+-", call_info),
        # SV parsing
        as.numeric(str_split(call_info, "-|\\(")[[1]][3]),
        # SNV parsing
        as.numeric(gsub(".*.\\/|\\(.*.", "", call_info))
      )
    ) %>%
    data.table()
  return(processed.tmp.table)
}

# main graphing function --------------------------------------------------

#' @export
plot_all_events <- function(
                            master.ref, results.dir,
                            criteria = "stringent") {
  # # test input section -----------------------------------------------------------
  # master.ref = fread('/juno/work/bergerm1/bergerlab/zhengy1/access_data_analysis/data/example_master_file.csv')
  # results.dir = paste0('/juno/work/bergerm1/MSK-ACCESS/ACCESS-Projects/test_access/access_data_analysis/output_042020/')
  # # criteria <- 'permissive'
  # criteria <- 'stringent'
  #
  # graph by patient --------------------------------------------------------
  output.dir <- paste0(results.dir, "/plots/")
  dir.create(output.dir)
  # for plotting consistency
  status_id <- c(
    "Called" = 19, "Not Called" = 4, "Signed out" = 15,
    "Not Signed out" = 13, "Not Covered" = 8, "Genotyped" = 17
  )

  # snv_sv_table = list.files(paste0(results.dir,'/results_',criteria,'_combined/'),full.names = T)
  lapply(unique(master.ref$cmo_patient_id), function(x) {
    # THIS PLOTS PLASMA SAMPLES ONLY
    # SNV
    tmp.table <- fread(list.files(paste0(results.dir, "/results_", criteria, "_combined/"), x, full.names = T))[
      call_confidence == "High" | grepl("Protein Fusion: in frame", HGVSp_Short)
    ]
    tmp.sample.sheets <- fread(paste0(results.dir, "/", x, "/", x, "_sample_sheet.tsv"))[, .(Sample_Barcode, cmo_patient_id, Sample_Type)]
    tmp.table <- table_to_maf(tmp.table, tmp.sample.sheets)
    tmp.table <- data.table(process_maf_for_graph(tmp.table))

    # CNA
    tmp.cna <- do.call(rbind, lapply(master.ref[cmo_patient_id == x]$cmo_sample_id_plasma, function(y) {
      fread(paste0(results.dir, "/CNA_final_call_set/", y, "_cna_final_call_set.txt"))
    }))

    # transform sample IDs into times
    if (all(!is.na(as.Date(master.ref[cmo_patient_id == x]$collection_date, "%m/%d/%y")))) {
      transform.vector <- structure(as.Date(master.ref[cmo_patient_id == x]$collection_date, "%m/%d/%y"),
        names = master.ref[cmo_patient_id == x]$cmo_sample_id_plasma
      )
      print(transform.vector)
    } else {
      transform.vector <- structure(as.character(master.ref[cmo_patient_id == x]$collection_date),
        names = master.ref[cmo_patient_id == x]$cmo_sample_id_plasma
      )
      print(transform.vector)
    }
    tmp.table$Tumor_Sample_Barcode <- transform.vector[tmp.table$Tumor_Sample_Barcode]
    factor.levels <- sort(unique(tmp.table$Tumor_Sample_Barcode))
    print(factor.levels)
    # tmp.table$Tumor_Sample_Barcode = factor(as.character(tmp.table$Tumor_Sample_Barcode),levels = factor.levels)
    tmp.table$Tumor_Sample_Barcode = as.character(tmp.table$Tumor_Sample_Barcode,format = "%Y-%b-%d")


    if (nrow(tmp.table) == 0 | all(tmp.table$t_alt_count == 0)) {
      print("skiping to the next")
      if (nrow(tmp.cna)) stop(paste0("Need to make CNA only file for: ", x))
      return()
    }

    colourCount <- nrow(unique(tmp.table[, .(Hugo_Symbol, HGVSp_Short)]))
    getPalette <- colorRampPalette(brewer.pal(8, "Set2"))
    SNV.SV.plot <- ggplot(tmp.table) +
      geom_line(aes(
        x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count / t_total_count)),
        color = paste0(Hugo_Symbol, " ", ifelse(grepl("^p\\.", HGVSp_Short), HGVSp_Short, "")), group = paste0(Hugo_Symbol, "_", HGVSp_Short)
      )) +
      geom_point(aes(
        x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count / t_total_count)),
        color = paste0(Hugo_Symbol, " ", ifelse(grepl("^p\\.", HGVSp_Short), HGVSp_Short, "")), shape = call_confidence
      ), size = 1.5) +
      labs(title = x, x = "Time Point", y = "log10(VAF)") +
      scale_x_discrete(breaks = sort(unique(tmp.table$Tumor_Sample_Barocde)),labels = sort(unique(tmp.table$Tumor_Sample_Barocde))) +
      #scale_x_date(date_labels = "%Y %b %d", breaks = "1 month") +
      scale_shape_manual(values = status_id, name = "Call Status") +
      scale_color_manual(values = getPalette(colourCount), name = "Alteration") +
      theme_minimal() +
      scale_y_log10() +
      theme(
        panel.grid.major = element_blank(), legend.position = "top", legend.box = "vertical",
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")
      )
    print(SNV.SV.plot)

    if (nrow(tmp.cna) > 0) {
      tmp.cna <- tmp.cna %>%
        mutate(Tumor_Sample_Barcode = factor(Tumor_Sample_Barcode, unique(tmp.sample.sheets[Sample_Type == "duplex"]$Sample_Barcode))) %>%
        # expand table on all empty samples without any calls
        data.table() %>%
        dcast.data.table(Hugo_Symbol + CNA ~ Tumor_Sample_Barcode, drop = c(TRUE, FALSE), fill = 0, value.var = "fc") %>%
        melt.data.table(id.vars = c("Hugo_Symbol", "CNA"), variable.name = "Tumor_Sample_Barcode", value.name = "fc") %>%
        data.table()
      tmp.table$Tumor_Sample_Barcode <- transform.vector[tmp.table$Tumor_Sample_Barcode]
      # factor.levels = sort(unique(tmp.table$Tumor_Sample_Barcode))
      # tmp.table$Tumor_Sample_Barcode = factor(as.character(tmp.table$Tumor_Sample_Barcode),levels = factor.levels)
      tmp.table$Tumor_Sample_Barcode = as.character(tmp.table$Tumor_Sample_Barcode, format = "%Y-%b-%d")
      colourCount <- nrow(unique(tmp.cna[, .(Hugo_Symbol, CNA)]))
      getPalette <- colorRampPalette(brewer.pal(8, "Set2"))
      CNA.plot <- ggplot(tmp.cna) +
        geom_bar(aes(x = Tumor_Sample_Barcode, y = abs(fc), fill = paste0(Hugo_Symbol, "_", CNA)), position = "dodge", stat = "identity") +
        labs(x = "Time Point", y = "Absolute fc") +
        scale_x_discrete(breaks = sort(unique(tmp.table$Tumor_Sample_Barocde)),labels = sort(unique(tmp.table$Tumor_Sample_Barocde))) +
        #scale_x_date(date_labels = "%Y %b %d", breaks = "1 month") +
        scale_fill_manual(values = getPalette(colourCount), name = "Alteration") +
        theme_minimal() +
        theme(panel.grid.major = element_blank(), legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))
      print(CNA.plot)

      pdf(paste0(output.dir, "/", x, "_all_events.pdf"), width = 10, height = 7)
      print(ggarrange(SNV.SV.plot, CNA.plot, ncol = 1, heights = c(2, 1)))
      dev.off()
    } else {
      pdf(paste0(output.dir, "/", x, "_all_events.pdf"), width = 10, height = 7)
      print(SNV.SV.plot)
      dev.off()
    }
  })
}

# Executable -----------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(tidyr)
  library(stringr)
  library(dplyr)
  library(argparse)
  library(ggplot2)
  library(ggpubr)
  library(RColorBrewer)
})

if (!interactive()) {
  parser <- ArgumentParser()
  parser$add_argument("-m", "--masterref", type = "character", help = "File path to master reference file")
  parser$add_argument("-o", "--resultsdir", type = "character", help = "Output directory")
  parser$add_argument("-c", "--criteria",
    type = "character", default = "stringent",
    help = "Calling criteria [default]"
  )
  args <- parser$parse_args()

  master.ref <- args$masterref
  results.dir <- args$resultsdir
  criteria <- args$criteria

  if (!criteria %in% c("stringent", "permissive")) {
    stop("Criteria argument should be either stringent or permissive")
  }

  suppressWarnings(plot_all_events(fread(master.ref), results.dir, criteria))
}