# VAF Overview Plot
# Author: Carmelina Charalambous, Alexander Ham
# Date: 11/30/2023

# load libraries
library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)
library(sqldf)
library(RSQLite)
library(readr)
library(argparse)
library(plotly)
library(htmlwidgets)
library(purrr) 


# Add arguments

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description = "Script creates VAF overview plots")

parser$add_argument("-o", "--resultsdir", type = "character", help = "Output directory")
parser$add_argument("-v", "--variants", type = "character", help = "File path to maf file")
parser$add_argument("-c", "--clinical", type = "character", help = "File path to clinical file")
parser$add_argument("-y", "--yaxis", type = "character", help = "Y-axis metric (mean, max, relative)", default = "mean")
parser$add_argument("-n", "--num_patients", type = "integer", default = 10, help = "Number of patients to include in each plot")
args = parser$parse_args()

# Access the arguments
results.dir <- args$resultsdir
variants.file <- args$variants
clinical.file <- args$clinical
yaxis.metric <- args$yaxis
num_patients <- args$num_patients

# Display inputs
cat("Results Directory:", results.dir, "\n")
cat("Variants File:", variants.file, "\n")
cat("Clinical File:", clinical.file, "\n")
cat("Y-axis Metric:", yaxis.metric, "\n")
cat("Number of Patients to display:", num_patients, "\n")

# load files
clinical_data <- read_tsv(clinical.file)
variants_data <- read_tsv(variants.file)


# test data
#clinical_data <- read_tsv("/Users/kcharanambou/Downloads/clinical_file_mock.txt")
#variants_data <- read_tsv("/Users/kcharanambou/Downloads/test_data.maf.txt")

# Check if Hugo_Symbol and HGVSp_Short columns exist
if("Hugo_Symbol" %in% colnames(variants_data) & "HGVSp_Short" %in% colnames(variants_data)) {
  variants_data$gene_variant <- paste(variants_data$Hugo_Symbol, variants_data$HGVSp_Short)
} else {
  stop("Required columns 'Hugo_Symbol' and/or 'HGVSp_Short' are missing in variants_data.")
}

# Concat gene-variant and categorize assay type
#variants_data$gene_variant = paste(variants_data$Hugo_Symbol, variants_data$HGVSp_Short)
variants_data$assay <- ifelse(grepl("^C-", variants_data$Tumor_Sample_Barcode), "ACCESS",
                    ifelse(grepl("*IM*", variants_data$Tumor_Sample_Barcode), "IMPACT", NA))


as.data.frame(variants_data)
# If covered column exists filter for covered if not continue with all ACCESS variants
if ("covered" %in% colnames(variants_data)) {
  access.variants <- variants_data %>% 
    filter (covered == "yes") %>%
    filter (assay == "ACCESS") 
} else {
  access.variants <- variants_data[variants_data$assay == "ACCESS", ]
}

#Create combined file

# combine
if("cmoSampleName" %in% colnames(clinical_data) && "Tumor_Sample_Barcode" %in% colnames(access.variants)) {
  variant.clinical.combined <- merge(clinical_data, access.variants, by.x = "cmoSampleName", by.y = "Tumor_Sample_Barcode")
} else {
  stop("Columns for merging not found or not unique in one of the data frames")
}

#clinical.combined <- merge(manifest_data, clinical_data, by.x = "cmoSampleName" , by.y = "cmo_sample_id_plasma")
#variant.clinical.combined <- merge(clinical_data, access.variants, by.x = "cmoSampleName" , by = "Tumor_Sample_Barcode")

keep.necessary <- variant.clinical.combined[, c("cmoSampleName", "cmoPatientId", "PatientId", "collection_date", "collection_in_days", "timepoint", "treatment_length", "treatmentName", "reason_for_tx_stop", "t_alt_freq", "DMP", "gene_variant")]



#get average vaf for each patient
keep.necessary$reason_for_tx_stop[is.na(keep.necessary$reason_for_tx_stop)] <- "NA"
vaf.av <- aggregate(t_alt_freq ~ cmoSampleName + cmoPatientId + collection_in_days + PatientId + treatment_length + reason_for_tx_stop, data = keep.necessary, FUN = mean)
vaf.av$reason_for_tx_stop[vaf.av$reason_for_tx_stop == "NA"] <- NA
#get initial and relative vaf values
vaf.av <- vaf.av %>% group_by(cmoPatientId) %>% mutate(init_vaf = t_alt_freq[which.min(collection_in_days)])
vaf.av$relative_vaf <- vaf.av$t_alt_freq / vaf.av$init_vaf


# Calculating the required statistics
vaf_stats <- keep.necessary %>%
  group_by(cmoSampleName) %>%
  summarise(
    AverageVAF = mean(t_alt_freq, na.rm = TRUE),
    MinVAF = min(t_alt_freq[t_alt_freq > 0], na.rm = TRUE),  # Exclude zeros in MinVAF calculation
    SDVAF = sd(t_alt_freq, na.rm = TRUE),
    MaxVAF = max(t_alt_freq, na.rm = TRUE)
  )

vaf_statistics <- merge(vaf.av, vaf_stats, by = "cmoSampleName", all.x = FALSE)
vaf_statistics <- vaf_stats2 %>% select(-t_alt_freq)


# Initialize yaxis.stats
yaxis.stats <- data.frame()

# Check for presence of t_alt_freq in keep.necessary
if (!"t_alt_freq" %in% names(keep.necessary)) {
  stop("Column 't_alt_freq' not found in keep.necessary")
}

# Calculate the Y-axis metric and merge with keep.necessary
if (yaxis.metric == "mean") {
  yaxis.stats <- aggregate(t_alt_freq ~ cmoSampleName, data = keep.necessary, FUN = mean)
  yaxis.label <- "Mean VAF"
} else if (yaxis.metric == "max") {
  yaxis.stats <- aggregate(t_alt_freq ~ cmoSampleName, data = keep.necessary, FUN = max)
  yaxis.label <- "Max VAF"
} else if (yaxis.metric == "relative") {
  yaxis.stats <- keep.necessary %>%
    group_by(cmoSampleName) %>%
    summarize(relative_vaf = mean(t_alt_freq / init_vaf))  # ensure aggregation to avoid dimension mismatch
  yaxis.label <- "Relative VAF"
} else {
  stop("Invalid Y-axis metric selected.")
}
# exception if 0 make it. 0.00001

# Merge the yaxis.stats with keep.necessary
yaxis.data <- merge(keep.necessary, yaxis.stats, by = "cmoSampleName", all.x = FALSE)
yaxis.data2 <- merge(yaxis.data, vaf.av, by = "cmoSampleName", all.x = FALSE)

# Function to split data into chunks based on the number of patients
split_data <- function(data, num) {
  if (nrow(data) == 0) {
    stop("Data frame is empty")
  }
  if (num <= 0) {
    stop("Number of patients per plot must be greater than 0")
  }
  unique_patients <- unique(data$cmoPatientId.x)
  if (length(unique_patients) < num) {
    stop("Number of patients per plot exceeds the number of unique patients")
  }
  # Split unique patients into chunks
  patient_chunks <- split(unique_patients, ceiling(seq_along(unique_patients) / num))
  # Return list of data frames, each containing all data for a chunk of patients
  lapply(patient_chunks, function(chunk) data[data$cmoPatientId.x %in% chunk, ])
}
# Splitting data for plots
patient_chunks <- split_data(yaxis.data2, num_patients)

# Generate and save plots for each chunk
# Generate and save plots for each chunk
for (i in seq_along(patient_chunks)) {
  chunk <- patient_chunks[[i]]
  # Create the necessary plots for this chunk
  # Replace the following placeholders with your actual plotting code
  init.vaf.plot <- ggplot(data = unique(chunk[c("cmoPatientId.x", "init_vaf")]), aes(x = init_vaf, y = NA)) + 
    geom_col(fill = "cadetblue1") + theme_classic() + xlab("VAF") + 
    geom_text(aes(x = init_vaf * 0.7, label = round(init_vaf, digits = 3)), fontface = "bold") + 
    facet_grid(rows = vars(cmoPatientId.x)) + ylab(element_blank()) + 
    scale_y_discrete(position = "right") + scale_x_reverse() + labs(title = "Initial VAF") +
    theme(panel.spacing = unit(0.5, "lines"), strip.text = element_blank(), axis.text.y = element_blank(), 
          axis.ticks = element_blank())
  vaf.grid <- ggplot(data = chunk, aes(x = collection_in_days.x, y = t_alt_freq.y, group = cmoPatientId.x)) + 
    geom_line(color = "darkslategray", size = 0.75) + geom_point(color = "darkslategray", size = 1.2) +
    geom_point(data = subset(chunk, t_alt_freq.y == 0), color = "white", size = 0.5) + theme_classic() + 
    theme(legend.position = "none") + 
    xlab("Treatment Day") + ylab(yaxis.label) +
    scale_y_continuous(limits = c(0, NA)) + facet_grid(rows = vars(cmoPatientId.x), scales = "free_y") + 
    theme(panel.spacing = unit(0.5, "lines"))
  treatment.length.plot <- ggplot(data = unique(chunk[c("cmoPatientId.x", "treatment_length.x")]), 
                                  aes(x = treatment_length.x, y = NA)) + geom_col(fill = "gold3") + theme_classic() + xlab("Days") + 
    geom_text(aes(x = treatment_length.x * 0.8, label = treatment_length.x, fontface = "bold")) + 
    facet_grid(rows = vars(cmoPatientId.x)) + ylab(element_blank()) + labs(title = "Treatment Duration") +
    theme(strip.text = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())
  tx.stop.plot <- ggplot(data = chunk, aes(y = 1, x = 0)) + 
    geom_tile(aes(fill = reason_for_tx_stop.x)) + scale_fill_brewer(palette = "Blues", direction = -1) +
    geom_text(aes(label = reason_for_tx_stop.x), size = 3.5) + theme_classic() + 
    labs(x = "", y = NULL, title = "Reason for Stopping Tx") + 
    facet_grid(rows = vars(cmoPatientId.x), labeller = labeller(cmoPatientId.x = dmp.cmo.key)) +
    theme(legend.position = "none", axis.text.y = element_blank(), axis.text.x = element_text(color = "white"), 
          axis.ticks = element_blank(), strip.text.y = element_text(size = 7, face = "bold", angle = 0), 
          plot.title = element_text(size = 12))

  
  # Combine the plots into a grid
  combined_plot <- grid.arrange(init.vaf.plot, vaf.grid, treatment.length.plot, tx.stop.plot, 
                                ncol = 4, widths = c(0.2, 0.45, 0.2, 0.15))
  
  # Save the combined plot
  plot_filename <- sprintf("VAF_overview_chunk_%d.pdf", i)
  ggsave(file.path(results.dir, plot_filename), combined_plot, width = 20, height = 10)

  # Convert to plotly and save as HTML
  combined_plot_html <- plotly::subplot(
    plotly::ggplotly(init.vaf.plot), 
    plotly::ggplotly(vaf.grid), 
    plotly::ggplotly(treatment.length.plot), 
    plotly::ggplotly(tx.stop.plot), 
    nrows = 1, margin = 0.01, widths = c(rep(1/4, 4))
  )

  plot_filename_html <- sprintf("VAF_overview_chunk_%d.html", i)
  htmlwidgets::saveWidget(combined_plot_html, file.path(results.dir, plot_filename_html), selfcontained = TRUE)
}

print("variant overview plot has been created in pdf and html format")

# save vaf average table
write.table(vaf_statistics, file = paste0(results.dir, "/vaf_statistics.txt"), sep = "\t", row.names = FALSE)
print("VAF statistics table saved")
