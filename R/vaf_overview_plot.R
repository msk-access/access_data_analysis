# VAF Overview Plot
# Author: Carmelina Charalambous, Alexander Ham
# Date: 11/30/2023

# load libraries
library(magrittr)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)
library(sqldf)
library(RSQLite)
library(readr)
library(argparse)
library(data.table)

if (!requireNamespace("argparse", quietly = TRUE)) {
  install.packages("argparse")
}

# Add arguments

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(description = "Script creates VAF overview plots")

parser$add_argument("-m", "--manifest", type = "character", help = "File path to manifest file")
parser$add_argument("-o", "--resultsdir", type = "character", help = "Output directory")
parser$add_argument("-v", "--variants", type = "character", help = "File path to maf file")
parser$add_argument("-c", "--clinical", type = "character", help = "File path to clinical file")

args = parser$parse_args()

# Access the arguments
manifest.file <- args$manifest
results.dir <- args$resultsdir
variants.file <- args$variants
clinical.file <- args$clinical

# Display inputs
cat("Manifest File:", manifest.file, "\n")
cat("Results Directory:", results.dir, "\n")
cat("Variants File:", variants.file, "\n")
cat("Clinical File:", clinical.file, "\n")



# load files
manifest_data <- read.csv(manifest.file)
clinical_data <- read_tsv(clinical.file)
variants_data <- read_csv(variants.file)

# Concat gene-variant and categorize assay type
variants_data$gene_variant = paste(variants_data$Hugo_Symbol, variants_data$HGVSp_Short)
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
clinical.combined <- merge(manifest_data, clinical_data, by.x = "cmoSampleName" , by.y = "cmo_sample_id_plasma")
variant.clinical.combined <- merge(clinical.combined, access.variants, by.x = "cmoSampleName" , by.y = "Tumor_Sample_Barcode")
print(colnames(variant.clinical.combined))
keep.necessary <- variant.clinical.combined[, c("cmoSampleName", "cmoPatientId", "investigatorSampleId", "DMP_PATIENT_ID", "collection_date", "collection_in_days", "timepoint", "treatment_length", "treatmentName", "reason_for_tx_stop", "t_alt_freq", "DMP", "gene_variant")]
print(keep.necessary)

#get average vaf for each patient
keep.necessary$reason_for_tx_stop[is.na(keep.necessary$reason_for_tx_stop)] <- "NA"
vaf.av <- aggregate(t_alt_freq ~ cmoSampleName + cmoPatientId + collection_in_days + DMP_PATIENT_ID + treatment_length + reason_for_tx_stop, data = keep.necessary, FUN = mean)
vaf.av$reason_for_tx_stop[vaf.av$reason_for_tx_stop == "NA"] <- NA

#get initial and relative vaf values
vaf.av <- vaf.av %>% group_by(cmoPatientId) %>% mutate(init_vaf = t_alt_freq[which.min(collection_in_days)])
vaf.av$relative_vaf <- vaf.av$t_alt_freq / vaf.av$init_vaf

#order patients based on earliest collection date
vaf.av <- vaf.av[order(vaf.av$treatment_length, decreasing = TRUE),]
vaf.av$cmoPatientId <- factor(vaf.av$cmoPatientId, levels = unique(vaf.av$cmoPatientId))

#CMO and DMP IDs for labeling
dmp.cmo.key <- paste0(unique(vaf.av$cmoPatientId), "\n\n", unique(vaf.av$DMP_PATIENT_ID)) %>% as.vector()
names(dmp.cmo.key) <- unique(vaf.av$cmoPatientId) %>% as.vector()

#plots
init.vaf.plot <- ggplot(data = unique(vaf.av[c("cmoPatientId", "init_vaf")]), aes(x = init_vaf, y = NA)) + 
  geom_col(fill = "cadetblue1") + theme_classic() + xlab("VAF") + 
  geom_text(aes(x = init_vaf * 0.7, label = round(init_vaf, digits = 3)), nudge_x = -0.02, fontface = "bold") + 
  facet_grid(rows = vars(cmoPatientId)) + ylab(element_blank()) + 
  scale_y_discrete(position = "right") + scale_x_reverse() + labs(title = "Initial VAF") +
  theme(panel.spacing = unit(0.5, "lines"), strip.text = element_blank(), axis.text.y = element_blank(), 
  axis.ticks = element_blank())

vaf.grid <- ggplot(data = vaf.av, aes(x = collection_in_days, y = t_alt_freq, group = cmoPatientId)) + 
  geom_line(color = "darkslategray", size = 0.75) + geom_point(color = "darkslategray", size = 1.2) +
  geom_point(data = subset(vaf.av, t_alt_freq == 0), color = "white", size = 0.5) + theme_classic() + 
  theme(legend.position = "none", panel.background = element_rect(fill = "aliceblue", colour = "aliceblue")) + 
  xlab("Treatment Day") + ylab("") + labs(title = "Average VAF") + 
  scale_y_continuous(limits = c(0, NA)) + facet_grid(rows = vars(cmoPatientId), scales = "free_y") + 
  theme(panel.spacing = unit(0.5, "lines"), strip.text = element_blank())

vaf.grid.relative <- ggplot(data = vaf.av, aes(x = collection_in_days, y = relative_vaf, group = cmoPatientId)) + 
  geom_line(color = "darkslategray", size = 0.75) + geom_point(color = "darkslategray", size = 1.2) +
  geom_point(data = subset(vaf.av, t_alt_freq == 0), color = "white", size = 0.5) + theme_classic() + 
  theme(legend.position = "none", panel.background = element_rect(fill = "aliceblue", colour = "aliceblue")) + 
  xlab("Treatment Day") + ylab("") + labs(title = "Relative Average VAF") + 
  scale_y_continuous(limits = c(0, NA)) + facet_grid(rows = vars(cmoPatientId), scales = "free_y") + 
  theme(panel.spacing = unit(0.5, "lines"), strip.text = element_blank())

treatment.length.plot <- ggplot(data = unique(vaf.av[c("cmoPatientId", "treatment_length")]), 
  aes(x = treatment_length, y = NA)) + geom_col(fill = "gold3") + theme_classic() + xlab("Days") + 
  geom_text(aes(x = treatment_length * 0.8, label = treatment_length, fontface = "bold")) + 
  facet_grid(rows = vars(cmoPatientId)) + ylab(element_blank()) + labs(title = "Treatment Duration") +
  theme(strip.text = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())

tx.stop.plot <- ggplot(data = vaf.av, aes(y = 1, x = 0)) + 
  geom_tile(aes(fill = reason_for_tx_stop)) + scale_fill_brewer(palette = "Blues", direction = -1) +
  geom_text(aes(label = reason_for_tx_stop), size = 3.5) + theme_classic() + 
  labs(x = "", y = NULL, title = "Reason for Stopping Tx") + 
  facet_grid(rows = vars(cmoPatientId), labeller = labeller(cmoPatientId = dmp.cmo.key)) +
  theme(legend.position = "none", axis.text.y = element_blank(), axis.text.x = element_text(color = "white"), 
  axis.ticks = element_blank(), strip.text.y = element_text(size = 7, face = "bold", angle = 0), 
  plot.title = element_text(size = 12))



#pdf(results.dir, "VAF_overview_plot.pdf",width=14, height = 10)
plot1 <- grid.arrange(init.vaf.plot, vaf.grid, treatment.length.plot, tx.stop.plot, ncol = 4, widths = c(0.2, 0.45, 0.2, 0.15))
dev.off()

# Create the output PDF file path 1
output_pdf1 <- file.path(results.dir, "VAF_overview_plot.pdf")


# Save the plot as a PDF in the specified output directory 1
ggsave(output_pdf1, plot1)

#pdf(results.dir,"VAF_overview_plot_relative.pdf",width=14, height = 10)
plot2 <- grid.arrange(init.vaf.plot, vaf.grid.relative, treatment.length.plot, tx.stop.plot, ncol = 4, widths = c(0.2, 0.45, 0.2, 0.15))
dev.off()

# Create the output PDF file path 2
output_pdf2 <- file.path(results.dir, "VAF_overview_plot_relative.pdf")

# Save the plot as a PDF in the specified output directory 2
ggsave(output_pdf2, plot2)

print("variant plots have been created")