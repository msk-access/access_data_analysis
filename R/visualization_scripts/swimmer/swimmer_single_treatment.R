# Load necessary libraries
library(data.table)
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(scales)
library(gridExtra)
library(png)
library(grid)
library(cowplot)
library(devtools)
library(tidyverse)
library(argparse)

# Create an Argument Parser
parser <- ArgumentParser(description = "Create swimmer plot for single treatment")

# Define arguments
parser$add_argument("-i", "--input", required = TRUE, help = "File path to input data file")
parser$add_argument("-o", "--output", required = TRUE, help = "File path for the output PDF file")
parser$add_argument("-t", "--timeunit", type = "character", help = "Time unit for x-axis (days, weeks, months, years)", default = "days")

# Parse the arguments
args <- parser$parse_args()

# Function to calculate time difference based on the specified unit
time_diff <- function(date1, date2, unit) {
  unit <- match.arg(unit, c("days", "weeks", "months", "years"))
  diff <- case_when(
    unit == "years" ~ as.numeric(interval(date2, date1) / years(1)),
    unit == "months" ~ as.numeric(interval(date2, date1) / months(1)),
    unit == "weeks" ~ as.numeric(interval(date2, date1) / weeks(1)),
    TRUE ~ as.numeric(difftime(date1, date2, units = unit))
  )
  return(ifelse(is.na(diff), 0, diff)) # Return 0 for NA values
}

# Read input file
tab <- fread(args$input) %>%
  mutate(
    collection_date = as.Date(collection_date),
    start = as.Date(start),
    endtouse = as.Date(endtouse),
    collection_date_numeric = time_diff(collection_date, start, args$timeunit),
    endtouse_numeric = time_diff(endtouse, start, args$timeunit)
  ) %>%
  mutate(within100days = factor(ifelse(collection_date_numeric >= -100 & collection_date_numeric <= 100, "yes", "no"), levels = c("yes", "no"))) %>%
  mutate(patient_id = paste0(cmoPatientId, " (", dmp_patient_id, ")"))

theme_set(theme_classic())

# Ensure endtouse_numeric is available for plot3
tab_plot3 <- tab %>%
  select(patient_id, endtouse, start, reason, endtouse_numeric) %>%
  unique()

# Plot 0
plot0 <- tab %>% ggplot(aes(y = fct_reorder(patient_id, endtouse_numeric))) +
  geom_bar(aes(fill = factor(paste(assay_type, clinical_or_research), levels = c("ACCESS research", "ACCESS clinical", "IMPACT clinical")), alpha = within100days), position = position_stack(reverse = TRUE))  +
  scale_alpha_discrete(range = c(1, 0.3), name = "within 100 days of tx?") +
  scale_x_reverse() + 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank()) +
  scale_y_discrete(position = "right") + 
  theme(legend.position = "left") + 
  scale_fill_brewer(name = "", palette = "Accent") + xlab("# Samples")

# Plot 1
plot1 <- tab %>% filter(collection_date_numeric >= -100 & collection_date_numeric <= 100) %>% 
  ggplot(aes(y = fct_reorder(patient_id, endtouse_numeric), x = collection_date_numeric)) +
  geom_segment(aes(x = 0, xend = endtouse_numeric, yend = patient_id), linewidth = 1) + 
  geom_point(aes(shape = assay_type, col = factor(paste(assay_type, clinical_or_research), levels = c("ACCESS research", "ACCESS clinical", "IMPACT clinical"))), size = 3, stroke = 1) + 
  theme(legend.title = element_blank(), axis.title.y = element_blank()) + xlab(paste("Time from Treatment Start (", args$timeunit, ")", sep = "")) + 
  scale_shape_manual(values = c(1, 8)) + guides(shape = "none") + 
  theme(legend.position = "none") + 
  scale_color_brewer(palette = "Accent") + 
  theme(panel.grid.major.y = element_line(color = "lightgray", linetype = "dashed")) + 
  geom_segment(data = tab %>% filter(is.na(reason)), aes(x = endtouse_numeric, xend = endtouse_numeric + 1, yend = patient_id), arrow = arrow(length = unit(0.2, "cm")))

# Plot 2
plot2 <- tab %>% filter(collection_date_numeric >= -100 & collection_date_numeric <= 100) %>% 
  ggplot(aes(y = fct_reorder(patient_id, endtouse_numeric), x = collection_date)) + 
  geom_segment(aes(x = start, xend = endtouse, yend = patient_id), linewidth = 1) + 
  geom_point(aes(shape = assay_type, col = factor(paste(assay_type, clinical_or_research), levels = c("ACCESS research", "ACCESS clinical", "IMPACT clinical"))), size = 3, stroke = 1) + 
  theme(legend.title = element_blank(), axis.title.y = element_blank()) + 
  xlab("Date") + 
  scale_shape_manual(values = c(1, 8)) + 
  guides(shape = "none") + 
  theme(legend.position = "none") + 
  scale_color_brewer(palette = "Accent") + 
  theme(panel.grid.major.y = element_line(color = "lightgray", linetype = "dashed")) + 
  geom_segment(data = tab %>% filter(is.na(reason)), aes(x = endtouse, xend = endtouse + 1, yend = patient_id), arrow = arrow(length = unit(0.2, "cm")))

# Plot 3
plot3 <- tab_plot3 %>%
  ggplot(aes(y = fct_reorder(patient_id, endtouse_numeric), x = " ")) +
  geom_tile(aes(fill = reason), col = "white") +
  geom_text(aes(label = reason), size = 2) +
  xlab("Tx\nReason") +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_blank()) +
  scale_fill_brewer(palette = "Blues", direction = -1)

# Save the plot to a pdf file
ggsave(args$output, plot = plot_grid(plot0, plot1, plot3, nrow = 1, rel_widths = c(2, 4.2, 0.8), align = "tb"), width = 11, height = 5.5)
