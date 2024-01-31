# Load necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(lubridate)
library(argparse)

# Create an Argument Parser
parser <- ArgumentParser(description = "Script creates swimmer plots")

parser$add_argument("-m", "--metadata", type = "character", help = "File path to metadata file")
parser$add_argument("-o", "--resultsdir", type = "character", help = "Output directory")
parser$add_argument("-c", "--colors", type = "character", help = "Comma-separated colors for treatment types", default = "blue,red,green,yellow")
parser$add_argument("-t", "--timeunit", type = "character", help = "Time unit for x-axis (days, weeks, months, years)", default = "days")
args <- parser$parse_args()

# Access the arguments
metadata.file <- args$metadata
results.dir <- args$resultsdir
time.unit <- args$timeunit
treatment_colors <- unlist(strsplit(args$colors, ","))

# Display inputs
cat("Metadata File:", metadata.file, "\n")
cat("Results Directory:", results.dir, "\n")
cat("Time Unit:", time.unit, "\n")

# Function to calculate time difference based on the specified unit
time_diff <- function(date1, date2, unit) {
  diff <- ifelse(unit == "years", 
                 as.numeric(interval(date2, date1) / years(1)), 
                 as.numeric(difftime(date1, date2, units = unit)))
  return(ifelse(is.na(diff), 0, diff)) # Return 0 for NA values
}

# Read the data from the file
data <- read_excel(metadata.file)

# Convert dates and calculate time differences
data <- data %>%
  mutate(
    start = as.Date(start),
    end = as.Date(end),
    collection_date = as.Date(collection_date),
    start_numeric = time_diff(start, min(start, na.rm = TRUE), time.unit),
    end_numeric = time_diff(end, min(start, na.rm = TRUE), time.unit)
  )

# Ensure that there are enough colors for the treatments
unique_treatments <- unique(data$treatment)
if(length(treatment_colors) < length(unique_treatments)) {
  stop("Not enough colors provided for all treatment types.")
}
color_mapping <- setNames(treatment_colors, unique_treatments)

# Create data frames for segments and points
data_segments <- data %>%
  select(dmp_patient_id, start_numeric, end_numeric, treatment) %>%
  arrange(dmp_patient_id, start_numeric) %>%
  filter(!is.na(treatment))

data_lines <- data_segments %>%
  group_by(dmp_patient_id) %>%
  mutate(next_start = lead(start_numeric)) %>%
  filter(!is.na(next_start))

data_points <- data %>%
  filter(ctdna_detection %in% c('positive', 'negative')) %>%
  select(dmp_patient_id, collection_date, ctdna_detection) %>%
  mutate(collection_date_numeric = time_diff(collection_date, min(data$start, na.rm = TRUE), time.unit))

# Create the swimmer plot
swimmer_plot <- ggplot() +
  geom_segment(data = data_lines, aes(x = start_numeric, xend = next_start, y = dmp_patient_id, yend = dmp_patient_id), linewidth = 1, color = "grey") +
  geom_segment(data = data_segments, aes(x = start_numeric, xend = end_numeric, y = dmp_patient_id, yend = dmp_patient_id, color = treatment), linewidth = 4) +
  geom_point(data = data_points, aes(x = collection_date_numeric, y = dmp_patient_id, shape = ctdna_detection), size = 1) +
  labs(title = paste("Swimmer Plot (", time.unit, ")", sep=""),
       x = paste("Time (", time.unit, ")", sep=""),
       y = "Patient ID",
       color = "Treatment Type",
       shape = "ctDNA Detection") +
  theme_minimal() +
  scale_color_manual(values = color_mapping) +
  scale_shape_manual(values = c(16, 17))

# Print or save the plot
plot_filename <- paste0(results.dir, "/swimmer_plot_", time.unit, ".png")
ggsave(plot_filename, swimmer_plot, width = 10, height = 6)
cat("Swimmer plot saved to:", plot_filename, "\n")
