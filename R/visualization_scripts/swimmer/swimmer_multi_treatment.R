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
args <- parser$parse_args()

# Access the arguments
metadata.file <- args$metadata
results.dir <- args$resultsdir
treatment_colors <- unlist(strsplit(args$colors, ","))

# Display inputs
cat("Metadata File:", metadata.file, "\n")
cat("Results Directory:", results.dir, "\n")

# Read the data from the file
data <- read_excel(metadata.file)

# Process data for plotting segments, excluding NA in treatment
data_segments <- data %>%
  select(dmp_patient_id, start, end, treatment) %>%
  mutate(start_numeric = as.numeric(difftime(start, min(start, na.rm = TRUE), units = "days")),
         end_numeric = as.numeric(difftime(end, min(start, na.rm = TRUE), units = "days"))) %>%
  arrange(dmp_patient_id, start) %>%
  filter(!is.na(treatment))  # Exclude NA values in treatment for plotting

# Ensure that there are enough colors for the treatments
unique_treatments <- unique(data_segments$treatment)
if(length(treatment_colors) < length(unique_treatments)) {
  stop("Not enough colors provided for all treatment types.")
}
color_mapping <- setNames(treatment_colors, unique_treatments)

# Create a dataframe for connecting lines
data_lines <- data_segments %>%
  group_by(dmp_patient_id) %>%
  arrange(dmp_patient_id, start_numeric) %>%
  mutate(next_start = lead(start_numeric)) %>%
  filter(!is.na(next_start)) %>%
  select(dmp_patient_id, end_numeric, next_start) %>%
  ungroup()

# Filter for ctDNA detection points
data_points <- data %>%
  filter(ctdna_detection %in% c('positive', 'negative')) %>%
  select(dmp_patient_id, collection_date, ctdna_detection) %>%
  mutate(collection_date_numeric = as.numeric(difftime(collection_date, min(data$start, na.rm = TRUE), units = "days")))

# Create the swimmer plot
swimmer_plot <- ggplot() +
  geom_segment(data = data_lines, aes(x = end_numeric, xend = next_start, y = dmp_patient_id, yend = dmp_patient_id), linewidth = 1, color = "grey") +
  geom_segment(data = data_segments, aes(x = start_numeric, xend = end_numeric, y = dmp_patient_id, yend = dmp_patient_id, color = treatment), linewidth = 4) +
  geom_point(data = data_points, aes(x = collection_date_numeric, y = dmp_patient_id, shape = ctdna_detection), size = 1) +
  labs(title = "Swimmer Plot (days)",
       x = "Days since first treatment",
       y = "Patient ID",
       color = "Treatment Type",
       shape = "ctDNA Detection") +
  theme_minimal() +
  scale_color_manual(values = color_mapping) +
  scale_shape_manual(values = c(16, 17))

# Print or save the plot
ggsave(paste0(results.dir, "/swimmer_plot_days.png"), swimmer_plot, width = 10, height = 6)


#Create second plot but days in months
# Function to calculate the difference in months
month_diff <- function(date1, date2) {
  interval <- as.period(interval(date1, date2), unit = "month")
  interval$month + interval$year * 12
}

# Process data for plotting segments, excluding NA in treatment
data_segments <- data %>%
  select(dmp_patient_id, start, end, treatment) %>%
  mutate(
    start_numeric = month_diff(min(start, na.rm = TRUE), start),
    end_numeric = month_diff(min(start, na.rm = TRUE), end)
  ) %>%
  arrange(dmp_patient_id, start) %>%
  filter(!is.na(treatment))  # Exclude NA values in treatment for plotting

# Ensure that there are enough colors for the treatments
unique_treatments <- unique(data_segments$treatment)
if(length(treatment_colors) < length(unique_treatments)) {
  stop("Not enough colors provided for all treatment types.")
}
color_mapping <- setNames(treatment_colors, unique_treatments)

# Create a dataframe for connecting lines
data_lines <- data_segments %>%
  group_by(dmp_patient_id) %>%
  arrange(dmp_patient_id, start_numeric) %>%
  mutate(next_start = lead(start_numeric)) %>%
  filter(!is.na(next_start)) %>%
  select(dmp_patient_id, end_numeric, next_start) %>%
  ungroup()

# Filter for ctDNA detection points
data_points <- data %>%
  filter(ctdna_detection %in% c('positive', 'negative')) %>%
  select(dmp_patient_id, collection_date, ctdna_detection) %>%
  mutate(collection_date_numeric = month_diff(min(data$start, na.rm = TRUE), collection_date))

# Create the swimmer plot
swimmer_plot <- ggplot() +
  geom_segment(data = data_lines, aes(x = end_numeric, xend = next_start, y = dmp_patient_id, yend = dmp_patient_id), linewidth = 1, color = "grey") +
  geom_segment(data = data_segments, aes(x = start_numeric, xend = end_numeric, y = dmp_patient_id, yend = dmp_patient_id, color = treatment), linewidth = 4) +
  geom_point(data = data_points, aes(x = collection_date_numeric, y = dmp_patient_id, shape = ctdna_detection), size = 1) + # Adjusted size here
  labs(title = "Swimmer Plot (months)",
       x = "Months since first treatment",
       y = "Patient ID",
       color = "Treatment Type",
       shape = "ctDNA Detection") +
  theme_minimal() +
  scale_color_manual(values = color_mapping) +
  scale_shape_manual(values = c(16, 17))

# Print or save the plot
ggsave(paste0(results.dir, "/swimmer_plot_months.png"), swimmer_plot, width = 10, height = 6)

print("Swimmer plots created")
