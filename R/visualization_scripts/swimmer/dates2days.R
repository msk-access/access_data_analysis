# Load necessary libraries
library(dplyr)
library(lubridate)
library(argparse)
library(readr) 

# Use to convert dates to days (input for swimmer_plots when running on juno

# Create an Argument Parser
parser <- ArgumentParser(description = "Process and save swimmer plot data")

# Define arguments
parser$add_argument("-i", "--input", required = TRUE, help = "File path to input .txt file")
parser$add_argument("-o", "--output", required = TRUE, help = "File path for the output .txt file")

# Parse the arguments
args <- parser$parse_args()

# Read the data from the input file (assuming a tabular format like CSV)
data <- read_delim(args$input, delim = "\t", show_col_types = FALSE) # Adjust delimiter if necessary

# Process the data
data_processed <- data %>%
  mutate(
    pre_tx_date = as.Date(pre_tx_date),
    start = as.Date(start),
    end = as.Date(end)
  ) %>%
  group_by(cmoPatientId) %>%
  mutate(
    start = as.numeric(difftime(start, pre_tx_date, units = "days")),
    end = as.numeric(difftime(end, pre_tx_date, units = "days")),
    collection_date = as.numeric(difftime(collection_date, pre_tx_date, units = "days")),
    endtouse = as.numeric(difftime(endtouse, pre_tx_date, units = "days")),
    pre_tx_date = as.numeric(difftime(pre_tx_date, pre_tx_date, units = "days")),
  ) %>%
  ungroup()

# Save the processed data to the output file 
write_delim(data_processed, args$output, delim = "\t") # Adjust delimiter if necessary

# Confirm saved file
cat("Processed data saved to:", args$output, "\n")
