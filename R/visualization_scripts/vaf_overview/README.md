# VAF Overview Plot Script

## Overview

This script, `vaf_overview_plot.R`, generates Variant Allele Frequency (VAF) overview plots for clinical and variant data. It creates visualizations in both PDF and HTML formats, providing insights into VAF trends, treatment durations, and reasons for stopping treatment for a specified number of patients.

## Features

- **Input Parsing**: Accepts clinical and variant data files as input.
- **Data Validation**: Ensures required columns are present in the input files.
- **Data Processing**:
  - Merges clinical and variant data.
  - Filters and categorizes data based on assay type.
  - Calculates VAF statistics (mean, max, relative VAF).
- **Visualization**:
  - Generates plots for initial VAF, VAF trends, treatment duration, and reasons for stopping treatment.
  - Combines plots into a grid for each patient chunk.
- **Output**:
  - Saves plots in both PDF and HTML formats.
  - Exports VAF statistics as a tab-delimited text file.

## Requirements

### R Packages

The script requires the following R packages:

- `ggplot2`
- `gridExtra`
- `tidyr`
- `dplyr`
- `sqldf`
- `RSQLite`
- `readr`
- `argparse`
- `plotly`
- `htmlwidgets`
- `purrr`

Install the required packages using the following command:

```R
install.packages(c("ggplot2", "gridExtra", "tidyr", "dplyr", "sqldf", "RSQLite", "readr", "argparse", "plotly", "htmlwidgets", "purrr"))
```

## Usage

### Command-Line Arguments

The script accepts the following arguments:

| Argument         | Type       | Description                                                                 | Default Value |
|-------------------|------------|-----------------------------------------------------------------------------|---------------|
| `-o, --resultsdir` | `character` | Output directory where plots and statistics will be saved.                  | None          |
| `-v, --variants`  | `character` | File path to the variant data (MAF file).                                   | None          |
| `-c, --clinical`  | `character` | File path to the clinical data file.                                        | None          |
| `-y, --yaxis`     | `character` | Y-axis metric for VAF plots (`mean`, `max`, or `relative`).                  | `mean`        |
| `-n, --num_patients` | `integer`  | Number of patients to include in each plot.                                 | `10`          |

### Example Command

```bash
Rscript vaf_overview_plot.R -o /path/to/output -v /path/to/variants.maf -c /path/to/clinical.tsv -y mean -n 10
```

### Input File Requirements

#### Clinical Data File

The clinical data file must be a tab-delimited file containing the following columns:

- `cmoSampleName`
- `cmoPatientId`
- `PatientId`
- `collection_date`
- `collection_in_days`
- `timepoint`
- `treatment_length`
- `treatmentName`
- `reason_for_tx_stop`

#### Variant Data File

The variant data file must be a tab-delimited file containing the following columns:

- `Hugo_Symbol`
- `HGVSp_Short`
- `Tumor_Sample_Barcode`
- `t_alt_freq`
- `covered` (optional)

## Outputs

1. **Plots**:
   - PDF files: One file per patient chunk (e.g., `VAF_overview_chunk_1.pdf`).
   - HTML files: Interactive plots for each patient chunk (e.g., `VAF_overview_chunk_1.html`).

2. **Statistics**:
   - A tab-delimited text file (`vaf_statistics.txt`) containing VAF statistics for all patients.

## Script Workflow

1. **Input Parsing**:
   - Reads the clinical and variant data files.
   - Validates the presence of required columns.

2. **Data Processing**:
   - Merges clinical and variant data.
   - Filters and categorizes variants based on assay type.
   - Calculates VAF statistics (mean, max, relative VAF).

3. **Visualization**:
   - Splits data into chunks based on the number of patients specified.
   - Generates the following plots for each chunk:
     - Initial VAF
     - VAF trends over time
     - Treatment duration
     - Reasons for stopping treatment
   - Combines the plots into a grid and saves them as PDF and HTML files.

4. **Output**:
   - Saves the combined plots and VAF statistics.

## Error Handling

The script includes error handling for the following scenarios:

- Missing required columns in the input files.
- Empty data frames after filtering.
- Invalid Y-axis metric.
- Number of patients per plot exceeding the total number of unique patients.

## Example Outputs

### PDF Plot

The PDF plot contains the following panels for each patient:

1. **Initial VAF**: Bar plot showing the initial VAF.
2. **VAF Trends**: Line plot showing VAF trends over time.
3. **Treatment Duration**: Bar plot showing the treatment duration in days.
4. **Reason for Stopping Treatment**: Tile plot showing the reason for stopping treatment.

### HTML Plot

The HTML plot is an interactive version of the PDF plot, allowing users to explore the data dynamically.

### VAF Statistics

The `vaf_statistics.txt` file contains the following columns:

- `cmoSampleName`
- `cmoPatientId`
- `collection_in_days`
- `PatientId`
- `treatment_length`
- `reason_for_tx_stop`
- `AverageVAF`
- `MinVAF`
- `SDVAF`
- `MaxVAF`

## Contact

For questions or issues, please contact:

- **Author**: Carmelina Charalambous, Alexander Ham
- **Date**: 11/30/2023
