# Swimmer Plot Scripts

## Overview

The `swimmer` folder contains R scripts designed to create swimmer plots for visualizing treatment timelines and related data. These scripts process input data, calculate time differences, and generate swimmer plots for single and multiple treatments. The plots are saved as PDF or PNG files for further analysis and reporting.

## Scripts

### 1. `swimmer_single_treatment.R`

#### Description
This script generates swimmer plots for single-treatment data. It processes input data, calculates time differences, and creates a swimmer plot with various visualizations, including treatment timelines and assay types.

#### Features
- Processes input data to calculate time differences.
- Generates swimmer plots for single-treatment data.
- Supports multiple time units (days, weeks, months, years).
- Saves the plot as a PDF file.

#### Arguments
| Argument       | Type       | Description                                           | Default Value |
|-----------------|------------|-------------------------------------------------------|---------------|
| `-i, --input`   | `character` | File path to the input data file.                     | None          |
| `-o, --output`  | `character` | File path for the output PDF file.                    | None          |
| `-t, --timeunit`| `character` | Time unit for the x-axis (days, weeks, months, years).| `days`        |

#### Example Command
```bash
Rscript swimmer_single_treatment.R -i input_data.txt -o output_plot.pdf -t days
```

---

### 2. `swimmer_multi_treatment.R`

#### Description
This script generates swimmer plots for multi-treatment data. It processes metadata, calculates time differences, and creates a swimmer plot with treatment timelines and ctDNA detection points.

#### Features
- Processes metadata to calculate time differences.
- Generates swimmer plots for multi-treatment data.
- Supports multiple time units (days, weeks, months, years).
- Allows customization of treatment colors.
- Saves the plot as a PNG file.

#### Arguments
| Argument       | Type       | Description                                           | Default Value |
|-----------------|------------|-------------------------------------------------------|---------------|
| `-m, --metadata`| `character` | File path to the metadata file.                       | None          |
| `-o, --resultsdir`| `character` | Output directory for the plot.                       | None          |
| `-c, --colors`  | `character` | Comma-separated colors for treatment types.           | `blue,red,green,yellow` |
| `-t, --timeunit`| `character` | Time unit for the x-axis (days, weeks, months, years).| `days`        |

#### Example Command
```bash
Rscript swimmer_multi_treatment.R -m metadata.xlsx -o /path/to/output -c blue,red,green -t weeks
```

---

### 3. `dates2days.R`

#### Description
This script converts date columns in the input data to numeric values representing time differences in specified units. The processed data is saved as a tab-delimited text file for use in swimmer plots.

#### Features
- Converts date columns to numeric time differences.
- Supports multiple time units (days, weeks, months, years).
- Saves the processed data as a tab-delimited text file.

#### Arguments
| Argument       | Type       | Description                                           | Default Value |
|-----------------|------------|-------------------------------------------------------|---------------|
| `-i, --input`   | `character` | File path to the input `.txt` file.                   | None          |
| `-o, --output`  | `character` | File path for the output `.txt` file.                 | None          |

#### Example Command
```bash
Rscript dates2days.R -i input_data.txt -o output_data.txt
```

---

## Requirements

### R Packages
The scripts require the following R packages:
- `dplyr`
- `ggplot2`
- `lubridate`
- `argparse`
- `readr`
- `readxl`
- `tidyr`
- `scales`
- `gridExtra`
- `cowplot`

Install the required packages using the following command:
```R
install.packages(c("dplyr", "ggplot2", "lubridate", "argparse", "readr", "readxl", "tidyr", "scales", "gridExtra", "cowplot"))
```

---

## Input File Requirements

### Single Treatment Input File
The input file for `swimmer_single_treatment.R` must contain the following columns:
- `collection_date`
- `start`
- `endtouse`
- `reason`
- `assay_type`
- `clinical_or_research`

### Multi-Treatment Metadata File
The metadata file for `swimmer_multi_treatment.R` must contain the following columns:
- `start`
- `end`
- `collection_date`
- `treatment`
- `ctdna_detection`

### Dates to Days Input File
The input file for `dates2days.R` must contain date columns such as:
- `pre_tx_date`
- `start`
- `end`

---

## Outputs

### Swimmer Plots
- **Single Treatment**: PDF file containing the swimmer plot.
- **Multi-Treatment**: PNG file containing the swimmer plot.

### Processed Data
- Tab-delimited text file with numeric time differences for use in swimmer plots.

---

## Example Workflow

1. **Convert Dates to Days**:
   ```bash
   Rscript dates2days.R -i input_data.txt -o processed_data.txt
   ```

2. **Generate Single Treatment Swimmer Plot**:
   ```bash
   Rscript swimmer_single_treatment.R -i processed_data.txt -o single_treatment_plot.pdf -t days
   ```

3. **Generate Multi-Treatment Swimmer Plot**:
   ```bash
   Rscript swimmer_multi_treatment.R -m metadata.xlsx -o /path/to/output -c blue,red,green -t weeks
   ```

---

## Contact

For questions or issues, please contact:

- **Author**: Carmelina Charalambous, Alexander Ham
- **Date**: 11/30/2023