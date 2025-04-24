# Manifest Update Script

## Overview

This Python script processes and updates an ACCESS manifest file by generating paths for various data types (e.g., BAM, MAF, CNA, SV files) and saves the updated manifest in both Excel and CSV formats. It supports both legacy and modern input formats and includes options for handling Protected Health Information (PHI).

## Features

- **Input Validation**:
  - Ensures required columns are present in the input manifest.
  - Validates date formats and handles missing values.
- **Path Generation**:
  - Automatically generates paths for BAM, MAF, CNA, and SV files based on sample type and assay type.
- **PHI Handling**:
  - Optionally removes collection dates to comply with privacy regulations.
- **Output**:
  - Saves the updated manifest in both Excel and CSV formats.
  - Supports custom output file prefixes.
- **Legacy Support**:
  - Handles legacy input file formats with specific path requirements.

## Requirements

### Python Packages

The script requires the following Python packages:

- `pandas`
- `typer`
- `rich`
- `arrow`
- `numpy`
- `openpyxl` (for Excel file handling)

Install the required packages using the following command:

```bash
pip install pandas typer rich arrow numpy openpyxl
```

## Usage

### Commands

The script provides two main commands:

1. **`make-manifest`**:
   Processes the input manifest file to generate paths for various data types and saves the updated manifest.

2. **`update-manifest`**:
   Updates a legacy ACCESS manifest file with specific paths.

### Command-Line Arguments

#### `make-manifest`

| Argument                  | Type       | Description                                                                 | Default Value |
|---------------------------|------------|-----------------------------------------------------------------------------|---------------|
| `-i, --input`             | `Path`     | Path to the input manifest file.                                            | None          |
| `-o, --output`            | `str`      | Prefix name for the output files (without extension).                       | None          |
| `--remove-collection-date`| `bool`     | Remove collection date from the output manifest (PHI).                      | `False`       |
| `-a, --assay-type`        | `str`      | Assay type, either `XS1` or `XS2`.                                          | `XS2`         |

#### `update-manifest`

| Argument      | Type       | Description                                                                 | Default Value |
|---------------|------------|-----------------------------------------------------------------------------|---------------|
| `-i, --input` | `Path`     | Path to the input manifest file.                                            | None          |
| `-o, --output`| `str`      | Prefix name for the output files (without extension).                       | None          |

### Example Commands

#### `make-manifest`

```bash
python manifest.py make-manifest -i input_manifest.xlsx -o updated_manifest --remove-collection-date -a XS2
```

#### `update-manifest`

```bash
python manifest.py update-manifest -i legacy_manifest.xlsx -o updated_legacy_manifest
```

## Input File Requirements

### Required Columns

The input manifest file must contain the following columns:

- `CMO Patient ID`
- `CMO Sample Name`
- `Sample Type`

For legacy input files, the following additional columns are required:

- `cmo_patient_id`
- `cmo_sample_id_normal`
- `cmo_sample_id_plasma`

### Date Format

The script supports the following date formats:

- `MM/DD/YY`
- `M/D/YY`
- `MM/D/YYYY`
- `YYYY/MM/DD`
- `YYYY-MM-DD`

Invalid or missing dates will raise an error unless the `--remove-collection-date` option is used.

## Outputs

The script generates two output files:

1. **Excel File**: `<output_prefix>.xlsx`
2. **CSV File**: `<output_prefix>.csv`

Both files contain the updated manifest with the following columns:

- `cmo_patient_id`
- `cmo_sample_id_plasma`
- `cmo_sample_id_normal`
- `bam_path_normal`
- `bam_path_plasma_duplex`
- `bam_path_plasma_simplex`
- `maf_path`
- `cna_path`
- `sv_path`
- `paired`
- `sex`
- `collection_date`
- `dmp_patient_id`

## Script Workflow

1. **Input Validation**:
   - Checks for required columns and missing values.
   - Validates date formats.

2. **Path Generation**:
   - Generates paths for BAM, MAF, CNA, and SV files based on sample type and assay type.

3. **DataFrame Creation**:
   - Creates separate DataFrames for normal and non-normal samples.
   - Merges the DataFrames to include paired and unpaired samples.

4. **Output**:
   - Saves the updated manifest in Excel and CSV formats.

## Error Handling

The script includes error handling for the following scenarios:

- Missing required columns.
- Missing or invalid date values.
- File read/write errors.

## Example Workflow

1. **Prepare Input Manifest**:
   Ensure the input manifest file contains the required columns and valid date formats.

2. **Run `make-manifest`**:
   ```bash
   python manifest.py make-manifest -i input_manifest.xlsx -o updated_manifest --remove-collection-date -a XS2
   ```

3. **Check Outputs**:
   Verify the generated Excel and CSV files in the specified output directory.

## Contact

For questions or issues, please contact:

- **Author**: Carmelina Charalambous, Ronak Shah (@rhshah)
- **Date**: June 21, 2024