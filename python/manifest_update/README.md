# ACCESS Manifest Update Script

This script processes an ACCESS manifest file to generate paths for various data types and saves the updated manifest in both Excel and CSV formats. It is designed to handle Protected Health Information (PHI) and ensures compliance with privacy regulations.

## Features

- Generates paths for BAM, MAF, CNA, and SV files based on the input manifest.
- Handles both normal and non-normal sample types.
- Allows removal of the `Collection Date` column to scrub PHI.
- Outputs the updated manifest in both `.xlsx` and `.csv` formats.

## Requirements

- Python 3.8 or higher
- Required Python packages:
  - `pandas`
  - `typer`
  - `rich`
  - `arrow`
  - `openpyxl` (for Excel file handling)

## Installation

1. Clone the repository:

   ```bash
   git clone https://github.com/your-repo/access_data_analysis.git
   cd access_data_analysis/python/manifest_update
    ```

2. Create a virtual environment and install dependencies:

   ```bash
    python3 -m venv venv
    source venv/bin/activate
    pip install -r requirements.txt
    ```

## Usage

Run the script using the following command:
```bash
python manifest_update.py --input <input_manifest.xlsx> --output <output_prefix> [--remove-collection-date]
```
Replace `<input_manifest.xlsx>` with the path to your input manifest file and `<output_prefix>` with the desired prefix for the output files.
The `--remove-collection-date` flag is optional. If specified, the script will remove the `Collection Date` column from the output manifest.
### Command Line Arguments

```bash
python manifest_update.py --help
```
```plaintext
Arguments
--input, -i: Path to the input manifest file (required).
--output, -o: Prefix name for the output files (required).
--remove-collection-date: Optional flag to remove the Collection Date column from the output manifest.
--help: Show help message and exit.
```
### Example
This command processes the manifest.xlsx file, removes the Collection Date column, and saves the updated manifest as updated_manifest.xlsx and updated_manifest.csv.
```bash
python manifest_update.py --input manifest.xlsx --output updated_manifest --remove-collection-date
```
## Code Overview
The script is organized into several functions to handle different tasks:
- `generate_paths`: Generates file paths for BAM, MAF, CNA, and SV files based on the sample type and other metadata.
- `create_new_dataframe`: Creates a new DataFrame for normal or non-normal samples. Handles flexible subsetting for non-normal sample types.
- `scrub_datetime`: Replaces datetime values with None to scrub PHI.
- `make_manifest`: Main function to process the input manifest and generate the updated manifest files.



### generate_paths(row)

Generates file paths for BAM, MAF, CNA, and SV files based on the sample type and other metadata.

### create_new_dataframe(df, sample_type="Normal")

Creates a new DataFrame for normal or non-normal samples. Handles flexible subsetting for non-normal sample types.

### scrub_datetime(value)

Replaces datetime values with None to scrub PHI.

### make_manifest

Main function to process the input manifest and generate the updated manifest files.

## Output

The script generates two output files:

* <output_prefix>.xlsx: Updated manifest in Excel format.
* <output_prefix>.csv: Updated manifest in CSV format.

## Error Handling

* The script checks for required columns in the input manifest and raises an error if any are missing.
* If the Collection Date column is missing and --remove-collection-date is specified, the column is created and filled with None.

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Authors

* Carmelina Charalambous
* Ronak Shah (@rhshah)