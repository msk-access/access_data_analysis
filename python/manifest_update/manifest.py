"""
Author: Carmelina Charalambous, Ronak Shah (@rhshah)
Date: June 21 2024
Description: Fills in ACCESS manifest file with specific paths ready for access_data_analysis
"""

import pandas as pd
import typer
from pathlib import Path
from rich.console import Console
from rich.panel import Panel
import arrow
from datetime import datetime
import numpy as np
from rich import print
app = typer.Typer()
console = Console()

BASE_PATH_BAMS = "/work/access/production/data/bams/"
BASE_PATH_SMALL_VARIANTS = "/work/access/production/data/small_variants/"
BASE_PATH_COPY_NUMBER_VARIANTS = "/work/access/production/data/copy_number_variants/"
BASE_PATH_STRUCTURAL_VARIANTS = "/work/access/production/data/structural_variants/"


def check_column(df: pd.DataFrame, column: str) -> None:
    """
    Checks if a required column exists in the DataFrame and contains no missing values.

    Args:
        df (pd.DataFrame): The DataFrame to check.
        column (str): The name of the column to check.

    Raises:
        typer.Exit: If the column is not found or contains missing values.
    """
    if column not in df.columns:
        console.print(
            Panel(
                f"Required column '{column}' not found in the manifest file.",
                title="Error",
                border_style="red",
            )
        )
        raise typer.Exit(code=1)
    if df[column].isnull().any():
        console.print(
            Panel(
                f"Required column '{column}' contains missing values.",
                title="Error",
                border_style="red",
            )
        )
        raise typer.Exit(code=1)


def validate_date(date_string):
    """
    Validates and converts a date string to a date object.

    Args:
        date_string: The date string to validate.

    Returns:
        A date object representing the validated date.

    Raises:
        ValueError: If the input date string is not in a valid format.
    """
    date_format = [
        "MM/DD/YY",
        "M/D/YY",
        "MM/D/YY",
        "M/DD/YY",
        "MM/DD/YYYY",
        "YYYY/MM/DD",
        "YYYY-MM-DD",
    ]
    for fmt in date_format:
        try:
            if isinstance(date_string, pd.Timestamp):
                return date_string.date()
            if isinstance(date_string, datetime):
                return date_string.date()
            return arrow.get(date_string, fmt).date()
        except ValueError:
            pass
        except Exception as e:
            print("Something else went wrong")
    raise ValueError("no valid date format found")


def generate_paths(row, assay_type="XS2"):
    """
    Generates file paths for BAM, MAF, CNA, and SV files based on the sample type and assay type.

    Args:
        row (pd.Series): A row of the DataFrame.
        assay_type (str): The assay type, either "XS1" or "XS2".

    Returns:
        pd.Series: A Series containing the generated paths.
    """
    sample_type = row["Sample Type"]
    cmo_patient_id = row["CMO Patient ID"]
    cmo_sample_name = row["CMO Sample Name"]
    base_bam_path = (
        f"{BASE_PATH_BAMS}{cmo_patient_id}/{cmo_sample_name}/current/{cmo_sample_name}"
    )

    if sample_type == "Normal":
        return pd.Series(
            {
                "bam_path_normal": f"{base_bam_path}_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX.bam",
                "bam_path_plasma_duplex": None,
                "bam_path_plasma_simplex": None,
                "maf_path": None,
                "cna_path": None,
                "sv_path": None,
            }
        )
    else:
        maf_path = (
            f"{BASE_PATH_SMALL_VARIANTS}{cmo_patient_id}/{cmo_sample_name}/current/{cmo_sample_name}.DONOR22-TP.combined-variants.vep_keptrmv_taggedHotspots_fillout_filtered.maf"
            if assay_type == "XS1"
            else f"{BASE_PATH_SMALL_VARIANTS}{cmo_patient_id}/{cmo_sample_name}/current/{cmo_sample_name}.Donor19F21c2206-TP01.combined-variants.vep_keptrmv_taggedHotspots_fillout_filtered.maf"
        )
        return pd.Series(
            {
                "bam_path_normal": None,
                "bam_path_plasma_duplex": f"{base_bam_path}_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.bam",
                "bam_path_plasma_simplex": f"{base_bam_path}_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-simplex.bam",
                "maf_path": maf_path,
                "cna_path": f"{BASE_PATH_COPY_NUMBER_VARIANTS}{cmo_patient_id}/{cmo_sample_name}/current/{cmo_sample_name}_copynumber_segclusp.genes.txt",
                "sv_path": f"{BASE_PATH_STRUCTURAL_VARIANTS}{cmo_patient_id}/{cmo_sample_name}/current/{cmo_sample_name}_AllAnnotatedSVs.txt",
            }
        )


def create_new_dataframe(df, sample_type="Normal"):
    """
    Creates a new DataFrame with specified columns for normal or non-normal samples.
    If sample_type is "Normal", it subsets for normal samples.
    If sample_type is not "Normal", it subsets for all non-normal samples.
    """
    if sample_type == "Normal":
        subset_df = df[df["Sample Type"] == sample_type].copy()
    else:
        subset_df = df[df["Sample Type"] != "Normal"].copy()  # Exclude "Normal"

    new_df = pd.DataFrame()
    if not subset_df.empty:  # Check if subset_df is not empty
        new_df["cmo_patient_id"] = subset_df["CMO Patient ID"]
        new_df["cmo_sample_id_plasma"] = subset_df["CMO Sample Name"]
        new_df["cmo_sample_id_normal"] = subset_df["CMO Sample Name"]
        new_df["bam_path_normal"] = subset_df["bam_path_normal"]
        new_df["sex"] = subset_df["Sex"]
        new_df["collection_date"] = subset_df["Collection Date"]
        new_df["dmp_patient_id"] = subset_df["DMP_PATIENT_ID"]

        if sample_type != "Normal":
            new_df["bam_path_plasma_duplex"] = subset_df["bam_path_plasma_duplex"]
            new_df["bam_path_plasma_simplex"] = subset_df["bam_path_plasma_simplex"]
            new_df["maf_path"] = subset_df["maf_path"]
            new_df["cna_path"] = subset_df["cna_path"]
            new_df["sv_path"] = subset_df["sv_path"]

    return new_df

# Define a function to convert datetime values to None
def scrub_datetime(value):
    # Convert all values to None to scrub PHI
    return None  # or return pd.NaT if you prefer

@app.command()
def make_manifest(
    input_file: Path = typer.Option(
        ...,
        "--input",
        "-i",
        help="Path to the input manifest file",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
    ),
    output_prefix: str = typer.Option(
        ...,
        "--output",
        "-o",
        help="Prefix name for the output files (without extension)",
    ),
    remove_collection_date: bool = typer.Option(
        False,
        "--remove-collection-date",
        help="Remove collection date from the output manifest (PHI).",
    ),
    assay_type: str = typer.Option(
        "XS2",
        "--assay-type",
        "-a",
        help="Assay type, either 'XS1' or 'XS2'. Default is 'XS2'.",
    ),
):
    """
    Processes the input manifest file to generate paths for various data types
    and saves the updated manifest in both Excel and CSV formats.

    Args:
        input_file (Path): Path to the input manifest file.
        output_prefix (str): Prefix name for the output files.
        remove_collection_date (bool): If True, the collection date column will be removed from the output.
        assay_type (str): Assay type, either "XS1" or "XS2".
    """
    console.rule("[bold blue]Making Manifest File[/]")

    console.print(
        Panel(
            "[bold yellow]Warning:[/bold yellow] The 'Collection Date' column contains Protected Health Information (PHI). This tool should only be used on systems that are compliant with HIPAA and other relevant privacy regulations.",
            title="PHI Warning",
            border_style="yellow",
        )
    )

    try:
        # Load input manifest
        console.log(f"Loading manifest file: {input_file}")
        if remove_collection_date:
            df = pd.read_excel(input_file, converters={"Collection Date":scrub_datetime})
        else:
            df = pd.read_excel(input_file, converters={"Collection Date": validate_date})
        # Check if required columns exist and have values
        required_columns = ["CMO Patient ID", "CMO Sample Name", "Sample Type"]
        for col in required_columns:
            check_column(df, col)

        # Generate paths and concatenate
        paths_df = df.apply(generate_paths, axis=1)
        df = pd.concat([df, paths_df], axis=1)

        # Create new DataFrames for normal and non-normal samples
        new_normal_df = create_new_dataframe(df, "Normal")
        new_non_normal_df = create_new_dataframe(df, "")

        # Create a list to store new rows for multiple normals
        new_rows = []

        # Iterate through each non-normal sample
        for index, row in new_non_normal_df.iterrows():
            # Find corresponding normal samples for the same patient
            normal_samples = new_normal_df[
                new_normal_df["cmo_patient_id"] == row["cmo_patient_id"]
            ]

            # If normal samples are found, create a new row for each normal sample
            if not normal_samples.empty:
                for normal_index, normal_row in normal_samples.iterrows():
                    new_row = row.copy()
                    new_row["cmo_sample_id_normal"] = normal_row["cmo_sample_id_normal"]
                    new_row["bam_path_normal"] = normal_row["bam_path_normal"]
                    new_rows.append(new_row)
            else:
                # If no normal sample is found, set normal-related columns to None
                new_non_normal_df.at[index, "cmo_sample_id_normal"] = None
                new_non_normal_df.at[index, "bam_path_normal"] = None

        # Append new rows to the DataFrame
        if new_rows:
            new_non_normal_df = pd.concat(
                [new_non_normal_df, pd.DataFrame(new_rows)], ignore_index=True
            )

        # Concatenate the two dataframes
        final_df = new_non_normal_df
        # Group and expand DataFrame for plasma samples with normal samples
        grouped = final_df.groupby(["cmo_sample_id_plasma", "cmo_patient_id"])

        # Create a list to store the processed rows
        processed_rows = []

        # Iterate through each group
        for (cmo_sample_id_plasma, cmo_patient_id), group_df in grouped:
            # Iterate through the columns of the group
            processed_row = {
                "cmo_patient_id": cmo_patient_id,
                "cmo_sample_id_plasma": cmo_sample_id_plasma,
            }
            for col in final_df.columns:
                if col not in ["cmo_sample_id_plasma", "cmo_patient_id"]:
                    # Get the unique values in the column for the group
                    unique_values = group_df[col].unique().tolist()

                    # If there's only one unique value, use it directly
                    if len(unique_values) == 1:
                        processed_row[col] = unique_values[0]
                    # If there are multiple unique values, take the second one
                    elif len(unique_values) > 1:
                        processed_row[col] = unique_values[1]
                    # If there are no unique values, set to None
                    else:
                        processed_row[col] = None
            processed_rows.append(processed_row)

        # Create the final DataFrame from the processed rows
        t_final_df = pd.DataFrame(processed_rows)
        # Add the "paired" column
        t_final_df["paired"] = t_final_df["cmo_sample_id_normal"].apply(
            lambda x: "Paired" if pd.notna(x) else "Unpaired"
        )

        # Define the desired column order
        column_order = [
            "cmo_patient_id",
            "cmo_sample_id_plasma",
            "cmo_sample_id_normal",
            "bam_path_normal",
            "paired",
            "sex",
            "collection_date",
            "dmp_patient_id",
            "bam_path_plasma_duplex",
            "bam_path_plasma_simplex",
            "maf_path",
            "cna_path",
            "sv_path",
        ]

        # Reorder the columns
        t_final_df = t_final_df[column_order]

        # Save the updated manifest to both Excel and CSV formats
        output_file_path_xlsx = f"{output_prefix}.xlsx"
        output_file_path_csv = f"{output_prefix}.csv"

        console.log(f"Saving updated manifest to: {output_file_path_xlsx}")
        t_final_df.to_excel(output_file_path_xlsx, index=False)
        console.log(f"Saving updated manifest to: {output_file_path_csv}")
        t_final_df.to_csv(output_file_path_csv, index=False)

        console.print(
            Panel(
                f"Saved updated manifest as [bold]{output_file_path_xlsx}[/] and [bold]{output_file_path_csv}[/]",
                title="Success",
                border_style="green",
            )
        )

    except Exception as e:
        console.print_exception(show_locals=True)
        raise typer.Exit(code=1) from e


@app.command()
def update_manifest(
    input_file: Path = typer.Option(
        ...,
        "--input",
        "-i",
        help="Path to the input manifest file",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
    ),
    output_prefix: str = typer.Option(
        ...,
        "--output",
        "-o",
        help="Prefix name for the output files (without extension)",
    ),
):
    """
    Fills in an ACCESS manifest file with specific paths ready for access_data_analysis.
    This version handles the legacy input file format.
    """

    console.rule("[bold blue]Updating Manifest File (Legacy)[/]")

    try:
        # Load input manifest
        console.log(f"Loading manifest file: {input_file}")
        df = pd.read_excel(input_file)

        # Check if required columns exist and have values
        required_columns = [
            "cmo_patient_id",
            "cmo_sample_id_normal",
            "cmo_sample_id_plasma",
        ]
        for col in required_columns:
            check_column(df, col)

        # Fill in bam normal
        df["bam_path_normal"] = (
            BASE_PATH_BAMS
            + df["cmo_patient_id"]
            + "/"
            + df["cmo_sample_id_normal"]
            + "/current/"
            + df["cmo_sample_id_normal"]
            + "_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX.bam"
        )

        # Fill in bam duplex
        df["bam_path_plasma_duplex"] = (
            BASE_PATH_BAMS
            + df["cmo_patient_id"]
            + "/"
            + df["cmo_sample_id_plasma"]
            + "/current/"
            + df["cmo_sample_id_plasma"]
            + "_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.bam"
        )

        #  Fill in bam  simplex
        df["bam_path_plasma_simplex"] = (
            BASE_PATH_BAMS
            + df["cmo_patient_id"]
            + "/"
            + df["cmo_sample_id_plasma"]
            + "/current/"
            + df["cmo_sample_id_plasma"]
            + "_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-simplex.bam"
        )

        # Fill maf
        df["maf_path"] = (
            BASE_PATH_SMALL_VARIANTS
            + df["cmo_patient_id"]
            + "/"
            + df["cmo_sample_id_plasma"]
            + "/current/"
            + df["cmo_sample_id_plasma"]
            + ".DONOR22-TP.combined-variants.vep_keptrmv_taggedHotspots_fillout_filtered.maf"
        )

        # Fill in cna
        df["cna_path"] = (
            BASE_PATH_COPY_NUMBER_VARIANTS
            + df["cmo_patient_id"]
            + "/"
            + df["cmo_sample_id_plasma"]
            + "/current/"
            + df["cmo_sample_id_plasma"]
            + "_copynumber_segclusp.genes.txt"
        )

        # Fill in sv
        df["sv_path"] = (
            BASE_PATH_STRUCTURAL_VARIANTS
            + df["cmo_patient_id"]
            + "/"
            + df["cmo_sample_id_plasma"]
            + "/current/"
            + df["cmo_sample_id_plasma"]
            + "_AllAnnotatedSVs.txt"
        )

        # Construct output file paths
        output_file_path_xlsx = f"{output_prefix}.xlsx"
        output_file_path_csv = f"{output_prefix}.csv"

        # Save the updated manifest to both Excel and CSV formats
        console.log(f"Saving updated manifest to: {output_file_path_xlsx}")
        df.to_excel(output_file_path_xlsx, index=False)
        console.log(f"Saving updated manifest to: {output_file_path_csv}")
        df.to_csv(output_file_path_csv, index=False)

        console.print(
            Panel(
                f"Saved updated manifest as [bold]{output_file_path_xlsx}[/] and [bold]{output_file_path_csv}[/]",
                title="Success",
                border_style="green",
            )
        )

    except Exception as e:
        console.print_exception(show_locals=True)
        raise typer.Exit(code=1) from e


if __name__ == "__main__":
    app()