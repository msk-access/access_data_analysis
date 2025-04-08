"""
Author: Carmelina Charalambous
Date: June 21 2024
Description: Fills in ACCESS manifest file with specific paths ready for access_data_analysis
"""

import pandas as pd
import typer
from pathlib import Path
from rich.console import Console
from rich.panel import Panel

app = typer.Typer()
console = Console()


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
        console.print(Panel(f"Required column '{column}' not found in the manifest file.", title="Error", border_style="red"))
        raise typer.Exit(code=1)
    if df[column].isnull().any():
        console.print(Panel(f"Required column '{column}' contains missing values.", title="Error", border_style="red"))
        raise typer.Exit(code=1)


@app.command()
def update_manifest(
    input_file: Path = typer.Option(..., "--input", "-i", help="Path to the input manifest file", exists=True, file_okay=True, dir_okay=False, readable=True),
    output_prefix: str = typer.Option(..., "--output", "-o", help="Prefix name for the output files (without extension)"),
):
    """
    Fills in an ACCESS manifest file with specific paths ready for access_data_analysis.
    """

    console.rule("[bold blue]Updating Manifest File[/]")

    try:
        # Load input manifest
        console.log(f"Loading manifest file: {input_file}")
        df = pd.read_excel(input_file)

        # Determine the correct patient ID column
        if 'cmo_patient_id' in df.columns:
            patient_id_col = 'cmo_patient_id'
        elif 'CMO Patient ID' in df.columns:
            patient_id_col = 'CMO Patient ID'
        else:
            console.print(Panel("Neither 'cmo_patient_id' nor 'CMO Patient ID' column found in the manifest file.", title="Error", border_style="red"))
            raise typer.Exit(code=1)

        # add base paths for normal/duplex/simplex bams and snv/cnv/sv files
        base_path_bams = "/work/access/production/data/bams/"
        base_path_small_variants = "/work/access/production/data/small_variants/"
        base_path_copy_number_variants = "/work/access/production/data/copy_number_variants/"
        base_path_structural_variants = "/work/access/production/data/structural_variants/"

        # Check if required columns exist and have values
        required_columns = [patient_id_col, 'cmo_sample_id_normal', 'cmo_sample_id_plasma']
        for col in required_columns:
            check_column(df, col)

        # Fill in bam normal
        df['bam_path_normal'] = base_path_bams + df[patient_id_col] + '/' + df['cmo_sample_id_normal'] + '/current/' + df['cmo_sample_id_normal'] + '_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX.bam'

        # Fill in bam duplex
        df['bam_path_plasma_duplex'] = base_path_bams + df[patient_id_col] + '/' + df['cmo_sample_id_plasma'] + '/current/' + df['cmo_sample_id_plasma'] + '_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.bam'

        #  Fill in bam  simplex
        df['bam_path_plasma_simplex'] = base_path_bams + df[patient_id_col] + '/' + df['cmo_sample_id_plasma'] + '/current/' + df['cmo_sample_id_plasma'] + '_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-simplex.bam'

        # Fill maf
        df['maf_path'] = base_path_small_variants + df[patient_id_col] + '/' + df['cmo_sample_id_plasma'] + '/current/' + df['cmo_sample_id_plasma'] + '.DONOR22-TP.combined-variants.vep_keptrmv_taggedHotspots_fillout_filtered.maf'

        # Fill in cna
        df['cna_path'] = base_path_copy_number_variants + df[patient_id_col] + '/' + df['cmo_sample_id_plasma'] + '/current/' + df['cmo_sample_id_plasma'] + '_copynumber_segclusp.genes.txt'

        # Fill in sv
        df['sv_path'] = base_path_structural_variants + df[patient_id_col] + '/' + df['cmo_sample_id_plasma'] + '/current/' + df['cmo_sample_id_plasma'] + '_AllAnnotatedSVs.txt'

        # Construct output file paths
        output_file_path_xlsx = f"{output_prefix}.xlsx"
        output_file_path_csv = f"{output_prefix}.csv"

        # Save the updated manifest to both Excel and CSV formats
        console.log(f"Saving updated manifest to: {output_file_path_xlsx}")
        df.to_excel(output_file_path_xlsx, index=False)
        console.log(f"Saving updated manifest to: {output_file_path_csv}")
        df.to_csv(output_file_path_csv, index=False)

        console.print(Panel(f"Saved updated manifest as [bold]{output_file_path_xlsx}[/] and [bold]{output_file_path_csv}[/]", title="Success", border_style="green"))

    except Exception as e:
        console.print_exception(show_locals=True)
        raise typer.Exit(code=1)

@app.command()
def make_manifest(
    input_file: Path = typer.Option(..., "--input", "-i", help="Path to the input manifest file", exists=True, file_okay=True, dir_okay=False, readable=True),
    output_prefix: str = typer.Option(..., "--output", "-o", help="Prefix name for the output files (without extension)"),
):
    """
    Fills in an ACCESS manifest file with specific paths ready for access_data_analysis.
    """

    console.rule("[bold blue]Updating Manifest File[/]")

    try:
        # Load input manifest
        console.log(f"Loading manifest file: {input_file}")
        df = pd.read_excel(input_file)

        # Determine the correct patient ID column
        if 'cmo_patient_id' in df.columns:
            patient_id_col = 'cmo_patient_id'
        elif 'CMO Patient ID' in df.columns:
            patient_id_col = 'CMO Patient ID'
        else:
            console.print(Panel("Neither 'cmo_patient_id' nor 'CMO Patient ID' column found in the manifest file.", title="Error", border_style="red"))
            raise typer.Exit(code=1)

        # add base paths for normal/duplex/simplex bams and snv/cnv/sv files
        base_path_bams = "/work/access/production/data/bams/"
        base_path_small_variants = "/work/access/production/data/small_variants/"
        base_path_copy_number_variants = "/work/access/production/data/copy_number_variants/"
        base_path_structural_variants = "/work/access/production/data/structural_variants/"

        # Check if required columns exist and have values
        required_columns = [patient_id_col, 'cmo_sample_id_normal', 'cmo_sample_id_plasma']
        for col in required_columns:
            check_column(df, col)

        # Fill in bam normal
        df['bam_path_normal'] = base_path_bams + df[patient_id_col] + '/' + df['cmo_sample_id_normal'] + '/current/' + df['cmo_sample_id_normal'] + '_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX.bam'

        # Fill in bam duplex
        df['bam_path_plasma_duplex'] = base_path_bams + df[patient_id_col] + '/' + df['cmo_sample_id_plasma'] + '/current/' + df['cmo_sample_id_plasma'] + '_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.bam'

        #  Fill in bam  simplex
        df['bam_path_plasma_simplex'] = base_path_bams + df[patient_id_col] + '/' + df['cmo_sample_id_plasma'] + '/current/' + df['cmo_sample_id_plasma'] + '_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-simplex.bam'

        # Fill maf
        df['maf_path'] = base_path_small_variants + df[patient_id_col] + '/' + df['cmo_sample_id_plasma'] + '/current/' + df['cmo_sample_id_plasma'] + '.DONOR22-TP.combined-variants.vep_keptrmv_taggedHotspots_fillout_filtered.maf'

        # Fill in cna
        df['cna_path'] = base_path_copy_number_variants + df[patient_id_col] + '/' + df['cmo_sample_id_plasma'] + '/current/' + df['cmo_sample_id_plasma'] + '_copynumber_segclusp.genes.txt'

        # Fill in sv
        df['sv_path'] = base_path_structural_variants + df[patient_id_col] + '/' + df['cmo_sample_id_plasma'] + '/current/' + df['cmo_sample_id_plasma'] + '_AllAnnotatedSVs.txt'

        # Construct output file paths
        output_file_path_xlsx = f"{output_prefix}.xlsx"
        output_file_path_csv = f"{output_prefix}.csv"

        # Save the updated manifest to both Excel and CSV formats
        console.log(f"Saving updated manifest to: {output_file_path_xlsx}")
        df.to_excel(output_file_path_xlsx, index=False)
        console.log(f"Saving updated manifest to: {output_file_path_csv}")
        df.to_csv(output_file_path_csv, index=False)

        console.print(Panel(f"Saved updated manifest as [bold]{output_file_path_xlsx}[/] and [bold]{output_file_path_csv}[/]", title="Success", border_style="green"))

    except Exception as e:
        console.print_exception(show_locals=True)
        raise typer.Exit(code=1)

if __name__ == "__main__":
    app()