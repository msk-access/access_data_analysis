from pathlib import Path
from typing import List, Optional
from bed_lookup import BedFile
import typer
import pandas as pd

app = typer.Typer()


@app.command()
def subset_cpt(
    cpt: Path = typer.Option(
        "/work/access/production/resources/cbioportal/current/msk_solid_heme/data_clinical_patient.txt",
        "--cpt",
        "-p",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Clinical Patient file generated by cBioportal repo",
    ),
    ids: Path = typer.Option(
        "",
        "--ids",
        "-i",
        help="List of ids to search for in the 'PATIENT_ID' column. Header of this file is 'sample_id'",
    ),
    sid: Optional[List[str]] = typer.Option(
        "",
        help="Identifiers to search for in the 'PATIENT_ID' column. Can be given multiple times",
    ),
    output_file: str = typer.Option(
        "output_clinical_patient.txt",
        "--name",
        "-n",
        help="Name of the output file",
    ),
    col_name: str = typer.Option(
        "PATIENT_ID",
        "--cname",
        "-c",
        help="Name of the column header to be used for sub-setting",
    ),
):
    """
    Subset data_clinical_patient.txt file for given set of patient ids.

    Tool to do the following operations:
    A. Get subset of clinical information for samples based on PATIENT_ID in data_clinical_patient.txt file

    Requirement:
    pandas; typing; typer; bed_lookup(https://github.com/msk-access/python_bed_lookup)

    """
    typer.echo("Subset Clinical Patients...")
    if not ids:
        typer.echo("Identifiers were not provided in a text file")
        if not sid:
            typer.echo("Identifiers were not provided via command line as well")
            raise typer.Abort()

    cpt_df = read_tsv(cpt)
    ids_to_subset = read_ids(sid, ids)
    subset_tsv = filter_by_rows(ids_to_subset, cpt_df, col_name)
    subset_tsv.drop_duplicates().to_csv(output_file, sep="\t", index=False)
    typer.echo("Done!")


@app.command()
def subset_cst(
    cst: Path = typer.Option(
        "/work/access/production/resources/cbioportal/current/msk_solid_heme/data_clinical_sample.txt",
        "--cst",
        "-s",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Clinical Sample file generated by cBioportal repo",
    ),
    ids: Path = typer.Option(
        "",
        "--ids",
        "-i",
        help="List of ids to search for in the 'SAMPLE_ID' column. Header of this file is 'sample_id'",
    ),
    sid: Optional[List[str]] = typer.Option(
        "",
        help="Identifiers to search for in the 'SAMPLE_ID' column. Can be given multiple times",
    ),
    output_file: str = typer.Option(
        "output_clinical_samples.txt",
        "--name",
        "-n",
        help="Name of the output file",
    ),
    col_name: str = typer.Option(
        "SAMPLE_ID",
        "--cname",
        "-c",
        help="Name of the column header to be used for sub-setting",
    ),
):
    """
    Subset data_clinical_samples.txt file for given set of sample ids.

    Tool to do the following operations:
    A. Get subset of clinical information for samples based on SAMPLE_ID in data_clinical_sample.txt file

    Requirement:
    pandas; typing; typer; bed_lookup(https://github.com/msk-access/python_bed_lookup)

    """
    typer.echo("Subset Clinical Samples...")
    if not ids:
        typer.echo("Identifiers were not provided in a text file")
        if not sid:
            typer.echo("Identifiers were not provided via command line as well")
            raise typer.Abort()

    cst_df = read_tsv(cst)
    ids_to_subset = read_ids(sid, ids)
    subset_tsv = filter_by_rows(ids_to_subset, cst_df, col_name)
    subset_tsv.drop_duplicates().to_csv(output_file, sep="\t", index=False)
    typer.echo("Done!")


@app.command()
def subset_cna(
    cna: Path = typer.Option(
        "/work/access/production/resources/cbioportal/current/msk_solid_heme/data_CNA.txt",
        "--cna",
        "-c",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Copy Number Variant file generated by cBioportal repo",
    ),
    ids: Path = typer.Option(
        "",
        "--ids",
        "-i",
        help="List of ids to search for in the 'header' of the file. Header of this file is 'sample_id'",
    ),
    sid: Optional[List[str]] = typer.Option(
        "",
        help="Identifiers to search for in the 'header' of the file. Can be given multiple times",
    ),
    output_file: str = typer.Option(
        "output_CNA.txt",
        "--name",
        "-n",
        help="Name of the output file",
    ),
):
    """
    Subset data_CNA.txt file for given set of sample ids.

    Tool to do the following operations:
    A. Get subset of samples based on column header in data_CNA.txt file

    Requirement:
    pandas; typing; typer; bed_lookup(https://github.com/msk-access/python_bed_lookup)

    """
    typer.echo("Subset Copy Number Variants...")
    if not ids:
        typer.echo("Identifiers were not provided in a text file")
        if not sid:
            typer.echo("Identifiers were not provided via command line as well")
            raise typer.Abort()

    cna_df = read_tsv(cna)
    ids_to_subset = read_ids(sid, ids)
    ids_to_subset.insert(0, "Hugo_Symbol")
    subset_tsv = filter_by_columns(ids_to_subset, cna_df)
    subset_tsv.drop_duplicates().to_csv(output_file, sep="\t", index=False)
    typer.echo("Done!")


@app.command()
def subset_sv(
    sv: Path = typer.Option(
        "/work/access/production/resources/cbioportal/current/msk_solid_heme/data_sv.txt",
        "--sv",
        "-s",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Structural Variant file generated by cBioportal repo",
    ),
    ids: Path = typer.Option(
        "",
        "--ids",
        "-i",
        help="List of ids to search for in the 'Sample_ID' column. Header of this file is 'sample_id'",
    ),
    sid: Optional[List[str]] = typer.Option(
        "",
        help="Identifiers to search for in the 'Sample_ID' column. Can be given multiple times",
    ),
    output_file: str = typer.Option(
        "output_sv.txt",
        "--name",
        "-n",
        help="Name of the output file",
    ),
    col_name: str = typer.Option(
        "Sample_ID",
        "--cname",
        "-c",
        help="Name of the column header to be used for sub-setting",
    ),
):
    """
    Subset data_sv.txt file for given set of sample ids.

    Tool to do the following operations:
    A. Get subset of structural variants based on Sample_ID in data_sv.txt file

    Requirement:
    pandas; typing; typer; bed_lookup(https://github.com/msk-access/python_bed_lookup)

    """
    typer.echo("Subset Structural Variants...")
    if not ids:
        typer.echo("Identifiers were not provided in a text file")
        if not sid:
            typer.echo("Identifiers were not provided via command line as well")
            raise typer.Abort()

    sv_df = read_tsv(sv)
    ids_to_subset = read_ids(sid, ids)
    subset_tsv = filter_by_rows(ids_to_subset, sv_df, col_name)
    subset_tsv.drop_duplicates().to_csv(output_file, sep="\t", index=False)
    typer.echo("Done!")


@app.command()
def subset_maf(
    maf: Path = typer.Option(
        "/work/access/production/resources/cbioportal/current/msk_solid_heme/data_mutations_extended.txt",
        "--maf",
        "-m",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="MAF file generated by cBioportal repo",
    ),
    ids: Path = typer.Option(
        "",
        "--ids",
        "-i",
        help="List of ids to search for in the 'Tumor_Sample_Barcode' column. Header of this file is 'sample_id'",
    ),
    sid: Optional[List[str]] = typer.Option(
        "",
        help="Identifiers to search for in the 'Tumor_Sample_Barcode' column. Can be given multiple times",
    ),
    bed: Path = typer.Option(
        "/work/access/production/resources/msk-access/current/regions_of_interest/current/MSK-ACCESS-v1_0-probe-A.sorted.bed",
        "--bed",
        "-b",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="BED file to find overlapping variants",
    ),
    output_file: str = typer.Option(
        "output.maf",
        "--name",
        "-n",
        help="Name of the output file",
    ),
    col_name: str = typer.Option(
        "Tumor_Sample_Barcode",
        "--cname",
        "-c",
        help="Name of the column header to be used for sub-setting",
    ),
):

    """
    Subset MAF/TSV file and mark if an alteration is covered by BED file or not

    Tool to do the following operations:
    A. Get subset of variants based on Tumor_Sample_Barcode in data_mutations_extended.txt file
    B. Mark the variants as overlapping with BED file as covered [yes/no], by appending "covered" column to the subset MAF

    Requirement:
    pandas; typing; typer; bed_lookup(https://github.com/msk-access/python_bed_lookup)

    """
    typer.echo("Subset and Annotate MAF...")
    if not ids:
        typer.echo("Identifiers were not provided in a text file")
        if not sid:
            typer.echo("Identifiers were not provided via command line as well")
            raise typer.Abort()

    maf_df = read_tsv(maf)
    ids_to_subset = read_ids(sid, ids)
    subset_maf = filter_by_rows(ids_to_subset, maf_df, col_name)
    bed_obj = read_bed(bed)
    maf_with_covered_tag = check_if_covered(bed_obj, subset_maf)
    maf_with_covered_tag["Chromosome"] = maf_with_covered_tag["Chromosome"].apply(str)
    typer.echo("Write Filtered and Annotated MAF...")
    maf_with_covered_tag.drop_duplicates().to_csv(output_file, sep="\t", index=False)
    typer.echo("Done!")


def read_tsv(tsv):
    """Read a tsv file

    Args:
        maf (File): Input MAF/tsv like format file

    Returns:
        data_frame: Output a data frame containing the MAF/tsv
    """
    typer.echo("Read TSV file...")
    skip = get_row(tsv)
    return pd.read_csv(tsv, sep="\t", skiprows=skip, low_memory=False, chunksize=1000000)


def read_ids(sid, ids):
    """make a list of ids

    Args:
        sid (tuple): Multiple ids as tuple
        ids (File): File containing multiple ids

    Returns:
        list: List containing all ids
    """
    typer.echo("Make a list of IDs...")
    if not sid:
        with open(ids) as file:
            sid = file.read().splitlines()[1:]
    return sid


def filter_by_columns(sid, tsv_df):
    """Filter data by columns

    Args:
        sid (list): list of columns to subset over
        tsv_df (data_frame): data_frame to subset from

    Returns:
        data_frame: A copy of the subset of the data_frame
    """
    typer.echo("Subset based on columns...")
    subset_tsv = tsv_df.filter(items=sid)
    return subset_tsv.copy(deep=True)


def filter_by_rows(sid, tsv_df, col_name):
    """Filter the data by rows

    Args:
        sid (list): list of row names to subset over
        tsv_df (data_frame): data_frame to subset from
        col_name (string): name of the column to filter using names in the sid

    Returns:
        data_frame: A copy of the subset of the data_frame
    """
    typer.echo("Subset based on rows...")
    ns = set(sid)
    pattern = "|".join([r"\b{}\b".format(i) for i in ns])
    result = tsv_df[tsv_df[col_name].str.contains(pattern, regex=True)]
    return result.copy(deep=True)


def read_bed(bed):
    """Read BED file using bed_lookup

    Args:
        bed (file): File ins BED format to read

    Returns:
        object : bed file object to use for filtering
    """
    typer.echo("Reading the BED file...")
    return BedFile(bed.as_posix())


def check_if_covered(bedObj, mafObj):
    """Function to check if a variant is covered in a given bed file

    Args:
        bedObj (object): BED file object to check coverage
        mafObj (data_frame): data frame to check coverage against coordinates using column 'Chromosome' and position column is 'Start_Position'

    Returns:
        data_frame: _description_
    """
    # Our chromosome column is 'Chromosome' and position column is 'Start_Position'.
    typer.echo("Annotating the MAF file...")
    mafObj["covered"] = bedObj.lookup_df(mafObj, "Chromosome", "Start_Position")
    mafObj.loc[mafObj["covered"].notnull(), "covered"] = "yes"
    mafObj.loc[mafObj["covered"].notna(), "covered"] = "yes"
    mafObj.loc[mafObj["covered"].isnull(), "covered"] = "no"
    mafObj.loc[mafObj["covered"].isna(), "covered"] = "no"
    return mafObj


# preprocessing
def get_row(tsv_file):
    """Function to skip rows

    Args:
        tsv_file (file): file to be read

    Returns:
        list: lines to be skipped
    """
    skipped = []
    with open(tsv_file, "r") as FH:
        skipped.extend(i for i, line in enumerate(FH) if line.startswith("#"))
    return skipped


if __name__ == "__main__":
    app()
