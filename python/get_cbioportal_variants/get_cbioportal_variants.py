from pathlib import Path
from typing import List, Optional
from bed_lookup import BedFile
import typer
import pandas as pd

def main(
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
):

    """
    Tool to do the following operations:
    A. Get subset of variants based on Tumor_Sample_Barcode in MAF file
    B. Mark the variants as overlapping with BED file as covered [yes/no], by appending "covered" column to the subset MAF

    Requirement:
    pandas; typing; typer; bed_lookup(https://github.com/msk-access/python_bed_lookup)

    """
    if not ids:
        typer.echo("Identifiers were not provided in a text file")
        if not sid:
            typer.echo("Identifiers were not provided via command line as well")
            raise typer.Abort()

    # Read maf files
    skip = get_row(maf)
    typer.echo("Skipping Rows:", skip)
    maf_df = pd.read_csv(maf, sep="\t", skiprows=skip, low_memory=False)
    # Read Identifiers
    if not sid:
        file = open(ids)
        sid = file.read().splitlines()[1:]
        file.close()
    # filter for ids
    ns = set(sid)
    pattern = "|".join([r"\b{}\b".format(i) for i in ns])
    result = maf_df[maf_df["Tumor_Sample_Barcode"].str.contains(pattern, regex=True)]
    results_covered = result.copy(deep=True)
    results_covered["Chromosome"] = results_covered["Chromosome"].apply(str)
    # Read bed file
    b = BedFile(bed.as_posix())
    # Our chromosome column is 'Chromosome' and position column is 'Start_Position'.
    results_covered["covered"] = b.lookup_df(
        results_covered, "Chromosome", "Start_Position"
    )
    results_covered.loc[results_covered["covered"].notnull(), "covered"] = "yes"
    results_covered.loc[results_covered["covered"].notna(), "covered"] = "yes"
    results_covered.loc[results_covered["covered"].isnull(), "covered"] = "no"
    results_covered.loc[results_covered["covered"].isna(), "covered"] = "no"
    results_covered.drop_duplicates().to_csv(output_file, sep="\t", index=False)


# preprocessing
def get_row(file):
    skipped = []
    with open(file, "r") as csv_file:
        skipped.extend(i for i, line in enumerate(csv_file) if line.startswith("#"))
    return skipped


if __name__ == "__main__":
    typer.run(main)
