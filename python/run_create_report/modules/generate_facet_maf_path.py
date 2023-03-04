import typer
import glob
import re
from pathlib import Path


def generate_facet_maf_path(facet_path, patient_id, sample_id=None):
    """Get path of maf associated with facet-suite output

    Args:
        facet_path (pathlib.PATH|str): path to search for the facet file
        patient_id (str): patient id to be used to search, default is set to None
        sample_id (str): sample id to be used to search, default is set to None

    Returns:
        str: path of the facets maf
    """

    if sample_id:
        maf_path = facet_path.joinpath(
            patient_id[:7], f"{sample_id}*", "default", "*[0-9].ccf.maf"
        )
    else:
        maf_path = facet_path.joinpath(
            patient_id[:7], f"{patient_id}*", "default", "*[0-9].ccf.maf"
        )

    print(maf_path)
    maf_list = glob.glob(maf_path.as_posix())
    print(maf_list)
    if len(maf_list) == 0:
        if patient_id:
            typer.secho(
                f"Could not find the facets-suite MAF file using patient id. {patient_id}",
                err=True,
                fg=typer.colors.BRIGHT_RED,
            )
            raise typer.Abort()
        if sample_id:
            typer.secho(
                f"Could not find the facets-suite MAF file using sample id. {sample_id}",
                err=True,
                fg=typer.colors.BRIGHT_RED,
            )
            raise typer.Abort()
    elif len(maf_list) > 1:
        maf_list = [Path(i) for i in maf_list]
        maf_list_sorted = sorted(maf_list, key=lambda i: int(i.stem))
        maf_list_sorted = [str(i) for i in maf_list]
    else:
        maf_list_sorted = maf_list

    return maf_list_sorted[0]


def num_sort(test_string):
    """Numeric sort of a list with 2-digit in the string

    Args:
        test_string (list[str]): list of string

    Returns:
        list: return a sorted list
    """
    return list(map(int, re.findall(r"\d[2]", test_string)))[0]
