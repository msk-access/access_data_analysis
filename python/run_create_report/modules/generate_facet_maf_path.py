import typer
import glob
from pathlib import Path
from utils import num_sort
def generate_facet_maf_path(facet_path,patient_id,sample_id=None):
    """Get path of maf associated with facet-suite output

    Args:
        facet_path (pathlib.PATH): path to search for the facet file
        patient_id (str): patient id to be used to search, default is set to None
        sample_id (str): sample id to be used to search, default is set to None

    Returns:
        str: path of the facets maf
    """
    maf_path = (
        facet_path.joinpath(
            patient_id[:7], f'{sample_id}*', 'default', '/*[0-9].ccf.maf'
        )
        if sample_id
        else facet_path.joinpath(
            patient_id[:7], f'{patient_id}*', 'default', '/*[0-9].ccf.maf'
        )
    )
    maf_list = glob.glob(maf_path.as_posix())
    if len(maf_list)==0:
        if patient_id:
            typer.secho(
                "Could not find the facets-suite MAF file using patient id. %s",
                patient_id,
                err=True,
                fg=typer.colors.BRIGHT_RED,
                )
            raise typer.Abort()
        if sample_id:
            typer.secho(
                "Could not find the facets-suite MAF file using sample id. %s",
                sample_id,
                err=True,
                fg=typer.colors.BRIGHT_RED,
                )
            raise typer.Abort()
    else:
        maf_list_sorted = maf_list.sort(key=num_sort) if len(maf_list > 1) else maf_list
        return maf_list_sorted[0]
