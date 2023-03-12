import typer
import glob
import pandas as pd
from pathlib import Path
from modules.read_manifest import read_manifest


def generate_facet_maf_path(facet_path, patient_id, sample_id, best_fit):
    """Get path of maf associated with facet-suite output

    Args:
        facet_path (pathlib.PATH|str): path to search for the facet file
        patient_id (str): patient id to be used to search, default is set to None
        sample_id (str): sample id to be used to search, default is set to None
        best_fit(bool) : if true attempt to get get best fit from facet repo

    Returns:
        str: path of the facets maf
    """

    if best_fit:
        manifest_path = (
            facet_path.joinpath(
                patient_id[:7], f"{sample_id}*", "facets_review.manifest"
            )
            if sample_id
            else facet_path.joinpath(
                patient_id[:7], f"{patient_id}*", "facets_review.manifest"
            )
        )
        manifest_path = glob.glob(manifest_path.as_posix())
        if len(manifest_path) == 0:        
            if sample_id:
                typer.secho(
                    f"Could not find the facets-suite `facets_review.manifest` file using sample id. {sample_id}",
                    err=True,
                    fg=typer.colors.BRIGHT_RED,
                )
                maf_path = get_maf_path(facet_path.joinpath(patient_id[:7], f"{sample_id}*", "default", "*[0-9].ccf.maf"), patient_id, sample_id)
            else:
                typer.secho(
                    f"Could not find the facets-suite `facets_review.manifest` file using patient id. {patient_id}",
                    err=True,
                    fg=typer.colors.BRIGHT_RED,
                )
                maf_path = get_maf_path(facet_path.joinpath(patient_id[:7], f"{patient_id}*", "default", "*[0-9].ccf.maf"), patient_id, None)
        elif len(manifest_path) > 1:
            manifest_path = [Path(i) for i in manifest_path]
            manifest_path = sorted(manifest_path, key=lambda i: str(i.stem))
            manifest_path_sorted = [str(i) for i in manifest_path]
            maf_path = (
                get_maf_path(best_fit_folder, patient_id, None)
                if (
                    best_fit_folder := get_best_fit_folder(
                        manifest_path_sorted[0]
                    )
                )
                else None
            )
        elif best_fit_folder := get_best_fit_folder(manifest_path[0]):
            maf_path = get_maf_path(best_fit_folder, patient_id, None)
        else:
            maf_path = None


    elif sample_id:
        maf_path = facet_path.joinpath(
            patient_id[:7], f"{sample_id}*", "default", "*[0-9].ccf.maf"
        )
        maf_path = get_maf_path(maf_path, patient_id, sample_id)
    else:
        maf_path = facet_path.joinpath(
            patient_id[:7], f"{patient_id}*", "default", "*[0-9].ccf.maf"
        )
        maf_path = get_maf_path(maf_path, patient_id, None)
    return maf_path

def get_maf_path(maf_path, patient_id, sample_id):
    """Get the path to the maf file

    Args:
        maf_path (pathlib.Path): Base path of the maf file
        patient_id (str): DMP Patient ID for facets
        sample_id (str): DMP Sample ID if any for facets

    Returns:
        str: Path to the maf file
    """
    maf_list = glob.glob(maf_path.as_posix())
    if len(maf_list) == 0:
        if patient_id:
            typer.secho(
                f"Could not find the facets-suite MAF file using patient id. {patient_id}",
                err=True,
                fg=typer.colors.BRIGHT_RED,
            )
        if sample_id:
            typer.secho(
                f"Could not find the facets-suite MAF file using sample id. {sample_id}",
                err=True,
                fg=typer.colors.BRIGHT_RED,
            )
        return None
    elif len(maf_list) > 1:
        maf_list = [Path(i) for i in maf_list]
        maf_list_sorted = sorted(maf_list, key=lambda i: str(i.stem))
        maf_list_sorted = [str(i) for i in maf_list]
        return maf_list_sorted[0]
    else:
        maf_list_sorted = maf_list
        return maf_list_sorted[0]
    
def get_best_fit_folder(facet_manifest_path):
    """Get the best fit folder for the given facet manifest path

    Args:
        facet_manifest_path (str): manifest path to be used for determining best fit

    Returns:
        pathlib.Path: path to the folder containing best fit maf files
    """
    facet_manifest_path = Path(facet_manifest_path)
    base_path = facet_manifest_path.parent
    facet_manifest = read_manifest(facet_manifest_path)
    facet_manifest[['date_reviewed', 'time_reviewed']] = facet_manifest.date_reviewed.str.split(" ", expand = True)
    facet_manifest['date_reviewed'] = pd.to_datetime(facet_manifest['date_reviewed'])
    facet_manifest.sort_values(by='date_reviewed',ascending=False)
    folder_name = facet_manifest['fit_name'].iloc[0]
    return (
        base_path.joinpath(folder_name, "*[0-9].ccf.maf")
        if "default" in folder_name or "alt" in folder_name
        else base_path.joinpath("default", "*[0-9].ccf.maf")
    )
    
    