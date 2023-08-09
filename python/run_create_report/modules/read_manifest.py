import pandas as pd


def read_manifest(manifest):
    """_summary_

    Args:
        manifest (pathlib.PATH): _description_

    Returns:
        data_frame: _description_
    """
    skip_rows = get_row(manifest)
    if manifest.suffix == ".csv":
        return pd.read_csv(manifest, sep=",", skiprows=skip_rows, low_memory=False, keep_default_na=False)
    else:
        return pd.read_csv(manifest, sep="\t", skiprows=skip_rows, low_memory=False, keep_default_na=False)


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
