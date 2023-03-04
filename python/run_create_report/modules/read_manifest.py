import pandas as pd


def read_manifest(manifest):
    """_summary_

    Args:
        manifest (pathlib.PATH): _description_

    Returns:
        data_frame: _description_
    """
    return pd.read_csv(manifest, sep="\t", low_memory=False)
