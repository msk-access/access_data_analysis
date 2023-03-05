import pandas as pd


def read_manifest(manifest):
    """_summary_

    Args:
        manifest (pathlib.PATH): _description_

    Returns:
        data_frame: _description_
    """
    if(manifest.suffix == ".csv"):
        return pd.read_csv(manifest, sep=',', low_memory=False)
    else:
        return pd.read_csv(manifest, sep='\t', low_memory=False)
