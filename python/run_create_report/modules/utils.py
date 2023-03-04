import pandas as pd
import re


def num_sort(test_string):
    """Numeric sort of a list with 2-digit in the string

    Args:
        test_string (list[str]): list of string

    Returns:
        list: return a sorted list
    """
    return list(map(int, re.findall(r"\d[2]", test_string)))[0]


def read_manifest(manifest):
    """_summary_

    Args:
        manifest (pathlib.PATH): _description_

    Returns:
        data_frame: _description_
    """
    return pd.read_csv(manifest, sep="\t", low_memory=False)
