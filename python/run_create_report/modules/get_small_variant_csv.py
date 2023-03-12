import glob

import typer


def get_small_variant_csv(patient_id, csv_path):
    """Get the path to CSV file to be used for a given patient containing all variants

    Args:
        patient_id (str): patient id used to identify the csv file
        csv_path (pathlib.path): base path where the csv file is expected to be present

    Raises:
        typer.Abort: if no csv file is returned
        typer.Abort: if more then one csv file is returned

    Returns:
        str: path to csv file containing the variants
    """

    csv_file = glob.glob(csv_path.joinpath(f"{patient_id}*.csv").as_posix())
    if len(csv_file) == 0:
        typer.secho(
            f"get_small_variant_csv:Could not find the csv_file using patient id {patient_id} in directory {csv_path}",
            err=True,
            fg=typer.colors.BRIGHT_RED,
        )
        raise typer.Abort()
    elif len(csv_file) > 1:
        typer.secho(
            f"get_small_variant_csv:Found more then one csv_file using patient id {patient_id} in directory {csv_path}",
            err=True,
            fg=typer.colors.BRIGHT_RED,
        )
        raise typer.Abort()
    else:
        return csv_file[0]
