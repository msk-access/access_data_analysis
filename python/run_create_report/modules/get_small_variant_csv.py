import typer
import glob


def get_small_variant_csv(patient_id, csv_path):

    csv_file = glob.glob(csv_path.joinpath(f"{patient_id}*.csv").as_posix())
    if len(csv_file) == 0:
        typer.secho(
            "get_small_variant_csv:Could not find the csv_file using patient id %s in directory %s",
            patient_id,
            csv_path,
            err=True,
            fg=typer.colors.BRIGHT_RED,
        )
        raise typer.Abort()
    elif len(csv_file) > 1:
        typer.secho(
            "get_small_variant_csv:Found more then one csv_file using patient id %s in directory %s",
            patient_id,
            csv_path,
            err=True,
            fg=typer.colors.BRIGHT_RED,
        )
        raise typer.Abort()
    else:
        return csv_file[0]
