import typer


def check_required_columns(manifest,template_days=None):
    """Check if all required columns are present in the sample manifest file

    Args:
        manifest (data_frame): meta information file with information for each sample
        template_days (bool): True|False if template days RMarkdown will be used

    Raises:
        typer.Abort: if "cmo_patient_id" column not provided
        typer.Abort: if "cmo_sample_id/sample_id" column not provided
        typer.Abort: if "dmp_patient_id" column not provided
        typer.Abort: if "collection_date/collection_day" column not provided
        typer.Abort: if "timepoint" column not provided

    Returns:
        list: column name for the manifest file
    """
    # Get the list of all column names from headers
    column_headers = list(manifest.columns.values)
    if "cmo_patient_id" not in column_headers:
        typer.secho(
            "check_required_columns:Could not find cmo_patient_id in %s",
            ",".join(column_headers),
            err=True,
            fg=typer.colors.BRIGHT_RED,
        )
        raise typer.Abort()
    if "cmo_sample_id" not in column_headers or "sample_id" not in column_headers:
        typer.secho(
            "check_required_columns:Could not find cmo_sample_id or sample_id in %s",
            ",".join(column_headers),
            err=True,
            fg=typer.colors.BRIGHT_RED,
        )
        raise typer.Abort()
    if "dmp_patient_id" not in column_headers:
        typer.secho(
            "check_required_columns:Could not find dmp_patient_id in %s",
            ",".join(column_headers),
            err=True,
            fg=typer.colors.BRIGHT_RED,
        )
        raise typer.Abort()
    if template_days:
        if "collection_day" not in column_headers:
            typer.secho(
                "check_required_columns:Could not find collection_day in %s",
                ",".join(column_headers),
                err=True,
                fg=typer.colors.BRIGHT_RED,
            )
            raise typer.Abort()
    elif "collection_date" not in column_headers:
        typer.secho(
            "check_required_columns:Could not find collection_date in %s",
            ",".join(column_headers),
            err=True,
            fg=typer.colors.BRIGHT_RED,
        )
        raise typer.Abort()
    if "timepoint" not in column_headers:
        typer.secho(
            "check_required_columns:Could not find timepoint in %s",
            ",".join(column_headers),
            err=True,
            fg=typer.colors.BRIGHT_RED,
        )
        raise typer.Abort()
    if "dmp_sample_id" in column_headers:
        typer.secho(
            "check_required_columns:dmp_sample_id is present, thus dmp_sample_id will be used for finding facet results in facet-repo",
            fg=typer.colors.GREEN,
        )
    else:
        typer.secho(
            "check_required_columns:dmp_sample_id is not present, thus dmp_patient_id will be used for finding facet results in facet-rep",
            fg=typer.colors.GREEN,
        )


    return column_headers
