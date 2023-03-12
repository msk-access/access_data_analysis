import typer


def generate_repo_path(
    repo_path=None, script_path=None, template_path=None, template_days=None
):
    """Generate path to create_report.R and template RMarkdown file

    Args:
        repo_path (pathlib.Path, optional): Path to clone of git repo access_data_analysis. Defaults to None.
        script_path (pathlib.Path, optional): Path to create_report.R. Defaults to None.
        template_path (pathlib.Path, optional): Path to template RMarkdown file. Defaults to None.
        template_days (bool, optional): True|False to use days template if using repo_path. Defaults to None.

    Raises:
        typer.Abort: Abort if both repo_path and script_path are not given
        typer.Abort: Abort if both repo_path and template_path are not given

    Returns:
        str: Path to create_report.R and path to template markdown file
    """
    if not repo_path:
        typer.secho(
            "Path to access_data_analysis repo is not provided!",
            fg=typer.colors.BRIGHT_YELLOW,
        )
        if not script_path:
            typer.secho(
                "Path to create_report.R script is not provided",
                err=True,
                fg=typer.colors.BRIGHT_RED,
            )
            raise typer.Abort()
        else:
            typer.secho(
                f"Path to create_report.R script is provided as {script_path}",
                fg=typer.colors.BRIGHT_GREEN,
            )
        if not template_path:
            typer.secho(
                "Path to RMarkdown template file is not provided",
                err=True,
                fg=typer.colors.BRIGHT_RED,
            )
            raise typer.Abort()
        else:
            typer.secho(
                f"Path to RMarkdown template file is provided as {template_path}",
                fg=typer.colors.BRIGHT_GREEN,
            )
    else:
        typer.secho(
            f"Path to access_data_analysis repo is provided as {repo_path}",
            fg=typer.colors.BRIGHT_GREEN,
        )
        if script_path or template_path:
            typer.secho(
                "Both repo-path and script-path as well as template path cannot be provided, please use either repo-path or script-path along with template-path",
                err=True,
                fg=typer.colors.BRIGHT_RED,
            )
            raise typer.Abort()
        script_path = repo_path.joinpath("reports", "create_report.R")
        typer.secho(
            f"Path to create_report.R is {script_path}", fg=typer.colors.BRIGHT_GREEN
        )
        if template_days:
            template_path = repo_path.joinpath("reports", "template_days.Rmd")
            typer.secho(
                f"Path to template_days.Rmd is {template_path}",
                fg=typer.colors.BRIGHT_GREEN,
            )
        else:
            template_path = repo_path.joinpath("reports", "template.Rmd")
            typer.secho(
                f"Path to template.Rmd is {template_path}",
                fg=typer.colors.BRIGHT_GREEN,
            )
    return (script_path.as_posix(), template_path.as_posix())
