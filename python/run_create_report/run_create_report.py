from pathlib import Path
from typing import List, Optional
import typer
import pandas as pd


def main(
    template_path: Path = typer.Option(
        "",
        "--template",
        "-t",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Path to the template.Rmd or template_days.Rmd to be used with create_report.R",
    ),
    script_path: Path = typer.Option(
        "",
        "--script",
        "-s",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Path to the create_report.R script",
    ),
    repo_path: Path = typer.Option(
        "",
        "--repo",
        "-r",
        exists=True,
        file_okay=False,
        dir_okay=True,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Base path to where the git repository is located for access_data_analysis",
    ),
    template_days: bool = typer.Option(
        "",
        "--template-days",
        "-td",
        help="If the `--repo` option is specified and if this is set to True then we will use the template_days RMarkdown file as the template",
    ),
    manifest: Path = typer.Option(
        "",
        "--manifest",
        "-m",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="File containing meta information per sample.",
    ),
    facet_repo: Path = typer.Option(
        "/juno/work/ccs/shared/resources/impact/facets/all/",
        "--facet-repo",
        "-f",
        exists=True,
        file_okay=False,
        dir_okay=True,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Base path for all results of facets on Clinical MSK-IMPACT samples",
    ),
    variant_path: Path = typer.Option(
        ...,
        "--variant-results",
        "-c",
        exists=True,
        file_okay=False,
        dir_okay=True,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Base path for all results of small variants as generated by filter_calls.R script in access_data_analysis (Make sure only High Confidence calls are included)",
    ),
    cnv_path: Path = typer.Option(
        ...,
        "--cnv-results",
        "-c",
        exists=True,
        file_okay=False,
        dir_okay=True,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Base path for all results of CNV as generated by CNV_processing.R script in access_data_analysis",
    ),
    tumor_type: str = typer.Option(
        "",
        "--tumor-type",
        "-t",
        help="Tumor type label for the report",
    ),
):
    """_summary_

    Args:
        template_path (Path, optional): _description_. Defaults to typer.Option( "", "--template", "-t", exists=True, file_okay=True, dir_okay=False, writable=False, readable=True, resolve_path=True, help="Path to the template.Rmd or template_days.Rmd to be used with create_report.R", ).
        script_path (Path, optional): _description_. Defaults to typer.Option( "", "--script", "-s", exists=True, file_okay=True, dir_okay=False, writable=False, readable=True, resolve_path=True, help="Path to the create_report.R script", ).
        repo_path (Path, optional): _description_. Defaults to typer.Option( "", "--repo", "-r", exists=True, file_okay=False, dir_okay=True, writable=False, readable=True, resolve_path=True, help="Base path to where the git repository is located for access_data_analysis", ).
        manifest (Path, optional): _description_. Defaults to typer.Option( "", "--manifest", "-m", exists=True, file_okay=True, dir_okay=False, writable=False, readable=True, resolve_path=True, help="File containing meta information per sample.", ).
        facet_repo (Path, optional): _description_. Defaults to typer.Option( "/juno/work/ccs/shared/resources/impact/facets/all/", "--facet-repo", "-f", exists=True, file_okay=False, dir_okay=True, writable=False, readable=True, resolve_path=True, help="Base path for all results of facets on Clinical MSK-IMPACT samples", ).
        variant_path (Path, optional): _description_. Defaults to typer.Option( "--variant-results", "-c", exists=True, file_okay=False, dir_okay=True, writable=False, readable=True, resolve_path=True, help="Base path for all results of small variants as generated by filter_calls.R script in access_data_analysis (Make sure only High Confidence calls are included)", ).
        cnv_path (Path, optional): _description_. Defaults to typer.Option( "--cnv-results", "-c", exists=True, file_okay=False, dir_okay=True, writable=False, readable=True, resolve_path=True, help="Base path for all results of CNV as generated by CNV_processing.R script in access_data_analysis", ).
        tumor_type (str, optional): _description_. Defaults to typer.Option( "", "--tumor-type", "-t", help="Tumor type label for the report", ).

    Raises:
        typer.Warning: _description_
        typer.Abort: _description_
        typer.Abort: _description_
        typer.Abort: _description_
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
                "Path to create_report.R script is provided as %s",
                script_path,
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
                "Path to RMarkdown template file is provided as %s",
                template_path,
                fg=typer.colors.BRIGHT_GREEN,
            )
    else:
        typer.secho(
            "Path to access_data_analysis repo is provided as %s",
            repo_path,
            fg=typer.colors.BRIGHT_GREEN,
        )
        script_path = repo_path.joinpath("reports", "create_report.R")
        typer.secho(
            "Path to create_report.R is %s", script_path, fg=typer.colors.BRIGHT_GREEN
        )
        if template_days:
            template_path = repo_path.joinpath("reports", "template_days.Rmd")
            typer.secho(
                "Path to template_days.Rmd is %s",
                template_path,
                fg=typer.colors.BRIGHT_GREEN,
            )
        else:
            template_path = repo_path.joinpath("reports", "template.Rmd")
            typer.secho(
                "Path to template.Rmd is %s",
                template_path,
                fg=typer.colors.BRIGHT_GREEN,
            )


if __name__ == "__main__":
    typer.run(main)
