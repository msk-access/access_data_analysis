from pathlib import Path

import typer
import numpy as np
from modules.check_required_columns import check_required_columns
from modules.generate_create_report_cmd import generate_create_report_cmd
from modules.generate_facet_maf_path import generate_facet_maf_path
from modules.generate_repo_paths import generate_repo_path
from modules.get_small_variant_csv import get_small_variant_csv
from modules.read_manifest import read_manifest
from modules.run_cmd import run_cmd
from rich import print
from rich.progress import Progress, SpinnerColumn, TextColumn


def main(
    repo_path: Path = typer.Option(
        None,
        "--repo",
        "-r",
        help="Base path to where the git repository is located for access_data_analysis",
    ),
    script_path: Path = typer.Option(
        None,
        "--script",
        "-s",
        help="Path to the create_report.R script, fall back if `--repo` is not given",
    ),
    template_path: Path = typer.Option(
        None,
        "--template",
        "-t",
        help="Path to the template.Rmd or template_days.Rmd to be used with create_report.R when `--repo` is not given",
    ),
    manifest: Path = typer.Option(
        ...,
        "--manifest",
        "-m",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="File containing meta information per sample. Require following columns in the header: cmo_patient_id, sample_id, dmp_patient_id, collection_date or collection_day, timepoint. If dmp_sample_id column is given and has information that will be used to run facets. If dmp_sample_id is not given and dmp_patient_id is given than it will be used to get the Tumor sample with lowest number. If dmp_sample_id or dmp_patient_id is not given then it will run without the facet maf file",
    ),
    variant_path: Path = typer.Option(
        ...,
        "--variant-results",
        "-v",
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
    best_fit: bool = typer.Option(
        False,
        "--best-fit",
        "-bf",
        help="If this is set to True then we will attempt to parse `facets_review.manifest` file to pick the best fit for a given dmp_sample_id",
    ),
    tumor_type: str = typer.Option(
        ...,
        "--tumor-type",
        "-l",
        help="Tumor type label for the report",
    ),
    copy_facet: bool = typer.Option(
        False,
        "--copy-facet-maf",
        "-cfm",
        help="If this is set to True then we will copy the facet maf file in the directory specified in `copy_facet_dir`",
    ),
    copy_facet_dir: Path = typer.Option(
        None,
        "--copy-facet-dir",
        "-cfd",
        help="Directory path where the facet maf file should be copied.",
    ),
    template_days: bool = typer.Option(
        False,
        "--template-days",
        "-d",
        help="If the `--repo` option is specified and if this is set to True then we will use the template_days RMarkdown file as the template",
    ),
    markdown: bool = typer.Option(
        False,
        "--generate-markdown",
        "-gm",
        help="If given, the create_report.R will be run with `-md` flag to generate markdown",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        "-ff",
        help="If this is set to True then we will not stop if an error is encountered in a given sample while running create_report.R but keep on running for the next sample",
    ),
):
    """Wrapper script to run create_report.R

    Args:
        repo_path (Path, optional): "Base path to where the git repository is located for access_data_analysis".
        script_path (Path, optional): "Path to the create_report.R script, fall back if `--repo` is not given".
        template_path (Path, optional): "Path to the template.Rmd or template_days.Rmd to be used with create_report.R when `--repo` is not given".
        manifest (Path, required): "File containing meta information per sample. Require following columns in the header: `cmo_patient_id`, `sample_id`, `dmp_patient_id`, `collection_date` or `collection_day`, `timepoint`. If dmp_sample_id column is given and has information that will be used to run facets. if dmp_sample_id is not given and dmp_patient_id is given than it will be used to get the Tumor sample with lowest number.If dmp_sample_id or dmp_patient_id is not given then it will run without the facet maf file".
        variant_path (Path, required): "Base path for all results of small variants as generated by filter_calls.R script in access_data_analysis (Make sure only High Confidence calls are included)".
        cnv_path (Path, required): "Base path for all results of CNV as generated by CNV_processing.R script in access_data_analysis".
        facet_repo (Path, required): "Base path for all results of facets on Clinical MSK-IMPACT samples".
        best_fit (bool, optional): "If this is set to True then we will attempt to parse `facets_review.manifest` file to pick the best fit for a given dmp_sample_id".
        tumor_type (str, required): "Tumor type label for the report".
        copy_facet (bool, optional): "If this is set to True then we will copy the facet maf file in the directory specified in `copy_facet_dir`".
        copy_facet_dir (Path, optional): "Directory path where the facet maf file should be copied.".
        template_days (bool, optional): "If the `--repo` option is specified and if this is set to True then we will use the template_days RMarkdown file as the template".
        markdown (bool, optional): "If given, the create_report.R will be run with `-md` flag to generate markdown".
        force (bool, optional): "If this is set to True then we will not stop if an error is encountered in a given sample but keep on running for the next sample".
    """
    # Read the manifest file
    manifest_df = read_manifest(manifest)
    # check required columns
    column_header, manifest_to_traverse = check_required_columns(
        manifest_df, template_days
    )
    print(
        "\nTraversing through",
        len(manifest_to_traverse),
        "patients to run create_report.R\n",
    )
    # get general paths
    (script_path, template_path) = generate_repo_path(
        repo_path, script_path, template_path, template_days
    )
    # iterate through each row and select information needed to generate the command
    summary = [
        "\t".join(
            [
                "cmo_patient_id",
                "dmp_patient_id",
                "dmp_sample_id",
                "facet_path",
                "comments",
            ]
        )
    ]
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        transient=True,
    ) as progress:
        print("\n")
        progress.add_task(description="Processing\n", total=None)
        for i in range(len(manifest_to_traverse)):
            cmo_patient_id = manifest_to_traverse.loc[i, "cmo_patient_id"]
            dmp_patient_id = manifest_to_traverse.loc[i, "dmp_patient_id"]
            typer.secho(
                f"Running for patient with CMO ID {cmo_patient_id}, and DMP ID {dmp_patient_id}",
                fg=typer.colors.BRIGHT_GREEN,
            )
            small_variants_path = get_small_variant_csv(cmo_patient_id, variant_path)
            if "dmp_sample_id" in column_header:
                dmp_sample_id = manifest_to_traverse.loc[i, "dmp_sample_id"]
                facet_path = generate_facet_maf_path(
                    facet_repo, dmp_patient_id, dmp_sample_id, best_fit
                )
            else:
                dmp_sample_id = None
                if not dmp_sample_id or dmp_sample_id != None or dmp_sample_id != '':
                    facet_path = generate_facet_maf_path(
                        facet_repo, dmp_patient_id, dmp_sample_id, best_fit
                    )
                else:
                    facet_path = None
            if not facet_path:
                typer.secho(
                    f"Running for patient with CMO ID {cmo_patient_id}, and DMP ID {dmp_patient_id} without facets maf",
                    err=True,
                    fg=typer.colors.BRIGHT_RED,
                )
            if facet_path:
                typer.secho(
                    f"Running for patient with CMO ID {cmo_patient_id}, and DMP ID {dmp_patient_id} with facets maf: {facet_path}",
                    fg=typer.colors.BRIGHT_GREEN,
                )
                facet_path = Path(facet_path)
                maf_id = facet_path.stem
                dmp_sample_id = maf_id.split("_", 1)[0]
                if copy_facet:
                    if not copy_facet_dir:
                        copy_facet_dir = Path.cwd() / "facet_files"
                        copy_facet_dir.mkdir(parents=True, exist_ok=True)
                    cp_facet_cmd = f"cp {facet_path} {copy_facet_dir.as_posix()}"
                    p1 = run_cmd(cp_facet_cmd, force)
                    typer.secho(
                        f"Done copying facet maf file for patient with CMO ID {cmo_patient_id}, and DMP ID {dmp_patient_id} and output is written in {copy_facet_dir}",
                        fg=typer.colors.BRIGHT_GREEN,
                    )
                create_report_cmd, html_output = generate_create_report_cmd(
                    script_path,
                    markdown,
                    template_path,
                    cmo_patient_id,
                    small_variants_path,
                    manifest,
                    cnv_path,
                    dmp_patient_id,
                    dmp_sample_id,
                    facet_path,
                    tumor_type,
                )
                p2 = run_cmd(create_report_cmd, force)
                if "Error" in str(p2) or "error" in str(p2):
                    summary.append(
                        "\t".join(
                            [
                                cmo_patient_id,
                                dmp_patient_id,
                                str(dmp_sample_id),
                                facet_path.as_posix(),
                                "create_report.R failed",
                            ]
                        )
                    )
                else:
                    summary.append(
                        "\t".join(
                            [
                                cmo_patient_id,
                                dmp_patient_id,
                                str(dmp_sample_id),
                                facet_path.as_posix(),
                                "create_report.R ran with facet maf",
                            ]
                        )
                    )
                typer.secho(
                    f"Done running create_report.R for patient with CMO ID {cmo_patient_id}, and DMP ID {dmp_patient_id} and output is written in {html_output}",
                    fg=typer.colors.BRIGHT_GREEN,
                )
            else:
                if copy_facet:
                    typer.secho(
                        f"No maf file to copy for patient with CMO ID {cmo_patient_id}, and DMP ID {dmp_patient_id} and thus skipping",
                        fg=typer.colors.BRIGHT_RED,
                    )
                create_report_cmd, html_output = generate_create_report_cmd(
                    script_path,
                    markdown,
                    template_path,
                    cmo_patient_id,
                    small_variants_path,
                    manifest,
                    cnv_path,
                    dmp_patient_id,
                    dmp_sample_id,
                    facet_path,
                    tumor_type,
                )
                p3 = run_cmd(create_report_cmd, force)
                if "Error" in str(p3) or "error" in str(p3):
                    summary.append(
                        "\t".join(
                            [
                                cmo_patient_id,
                                dmp_patient_id,
                                str(dmp_sample_id),
                                "NA",
                                "create_report.R failed",
                            ]
                        )
                    )
                else:
                    summary.append(
                        "\t".join(
                            [
                                cmo_patient_id,
                                dmp_patient_id,
                                str(dmp_sample_id),
                                "NA",
                                "create_report.R ran without facet maf",
                            ]
                        )
                    )
                typer.secho(
                    f"Done running create_report.R for patient with CMO ID {cmo_patient_id} and output is written in {html_output}",
                    fg=typer.colors.BRIGHT_GREEN,
                )

    print("\nSummary for all patient processed..\n")
    print(summary)
    summary_file = Path.cwd().joinpath("run_create_report_summary.tsv")
    with open(summary_file, "w") as fp:
        for item in summary:
            # write each item on a new line
            fp.write("%s\n" % item)
    typer.secho("Done running run_create_report!", fg=typer.colors.BRIGHT_GREEN)


if __name__ == "__main__":
    typer.run(main)
