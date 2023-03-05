import typer
import pandas as pd
from pathlib import Path
from modules.run_cmd import run_cmd
from modules.read_manifest import read_manifest
from modules.check_required_columns import check_required_columns
from modules.generate_repo_paths import generate_repo_path
from modules.get_small_variant_csv import get_small_variant_csv
from modules.generate_facet_maf_path import generate_facet_maf_path
from modules.generate_create_report_cmd import generate_create_report_cmd


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
        help="File containing meta information per sample.",
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
    tumor_type: str = typer.Option(
        "",
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
):
    # Read the manifest file
    manifest_df = read_manifest(manifest)
    # check required columns
    column_header = check_required_columns(manifest_df, template_days)
    # get general paths
    (script_path, template_path) = generate_repo_path(
        repo_path, script_path, template_path, template_days
    )
    # iterate through each row and select information needed to generate the command
    for i in range(len(manifest_df)):
        cmo_patient_id = manifest_df.loc[i, "cmo_patient_id"]
        dmp_patient_id = manifest_df.loc[i, "dmp_patient_id"]
        typer.secho(
            f"Running for patient with CMO ID {cmo_patient_id}, and DMP ID {dmp_patient_id}",
            fg=typer.colors.BRIGHT_GREEN,
        )
        small_variants_path = get_small_variant_csv(cmo_patient_id, variant_path)
        if "dmp_sample_id" in column_header:
            dmp_sample_id = manifest_df.loc[i, "dmp_sample_id"]
            facet_path = generate_facet_maf_path(
                facet_repo, dmp_patient_id, dmp_sample_id
            )
        else:
            facet_path = generate_facet_maf_path(facet_repo, dmp_patient_id, None)
            # Get the sample id from the Facet file
            facet_path = Path(facet_path)
            maf_id = facet_path.stem
            dmp_sample_id = maf_id.split("_", 1)[0]
        if copy_facet:
            if not copy_facet_dir:
                copy_facet_dir = Path.cwd() / 'facet_files'
                copy_facet_dir.mkdir(parents=True, exist_ok=True)
            cp_facet_cmd = (f"cp {facet_path} {copy_facet_dir.as_posix()}")
            typer.secho(
                f"command: {cp_facet_cmd}",
                fg=typer.colors.BRIGHT_MAGENTA,
                )
            p1 = run_cmd(cp_facet_cmd)
            typer.secho(
            f"Done copying facet maf file for patient with CMO ID {cmo_patient_id}, and DMP ID {dmp_patient_id} and output is written in {copy_facet_dir}",
            fg=typer.colors.BRIGHT_MAGENTA,
            )
        create_report_cmd,html_output = generate_create_report_cmd(
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
        typer.secho(
            f"command: {create_report_cmd}",
            fg=typer.colors.BRIGHT_MAGENTA,
        )
        p2 = run_cmd(create_report_cmd)
        typer.secho(
            f"Done running create_report.R for patient with CMO ID {cmo_patient_id}, and DMP ID {dmp_patient_id} and output is written in {html_output}",
            fg=typer.colors.BRIGHT_MAGENTA,
        )

if __name__ == "__main__":
    typer.run(main)
