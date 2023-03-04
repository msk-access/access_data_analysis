from pathlib import Path


def generate_create_report_cmd(
    script,
    markdown,
    template_file,
    cmo_patient_id,
    csv_file,
    manifest,
    cnv_path,
    dmp_patient_id,
    dmp_sample_id,
    dmp_facet_maf,
    tumor_type=None
):
    """Create the system command that should be run for create_report.R

    Args:
        script (str): path for create_report.R
        markdown (bool): True|False to generate markdown output
        template_file (str): path for the template file
        cmo_patient_id (str): patient id from CMO
        csv_file (str): path to csv file containing variant information
        tumor_type (str): tumor type label
        manifest (pathlib.Path): path to the manifest containing meta data
        cnv_path (pathlib.Path): path to directory having cnv files
        dmp_patient_id (str): patient id of the clinical msk-impact sample
        dmp_sample_id (str): sample id of the clinical msk-impact sample
        dmp_facet_maf (str): path to the clinical msk-impact maf file annotated for facets results

    Returns:
        str: system command to run for create_report.R
    """
    html_output = Path.cwd().joinpath(f"{cmo_patient_id}_report.html")
    cmd = (
            "Rscript "
            + str(script)
            + " -t "
            + str(template_file)
            + " -p "
            + str(cmo_patient_id)
            + " -r "
            + str(csv_file)
            + " -rc "
            + str(cnv_path.as_posix())
            + " -m "
            + str(manifest.as_posix())
            + " -o "
            + str(html_output)
            + " -d "
            + str(dmp_patient_id)
            + " -ds "
            + str(dmp_sample_id)
            + " -dm "
            + str(dmp_facet_maf)
        )
    if markdown:
        cmd = (
            f"{cmd} --md --tt {str(tumor_type)}"
            if tumor_type is not None
            else f"{cmd} --md"
        )
    return (cmd)
