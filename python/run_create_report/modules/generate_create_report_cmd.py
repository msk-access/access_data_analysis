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
    tumor_type,
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
        cmd (str): system command to run for create_report.R
        html_output (pathlib.Path): where the output file should be written
    """
    cohort = 'KRAS'
    stage = 'final'
    html_output = f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/results/{cmo_patient_id}_results/{cmo_patient_id}_{stage}_report.html'
    if dmp_facet_maf:
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
            + " -tt "
            + f"'{tumor_type}'"
        )
    else:
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
            + " -tt "
            + f"'{tumor_type}'"
        )
    if markdown:
        cmd = f"{cmd} -md"
    return (cmd, html_output)
