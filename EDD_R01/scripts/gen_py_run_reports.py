import subprocess
import pandas as pd

# options 'BRAF', 'KRAS', 'HER2', 'RET'
cohort = 'KRAS'
# options 'initial' or 'final'
stage = 'final'

data = pd.read_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/patient_list.tsv', sep='\t')

report_path = f'/home/guturus1/repos/access_data_analysis/reports/template_{stage}_report_run.Rmd'
#report_path = f'/home/guturus1/repos/access_data_analysis/reports/template_{stage}_report_run_yield.Rmd'


for patient, group in data.groupby(['cmo_patient_id']):
    patient_string = str(patient).strip("()").replace(",", "").replace("'", "")
    cancer_type = group['cancer_type'].iloc[0]
    #dmp_id = group['dmp_patient_id'].iloc[0]

    report_cmd = f'python /home/guturus1/repos/access_data_analysis/python/run_create_report/run_create_report.py \
                -m /juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/patients/clinical/{patient_string}_masterfile.tsv \
                -s /home/guturus1/repos/access_data_analysis/reports/create_report.R \
                -t {report_path} \
                -v /juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/results/{patient_string}_results/results_stringent/ \
                -c /juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/results/{patient_string}_results/CNA_final_call_set/ \
                -l {cancer_type} \
                -d \
                -ff'

    print(f"Executing command: {report_cmd}")
    subprocess.run(report_cmd, shell=True)
