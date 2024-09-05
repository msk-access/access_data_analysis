import subprocess
import pandas as pd

# all cohorts
cohorts = ['RET', 'KRAS', 'HER2', 'BRAF']

# one cohort
#cohorts = ['BRAF']

for cohort in cohorts:
    # all patients
    data = pd.read_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/patient_list.tsv', sep='\t')
    # one patient
    #data = pd.read_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/patients/C-96JR4P_masterfile.tsv', sep='\t')

    for patient, group in data.groupby(['cmo_patient_id']):
        patient_string = str(patient).strip("()").replace(",", "").replace("'", "")

        filter_calls_cmd = f'Rscript /home/guturus1/repos/access_data_analysis/R/filter_calls.R -m /juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/patients/{patient_string}_masterfile.tsv -o /juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/results/{patient_string}_results'
        print(f"Executing command: {filter_calls_cmd}")
        subprocess.run(filter_calls_cmd, shell=True)


