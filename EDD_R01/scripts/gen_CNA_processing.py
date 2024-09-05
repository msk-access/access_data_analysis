import subprocess
import pandas as pd

cohorts = ['RET', 'KRAS', 'HER2', 'BRAF']
#cohorts = ['BRAF']

for cohort in cohorts:
    data = pd.read_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/patient_list.tsv', sep='\t')
    #data = pd.read_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/patients/C-96JR4P_masterfile.tsv', sep='\t')

    for patient, group in data.groupby(['cmo_patient_id']):
        patient_string = str(patient).strip("()").replace(",", "").replace("'", "")
    
        cna_cmd = f'Rscript /home/guturus1/repos/access_data_analysis/R/CNA_processing.R -m /juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/patients/{patient_string}_masterfile.tsv -o /juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/results/{patient_string}_results'
        print(f"Executing command: {cna_cmd}")
        subprocess.run(cna_cmd, shell=True)


