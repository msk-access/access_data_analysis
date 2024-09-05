import subprocess
import pandas as pd

# all cohorts
cohorts = ['RET', 'KRAS', 'HER2', 'BRAF']

# one cohort
# cohorts = ['BRAF']

for cohort in cohorts:

    # all patients
    data = pd.read_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/patient_list.tsv', sep='\t')

    # one patient
    # data = pd.read_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/patients/C-YVTHXU_masterfile.tsv', sep='\t')

    for patient, group in data.groupby(['cmo_patient_id']):
        patient_string = str(patient).strip("()").replace(",", "").replace("'", "")
        dmp_id = group['dmp_patient_id'].iloc[0]

        if pd.isna(dmp_id):
            compile_reads_all_cmd = f'bsub -J {patient_string}_compile_reads -W 12:00 -n 1 -R rusage[mem=4] -e /juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/results/{patient_string}_results/%J.err -o /juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/results/{patient_string}_results/%J.out Rscript /home/guturus1/repos/access_data_analysis/R/compile_reads.R -m /juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/patients/{patient_string}_masterfile.tsv -o /juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/results/{patient_string}_results'

        else:
            compile_reads_all_cmd = f'bsub -J {patient_string}_compile_reads -W 12:00 -n 1 -R rusage[mem=4] -e /juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/results/{patient_string}_results/%J.err -o /juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/results/{patient_string}_results/%J.out Rscript /home/guturus1/repos/access_data_analysis/R/compile_reads_all.R -m /juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/patients/{patient_string}_masterfile.tsv -o /juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/results/{patient_string}_results'
        
        print(f"Executing command: {compile_reads_all_cmd}")
        subprocess.run(compile_reads_all_cmd, shell=True)

