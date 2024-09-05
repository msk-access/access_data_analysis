import subprocess
import pandas as pd
import numpy as np

cohorts = ['RET', 'KRAS', 'HER2', 'BRAF']

for cohort in cohorts:
     # all samples
    sample_data = pd.read_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/sample_list.tsv', sep='\t')
    # all patients
    patient_data = pd.read_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/patient_list.tsv', sep='\t')
    dedup_patient_data = patient_data.drop(columns=['cohort', 'project', 'sex', 'dmp_patient_id', 'MRN', 'cmo_sample_id_normal', 'paired'])

    # create a data file with all the treatment information necessary from patient_data and sample_data
    merged_df = pd.merge(sample_data, dedup_patient_data, on="cmo_patient_id")

    merged_df = merged_df.rename(columns={'cmo_patient_id': 'cmoPatientId', 'cmo_sample_id_plasma': 'sample_id', 'treatment_date': 'start',
                                          'treatment_end': 'endtouse', 'reason_tx_stop': 'reason'})
    
    merged_df['assay_type'] = np.where(merged_df['sample_id'].str.contains('-IM') , 'IMPACT', 'ACCESS')
    merged_df['clinical_or_research'] = np.where(merged_df['sample_id'].str.startswith('P-') , 'Clinical', 'Research')

    merged_df.to_csv(f"/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/swimmer_info.tsv", sep='\t', index=False)

    #swimmer_cmd = f'Rscript /juno/cmo/bergerlab/access_projects/EDD_R01/scripts/swimmer_plots.R -i /juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/swimmer_info.csv -o /juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/visualizations/{cohort}_swimmer_plot.pdf'

    #print(f"Executing command: {swimmer_cmd}")
    #subprocess.run(swimmer_cmd, shell=True)

    # delete the input file