import pandas as pd

cohorts = ['HER2', 'BRAF', 'KRAS', 'RET']

for cohort in cohorts:

     mean_vaf_df = pd.read_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/VAFs/{cohort}_mean_VAFs.csv')
     samples_df = pd.read_csv('/juno/cmo/bergerlab/access_projects/EDD_R01/data/all_samples_deid.tsv', sep='\t')
     samples_df.rename(columns={"cmo_sample_id_plasma": "cmoSampleName"}, inplace=True)
     mean_vaf_df = mean_vaf_df.merge(samples_df[['cmoSampleName','cmo_patient_id', 'collection_day', 'dmp_patient_id']],on='cmoSampleName', how='left')

     patients_df = pd.read_csv('/juno/cmo/bergerlab/access_projects/EDD_R01/data/all_patients_deid.tsv', sep='\t')
     mean_vaf_df = mean_vaf_df.merge(patients_df[['cmo_patient_id','treatment_length', 'reason_tx_stop']],on='cmo_patient_id', how='left')
     
     mean_vaf_df.rename(columns={'dmp_patient_id': 'DMP.ID', 'cmo_patient_id': 'cmoPatientId', 'reason_tx_stop': 'reason_for_tx_stop'}, inplace=True)

     final_df = mean_vaf_df[['cmoSampleName', 'cmoPatientId', 'collection_day', 'DMP.ID', 'treatment_length', 'reason_for_tx_stop', 'VAF']]
     final_df.to_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/VAFs/{cohort}_volumetric_vaf_table.tsv', sep='\t', index=False)
