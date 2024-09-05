import pandas as pd

volume_df = pd.read_csv('/juno/cmo/bergerlab/access_projects/EDD_R01/data/volumetric_raw.tsv', sep='\t')
patients_df = pd.read_csv('/juno/cmo/bergerlab/access_projects/EDD_R01/data/all_patients.tsv', sep='\t')

merged_df = pd.merge(volume_df, patients_df, on='MRN', how='left')

merged_df['treatment_start'] = pd.to_datetime(merged_df['treatment_start'], format="%m/%d/%y")
merged_df['Study_Date'] = pd.to_datetime(merged_df['Study_Date'], format="%m/%d/%y")

merged_df['Study_Day'] = (merged_df['Study_Date'] - merged_df['treatment_start']).dt.days

clean_df = merged_df[['cohort', 'cmo_patient_id', 'dmp_patient_id', 'MRN',  'treatment_start', 'Study_Date', 'Study_Day',
                      'Target_Lesion', 'Lesion_Site', 'Uni_mm', 'Perp_mm', 'Bi_mm^2', 'Volume_mm^3', 'cancer_type']]

clean_df.to_csv('/juno/cmo/bergerlab/access_projects/EDD_R01/data/volumetric.tsv', sep='\t', index=False)

# split data by cohort
for cohort, group in clean_df.groupby(['cohort']):
     cohort_string = str(cohort).strip("()").replace(",", "").replace("'", "")
     group.to_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort_string}_Mutations/VAFs/{cohort_string}_volumetric.csv', index=False)

