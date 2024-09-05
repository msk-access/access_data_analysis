import pandas as pd

recist_df = pd.read_csv('/juno/cmo/bergerlab/access_projects/EDD_R01/recist_raw.csv')
samples_df = pd.read_csv('/juno/cmo/bergerlab/access_projects/EDD_R01/all_samples.tsv', sep='\t')

merged_df = pd.merge(recist_df, samples_df, on='cmo_sample_id_plasma', how='left')

merged_df['treatment_date'] = pd.to_datetime(merged_df['treatment_date'], format="%m/%d/%y")
merged_df['exam_date'] = pd.to_datetime(merged_df['collection_date'], format="%m/%d/%y")

merged_df['exam_day'] = (merged_df['exam_date'] - merged_df['treatment_date']).dt.days

# re-map timepoint response column to RECIST values
timepoint_mapping = {1: "CR", 2: "PR", 3: "PD", 4: "SD", 5: "N/A", 6: "PR/SD"}
merged_df.replace({'Timepoint.Response': timepoint_mapping}, inplace=True)
merged_df.rename(columns={"dmp_id": "DMP.ID", "cmo_patient_id": "cmoPatientId"}, inplace=True)

cols = ["cohort", "cmoPatientId", "cmo_sample_id_plasma", "Timepoint.Response", "Total.Diameter", "exam_date", "exam_day", "DMP.ID"]
final_df = merged_df[cols]
final_df2 = final_df.dropna(subset=['Timepoint.Response'])

final_df2.to_csv('/juno/cmo/bergerlab/access_projects/EDD_R01/recist.csv', index=False)

# split recist data by cohort
for cohort, group in final_df2.groupby(['cohort']):
     cohort_string = str(cohort).strip("()").replace(",", "").replace("'", "")
     group.to_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort_string}_Mutations/VAFs/{cohort_string}_RECIST.csv', index=False)

