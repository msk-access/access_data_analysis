import pandas as pd
from datetime import datetime

# calculate length of treatment for each patient
patients = pd.read_csv("/juno/cmo/bergerlab/access_projects/EDD_R01/all_patients.tsv", sep='\t')

patients['treatment_start'] = pd.to_datetime(patients['treatment_start'], format="%m/%d/%y")
patients['treatment_end'] = pd.to_datetime(patients['treatment_end'], format="%m/%d/%y")

patients['treatment_length'] = (patients['treatment_end'] - patients['treatment_start']).dt.days
patients.to_csv("/juno/cmo/bergerlab/access_projects/EDD_R01/all_patients_days.tsv", index=False, sep='\t')

# calculate day of collection for each sample
samples = pd.read_csv("/juno/cmo/bergerlab/access_projects/EDD_R01/all_samples.tsv", sep='\t')

samples['treatment_date'] = pd.to_datetime(samples['treatment_date'], format="%m/%d/%y")
samples['collection_date'] = pd.to_datetime(samples['collection_date'], format="%m/%d/%y")

samples['collection_day'] = (samples['collection_date'] - samples['treatment_date']).dt.days
samples.to_csv("/juno/cmo/bergerlab/access_projects/EDD_R01/all_samples_days.tsv", index=False, sep='\t')
