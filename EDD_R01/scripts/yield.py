import pandas as pd

samples_df = pd.read_csv('/juno/cmo/bergerlab/access_projects/EDD_R01/all_samples.tsv', sep='\t')
yield_df = pd.read_csv('/juno/cmo/bergerlab/access_projects/EDD_R01/yield.csv')

merged_df = pd.merge(samples_df, yield_df, on='cmo_sample_id_plasma', how='left')
merged_df.to_csv('/juno/cmo/bergerlab/access_projects/EDD_R01/merged_yield.tsv', sep='\t', index=False)


