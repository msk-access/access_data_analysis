import pandas as pd

# all cohorts
cohorts = ['RET', 'KRAS', 'HER2', 'BRAF']

# one cohort
#cohorts = ['HER2']

for cohort in cohorts:
    # all patients
    patient_data = pd.read_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/patient_list.tsv', sep='\t')
    # one patient
    #patient_data = pd.read_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/patients/C-96JR4P_masterfile.tsv', sep='\t')
    vafs = []

    for patient, group in patient_data.groupby(['cmo_patient_id']):
        patient_string = str(patient).strip("()").replace(",", "").replace("'", "")

        if patient_string in ('C-JP56U3', 'C-4KE2A7', 'C-6J30R4', 'C-217F4D', 'C-5NYJV2', 'C-911P8A', 'C-XP28CV', 'C-HAT9MW', 'C-YVTHXU'):
            print(patient_string)
            print("skipped")
            continue

        variants = f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/results/{patient_string}_results/{patient_string}_variant_table.csv'
        data = pd.read_csv(variants)
        # Filter for passing variants
        passing_variants = data[data['QC'] == 'PASS']

        # Only keep columns for research and clinical ACCESS samples
        #sample_columns = [col for col in data.columns if not col.str.contains('IM')]
        sample_columns = [col for col in data.columns if col.startswith('C-') or 'XS' in col]

        # Initialize an empty DataFrame to store the results
        results = []

        # Iterate through each passing variant
        for idx, row in passing_variants.iterrows():

            var_name = f"{row['Hugo_Symbol']}{'' if pd.isna(row['HGVSp_Short']) else row['HGVSp_Short']}"

            for col in sample_columns:
                sample_name = col

                # Extract VAF value
                vaf_value = row[col].split('(')[-1].replace(')', '') if '(' in row[col] else row[col]
                vaf = float(vaf_value)

                if row['clonality'] == 'CLONAL':
                    t = row['tcn']
                    exp = row['expected_alt_copies']
                    ncn = 2

                    # adj_vaf = (n* vaf_value) / (exp + (n - t))
                    vaf = vaf*ncn / (exp + (ncn - t)*vaf )

                results.append({'cmoSampleName': sample_name, 'VarName': var_name, 'VAF': vaf, 'cmoPatientId': patient_string})


        result_df = pd.DataFrame(results)
        vafs.append(result_df)
        # Save the result to a new CSV file
        result_df.to_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/VAFs/{patient_string}_VAFs.csv', index=False)

        print(f"The table for passing variants has been saved to {patient_string}_VAFs.csv.")

    # Concatenate all DataFrames in the list into a single DataFrame
    combined_vafs = pd.concat(vafs, ignore_index=True)
    samples_df = pd.read_csv('/juno/cmo/bergerlab/access_projects/EDD_R01/data/all_samples_deid.tsv', sep='\t')
    patients_df = pd.read_csv('/juno/cmo/bergerlab/access_projects/EDD_R01/data/all_patients_deid.tsv', sep='\t')
    samples_df.rename(columns={"cmo_sample_id_plasma": "cmoSampleName", "dmp_patient_id": 'DMP.ID'}, inplace=True)
    patients_df.rename(columns={"cmo_sample_id_plasma": "cmoSampleName", "reason_tx_stop": "reason_for_tx_stop"}, inplace=True)

    first_merge = pd.merge(samples_df, patients_df, on='cmo_patient_id', how='left')
    second_merge = pd.merge(combined_vafs, first_merge, on='cmoSampleName', how='left')

    cols = ["cmoSampleName", "VarName", "VAF", "cmoPatientId", "DMP.ID", "collection_day", "treatment_length", "treatment_ongoing", "reason_for_tx_stop"]
    final_df = second_merge[cols]
    average_df = final_df.groupby('cmoSampleName')['VAF'].mean().reset_index()
    average_df['VAF'] = average_df['VAF'] * 10
    
    # Save the combined DataFrame to a new CSV file
    final_df.to_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/VAFs/{cohort}_VAFs.csv', index=False)
    average_df.to_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/VAFs/{cohort}_mean_VAFs.csv', index=False)
    print(f"The table for combined vafs has been saved to EDD_R01/{cohort}_Mutations/VAFs/{cohort}_VAFs.csv'.")




