import pandas as pd
import argparse
import csv
import os
import sys
import re
import numpy as np

def generateResearchPaths(samples,
                          bam_path_normal,
                          bam_path_plasma_duplex,
                          bam_path_plasma_simplex,
                          bam_path_index_normal,
                          bam_path_index_duplex,
                          bam_path_index_simplex,
                          maf_path,
                          cna_path,
                          sv_path):
    
    # open samples.csv and parse patient ids
    samples_df = pd.read_csv(samples)
    samples_df['cmo_patient_id'] = samples_df['cmo_sample_id'].str.split('-').str[:2].str.join('-')
    
    # create a new dataframe for storing bam paths
    bam_paths_df = []

    # create df with rows that are tumor only and research access only
    tumor_df = samples_df[(samples_df['tumor_normal'] == "tumor") & (samples_df['cmo_sample_id'].str[0] == "C")].copy()
    normal_df = samples_df[samples_df['tumor_normal'] == "normal"].copy()

    # merge the tumor_df with normal_df on the cmo_patient_id column
    normal_df.rename(columns={'cmo_sample_id': 'cmo_sample_id_normal'}, inplace=True)
    tumor_df = tumor_df.merge(normal_df[['cmo_patient_id', 'cmo_sample_id_normal']], 
                           on='cmo_patient_id', 
                           how='left')

    # get normal, duplex, simplex, maf, cna, sv paths
    tumor_df['bam_path_normal'] = tumor_df.apply(lambda row: bam_path_normal.format(cmo_patient_id=row['cmo_patient_id'], cmo_sample_id_normal=row['cmo_sample_id_normal']), axis=1)
    tumor_df['bam_path_plasma_duplex'] = tumor_df.apply(lambda row: bam_path_plasma_duplex.format(cmo_patient_id=row['cmo_patient_id'], cmo_sample_id=row['cmo_sample_id']), axis=1)
    tumor_df['bam_path_plasma_simplex'] = tumor_df.apply(lambda row: bam_path_plasma_simplex.format(cmo_patient_id=row['cmo_patient_id'], cmo_sample_id=row['cmo_sample_id']), axis=1)
    
    tumor_df['bam_path_index_normal'] = tumor_df.apply(lambda row: bam_path_index_normal.format(cmo_patient_id=row['cmo_patient_id'], cmo_sample_id_normal=row['cmo_sample_id_normal']), axis=1)
    tumor_df['bam_path_index_plasma_duplex'] = tumor_df.apply(lambda row: bam_path_index_duplex.format(cmo_patient_id=row['cmo_patient_id'], cmo_sample_id=row['cmo_sample_id']), axis=1)
    tumor_df['bam_path_index_plasma_simplex'] = tumor_df.apply(lambda row: bam_path_index_simplex.format(cmo_patient_id=row['cmo_patient_id'], cmo_sample_id=row['cmo_sample_id']), axis=1)

    tumor_df['maf_path'] = tumor_df.apply(lambda row: maf_path.format(cmo_patient_id=row['cmo_patient_id'], cmo_sample_id=row['cmo_sample_id']), axis=1)
    tumor_df['cna_path'] = tumor_df.apply(lambda row: cna_path.format(cmo_patient_id=row['cmo_patient_id'], cmo_sample_id=row['cmo_sample_id']), axis=1)
    tumor_df['sv_path'] = tumor_df.apply(lambda row: sv_path.format(cmo_patient_id=row['cmo_patient_id'], cmo_sample_id=row['cmo_sample_id']), axis=1)

    # output to csv bam paths
    for patient, group in tumor_df.groupby(['cmo_patient_id']):
        patient_string = str(patient).strip("()").replace(",", "").replace("'", "")
        group.to_csv(f'./{patient_string}_masterfile.csv', index=False)
    
    checkPaths(tumor_df)
    return(tumor_df)

def checkPaths(bam_paths_df):

    paths_missing = False
    # open bam paths
    # bam_paths = pd.read_csv(bam_paths)
    path_columns = [col for col in bam_paths_df.columns if 'path' in col]
    real_bam_paths_df = bam_paths_df.copy()
    print(path_columns)
    #columns_containing_string = [col for col in df.columns if df[col].astype(str).str.contains(string_to_search, ignore_case=True).any()]

    with open("./check_paths.csv", 'w') as output:
        writer = csv.writer(output)
        for index,row in bam_paths_df.iterrows():
            sample = row['cmo_sample_id']
            for col in path_columns:
                if col == 'bam_path_normal' and pd.isna(row[col]):
                    continue
                if col == 'bam_index_normal' and pd.isna(row[col]):
                    continue
                path = row[col].strip()
                # check if path exists
                if not os.path.isfile(path):
                    # write to csv if missing
                    writer.writerow([sample, col, path])
                    paths_missing = True

    if paths_missing:
        print("BAM file contains missing paths.")
        sys.exit(1)
    else:
        print("BAM paths all exist.")

    for index, row in real_bam_paths_df.iterrows():
        for col in path_columns:
            # skip NA paths if normal is missing
            if pd.notna(row[col]):
                # Replace relative or softlinked path with absolute path
                absolute_path = os.path.realpath(row[col]) 
                real_bam_paths_df.at[index, col] = absolute_path

    # Write the updated DataFrame back to a new CSV file
    real_bam_paths_df.to_csv("./real_bam_paths.csv", index=False)

def generateClinicalPaths(samples,
                            id_mapping, 
                            dmp_dir, 
                            mirror_bam_dir,
                            mirror_access_bam_dir,
                            dmp_key_path,
                            access_key_path):

    # get samples and patient ids for samples
    samples_df = pd.read_csv(samples)
    samples_df['cmo_patient_id'] = samples_df['cmo_sample_id'].str.split('-').str[:2].str.join('-')
    samples_df = samples_df[(samples_df['cmo_sample_id'].str[0] == "C")].copy()
    unique_patient_ids = samples_df["cmo_patient_id"].unique()
    ids_df = pd.read_csv(id_mapping)

    sample_sheet = []

    for cmo_patient in unique_patient_ids:
        # get the DMP ID from the mapping
        matching_rows = ids_df.loc[ids_df['cmo_patient_id'] == cmo_patient, 'dmp_patient_id']
        dmp_id = matching_rows.iloc[0]
        if matching_rows.isna().any():
            print(f'Patient {cmo_patient} has no DMP ID. Skipping.')
            continue

        with open(dmp_key_path, 'r') as impact_key:
            for line in impact_key:
                if re.search(rf"{dmp_id}-T.*-IM.*", line):
                    sample_id, standard_bam = getClinicalBam(line, mirror_bam_dir)
                    sample_sheet.append({'cmo_patient_id': cmo_patient, 'dmp_patient_id': dmp_id, 'sample_id': sample_id, 'standard_bam': standard_bam, 'duplex_bam': 'NA', 'simplex_bam': 'NA', 'maf': 'NA'})

                if re.search(rf"{dmp_id}-N.*-IM.*", line):
                    sample_id, standard_bam = getClinicalBam(line, mirror_bam_dir)
                    sample_sheet.append({'cmo_patient_id': cmo_patient, 'dmp_patient_id': dmp_id, 'sample_id': sample_id, 'standard_bam': standard_bam, 'duplex_bam': 'NA', 'simplex_bam': 'NA', 'maf': 'NA'})

        with open(access_key_path, 'r') as access_key:
            for line in access_key:
                if re.search(rf"{dmp_id}.*-T.*-XS.*simplex.*", line):
                    sample_id, simplex_bam = getClinicalBam(line, mirror_access_bam_dir)
                    sample_sheet.append({'cmo_patient_id': cmo_patient, 'dmp_patient_id': dmp_id, 'sample_id': sample_id, 'standard_bam': 'NA', 'duplex_bam': 'NA', 'simplex_bam': simplex_bam, 'maf': 'NA'})

                if re.search(rf"{dmp_id}.*-T.*-XS.*duplex.*", line):
                    sample_id, duplex_bam = getClinicalBam(line, mirror_access_bam_dir)
                    sample_sheet.append({'cmo_patient_id': cmo_patient, 'dmp_patient_id': dmp_id, 'sample_id': sample_id, 'standard_bam': 'NA', 'duplex_bam': duplex_bam, 'simplex_bam': 'NA', 'maf': 'NA'})

                if re.search(rf"{dmp_id}.*-N.*-XS.*unfilter.*", line):
                    sample_id, standard_bam = getClinicalBam(line, mirror_access_bam_dir)
                    sample_sheet.append({'cmo_patient_id': cmo_patient, 'dmp_patient_id': dmp_id, 'sample_id': sample_id, 'standard_bam': standard_bam, 'duplex_bam': 'NA', 'simplex_bam': 'NA', 'maf': 'NA'})
    
    sample_sheet_df = pd.DataFrame(sample_sheet)
    sample_sheet_df = (sample_sheet_df.groupby('sample_id', as_index=False)
                     .agg(lambda x: next((val for val in x if val != 'NA'), 'NA')))
    sample_sheet_df.to_csv("sample_sheet.csv", index=False)
    return(sample_sheet_df)

def getClinicalBam(line, mirror):

    cols = line.split(sep=',')
    sample_id = cols[0]
    base_bam = cols[1]
    fl = base_bam[0]
    sl = base_bam[1]
    bam_path = f'{mirror}/{fl}/{sl}/{base_bam}.bam'

    return sample_id, bam_path

def getPlasmaCalls(tumor_df):
    plasma_calls = []
    for path in tumor_df["maf_path"]:
        with open(path, 'r') as maf:
            reader = csv.DictReader(maf, delimiter='\t')
            for row in reader:
                if not row['Status']:
                    continue
                plasma_calls.append({'Hugo_Symbol': row['Hugo_Symbol'], 
                    'Chromosome': row['Chromosome'], 
                    'Start_Position': row['Start_Position'], 
                    'End_Position': row['End_Position'], 
                    'Reference_Allele': row['Reference_Allele'], 
                    'Tumor_Seq_Allele1': row['Tumor_Seq_Allele1'], 
                    'Tumor_Seq_Allele2': row['Tumor_Seq_Allele2'],
                    'Tumor_Sample_Barcode': row['Tumor_Sample_Barcode'],
                    'Matched_Norm_Sample_Barcode': row['Matched_Norm_Sample_Barcode'],
                    't_ref_count': row['D_t_ref_count_fragment'],
                    't_alt_count': row['D_t_alt_count_fragment'],
                    'n_ref_count': row['n_ref_count_fragment'],
                    'n_alt_count': row['n_alt_count_fragment'],
                    'Variant_Classification': row['Variant_Classification']})
                #plasma_calls.append({"patient_id": patient, "sample_id": sample_id, "hugo_symbol": hugo, "site": site, "fold_change": fc, "p_val": p_val})

    plasma_calls_df = pd.DataFrame(plasma_calls)
    #plasma_calls_df["CNA_tumor"] = np.where(plasma_calls_df["fold_change"] > "0", "AMP", "HOMDEL")
    plasma_calls_df.to_csv("plasma_calls.csv", index=False)
    return(plasma_calls_df)

def getDMPCalls(samples,
                id_mapping,
                dmp_dir):

    # get samples and patient ids for samples
    samples_df = pd.read_csv(samples)
    samples_df['cmo_patient_id'] = samples_df['cmo_sample_id'].str.split('-').str[:2].str.join('-')
    samples_df = samples_df[(samples_df['cmo_sample_id'].str[0] == "C")].copy()
    unique_patient_ids = samples_df["cmo_patient_id"].unique()
    ids_df = pd.read_csv(id_mapping)

    dmp_calls = []
    duplex_plasma_calls = []
    dmp_SV_calls = []
    plasma_SV_calls = []
    dmp_CNA_calls = []
    plasma_CNA_calls = []

    for cmo_patient in unique_patient_ids:
        # get the DMP ID from the mapping
        matching_rows = ids_df.loc[ids_df['cmo_patient_id'] == cmo_patient, 'dmp_patient_id']
        dmp_id = matching_rows.iloc[0]
        if matching_rows.isna().any():
            print(f'Patient {cmo_patient} has no DMP ID. Skipping.')
            continue

        with open(f'{dmp_dir}/data_mutations_extended.txt', 'r') as dmp:
            print(dmp_id)
            for line in dmp:
                if re.search(rf"{dmp_id}", line) and not re.search("GERMLINE", line) and not re.search('sequenced_samples:', line):
                    cols = line.split(sep='\t')

                    dmp_calls.append({'Hugo_Symbol': cols[0], 
                                            'Chromosome': cols[4], 
                                            'Start_Position': cols[5], 
                                            'End_Position': cols[6], 
                                            'Reference_Allele': cols[11], 
                                            'Tumor_Seq_Allele1': cols[12], 
                                            'Tumor_Seq_Allele2': cols[13],
                                            'Tumor_Sample_Barcode': cols[16],
                                            'Matched_Norm_Sample_Barcode': 'NA',
                                            't_ref_count': 0,
                                            't_alt_count': 0,
                                            'n_ref_count': 0,
                                            'n_alt_count': 0,
                                            'Variant_Classification': cols[9]})
        

                
    dmp_calls_df = pd.DataFrame(dmp_calls)
    dmp_calls_df.to_csv("dmp_calls.csv", index=False)

    return(dmp_calls_df)

def getDMPCNACalls(samples, id_mapping, dmp_dir):
    
    # get samples and patient ids for samples
    samples_df = pd.read_csv(samples)
    samples_df['cmo_patient_id'] = samples_df['cmo_sample_id'].str.split('-').str[:2].str.join('-')
    samples_df = samples_df[(samples_df['cmo_sample_id'].str[0] == "C")].copy()
    unique_patient_ids = samples_df["cmo_patient_id"].unique()
    ids_df = pd.read_csv(id_mapping)

    dmp_CNA_calls = []

    for cmo_patient in unique_patient_ids:
        # get the DMP ID from the mapping
        matching_rows = ids_df.loc[ids_df['cmo_patient_id'] == cmo_patient, 'dmp_patient_id']
        dmp_id = matching_rows.iloc[0]
        if matching_rows.isna().any():
            print(f'Patient {cmo_patient} has no DMP ID. Skipping.')
            continue

        with open(f'{dmp_dir}/data_CNA.txt', "r") as dmp_CNA:
            reader = csv.reader(dmp_CNA, delimiter="\t")
            dmp_data = list(reader)

            sample_cols = dmp_data[0]
            
            for index, sample_col in enumerate(sample_cols[1:], start=1):
                if dmp_id in sample_col:
                    for row in dmp_data[1:]:
                        fc = row[index]
                        if fc and fc!= '0':
                            print(fc)
                            hugo_symbol = row[0]
                            dmp_CNA_calls.append({"dmp_patient_id": dmp_id, "dmp_sample_id": sample_col, "hugo_symbol": hugo_symbol, "fold_change": fc})

    dmp_CNA_calls_df = pd.DataFrame(dmp_CNA_calls)
    dmp_CNA_calls_df["fold_change"] = pd.to_numeric(dmp_CNA_calls_df["fold_change"])
    dmp_CNA_calls_df["CNA_tumor"] = np.where(dmp_CNA_calls_df["fold_change"] > 0, "AMP", "HOMDEL")
    dmp_CNA_calls_df.to_csv("dmp_CNA_calls.csv", index=False)

def getPlasmaCNACalls(tumor_df):

    plasma_CNA_calls = []
    for path in tumor_df["cna_path"]:
        with open(path, 'r') as cna_data:
            reader = csv.DictReader(cna_data, delimiter='\t')
            for row in reader:
                if row["p.adj"] < "0.05":
                    sample_id = row["sample"]
                    hugo = row["region"]
                    fc = row['fc']
                    p_val = row["p.adj"]
                    site = row["Cyt"]
                    patient = '-'.join(row['sample'].split('-')[:2])
                    plasma_CNA_calls.append({"patient_id": patient, "sample_id": sample_id, "hugo_symbol": hugo, "site": site, "fold_change": fc, "p_val": p_val})

    plasma_CNA_calls_df = pd.DataFrame(plasma_CNA_calls)
    plasma_CNA_calls_df["CNA_tumor"] = np.where(plasma_CNA_calls_df["fold_change"] > "0", "AMP", "HOMDEL")
    plasma_CNA_calls_df.to_csv("plasma_CNA_calls.csv", index=False)

    
    #for patient, group in tumor_df.groupby(['cmo_patient_id']):
     #   patient_string = str(patient).strip("()").replace(",", "").replace("'", "")
      #  group.to_csv(f'./{patient_string}_CNA.csv', index=False)

def getDMPSVCalls(samples, id_mapping, dmp_dir):
    # get samples and patient ids for samples
    samples_df = pd.read_csv(samples)
    samples_df['cmo_patient_id'] = samples_df['cmo_sample_id'].str.split('-').str[:2].str.join('-')
    samples_df = samples_df[(samples_df['cmo_sample_id'].str[0] == "C")].copy()
    unique_patient_ids = samples_df["cmo_patient_id"].unique()
    ids_df = pd.read_csv(id_mapping)

    dmp_SV_calls = []

    for cmo_patient in unique_patient_ids:
        # get the DMP ID from the mapping
        matching_rows = ids_df.loc[ids_df['cmo_patient_id'] == cmo_patient, 'dmp_patient_id']
        dmp_id = matching_rows.iloc[0]
        print(dmp_id)
        if matching_rows.isna().any():
            print(f'Patient {cmo_patient} has no DMP ID. Skipping.')
            continue

        with open(f'{dmp_dir}/data_sv.txt', "r") as dmp_SV:
            reader = csv.DictReader(dmp_SV, delimiter="\t")
            dmp_data = list(reader)
            
            for row in dmp_data:
                if dmp_id in row["Sample_ID"]:
                    print(row)
                    dmp_SV_calls.append(row)
    dmp_SV_calls_df = pd.DataFrame(dmp_SV_calls)
    dmp_SV_calls_df.to_csv("dmp_SV_calls.csv", index=False)
                        
def getPlasmaSVCalls(tumor_df):

    plasma_SV_calls = []
    for path in tumor_df["sv_path"]:
        with open(path, 'r') as sv_data:
            reader = csv.DictReader(sv_data, delimiter='\t')
            for row in reader:
                plasma_SV_calls.append(row)
    plasma_SV_calls_df = pd.DataFrame(plasma_SV_calls)
    plasma_SV_calls_df.to_csv("plasma_SV_calls.csv", index=False)

def getMSIScores(samples, id_mapping):
    # get samples and patient ids for samples
    samples_df = pd.read_csv(samples)
    samples_df['cmo_patient_id'] = samples_df['cmo_sample_id'].str.split('-').str[:2].str.join('-')
    samples_df = samples_df[(samples_df['cmo_sample_id'].str[0] == "C")].copy()
    unique_patient_ids = samples_df["cmo_patient_id"].unique()
    ids_df = pd.read_csv(id_mapping)

    MSI_scores = []

    for cmo_patient in unique_patient_ids:
        # get the DMP ID from the mapping
        matching_rows = ids_df.loc[ids_df['cmo_patient_id'] == cmo_patient, 'dmp_patient_id']
        dmp_id = matching_rows.iloc[0]
        print(dmp_id)
        if matching_rows.isna().any():
            print(f'Patient {cmo_patient} has no DMP ID. Skipping.')
            continue

        with open('/juno/work/access/production/resources/dmp_msi_admie_scores/current/all_admie_results_from_database.csv', "r") as dmp_MSI:
            reader = csv.DictReader(dmp_MSI)
            dmp_data = list(reader)
            
            for row in dmp_data:
                if dmp_id in row["DMP_PATIENT_ID"]:
                    print(row)
                    MSI_scores.append(row)
    MSI_scores_df = pd.DataFrame(MSI_scores)
    MSI_scores_df.to_csv("MSI_scores.csv", index=False)

def aggregateSampleSheet(research_df, clinical_df, maf_paths):
    
    # remove some columns in the research
    research_df = research_df[['cmo_patient_id', "bam_path_plasma_duplex", "bam_path_plasma_simplex"]].copy()
    research_df.rename(columns={'bam_path_plasma_duplex': 'duplex_bam', 'bam_path_plasma_simplex': 'simplex_bam'}, inplace=True)
    research_df.insert(1, 'standard_bam', 'NA')
    research_df['maf'] = 'NA'

    # keep cmo_patient_id, standard bam is na, duplex bam is bam_path_plasma_duplex, simplex is bam_path_plasma simplex 
    # keep the clincal as is, or remove the dmp_id column 
    clinical_df.drop('dmp_patient_id', axis=1, inplace=True)
    clinical_df.drop('sample_id', axis=1, inplace=True)

    combined_sample_sheet = pd.concat([research_df, clinical_df], ignore_index=True)
    for patient, group in combined_sample_sheet.groupby(['cmo_patient_id']):
        patient_string = str(patient).strip("()").replace(",", "").replace("'", "")
        print(f'Processing {patient_string} sample sheet')
        maf_path = maf_paths[patient_string]
        group['maf'] =  os.path.abspath(maf_path)
        group.to_csv(f'./{patient_string}_sample_sheet.csv', index=False)
    
    return(combined_sample_sheet)

def generateMAFs(plasma_calls, dmp_calls):
    all_small_calls = pd.concat([plasma_calls, dmp_calls], ignore_index=True)
    all_small_calls['cmo_patient_id'] = all_small_calls['Tumor_Sample_Barcode'].str.split('-').str[:2].str.join('-')
    all_small_calls = all_small_calls.drop_duplicates(keep='first')
    mafs = {}
    

    for patient, group in all_small_calls.groupby(['cmo_patient_id']):
        patient_string = str(patient).strip("()").replace(",", "").replace("'", "")
        print(f'Processing {patient_string} sample sheet')
        maf_path = f'./{patient_string}_all_calls.maf'
        group.to_csv(maf_path, index=False)
        mafs[patient_string] = maf_path

    return(mafs)





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate BAM paths.")
    parser.add_argument("--samples", required=True, help="Path to samples CSV file.")
    parser.add_argument("--id_mapping", required=True)
    parser.add_argument("--bam_path_normal", required=True)
    parser.add_argument("--bam_path_plasma_duplex", required=True)
    parser.add_argument("--bam_path_plasma_simplex", required=True)
    parser.add_argument("--bam_path_index_normal", required=True)
    parser.add_argument("--bam_path_index_duplex", required=True)
    parser.add_argument("--bam_path_index_simplex", required=True)
    parser.add_argument("--maf_path", required=True)
    parser.add_argument("--cna_path", required=True)
    parser.add_argument("--sv_path", required=True)
    parser.add_argument("--dmp_dir", required=True)
    parser.add_argument("--mirror_bam_dir", required=True)
    parser.add_argument("--mirror_access_bam_dir", required=True)
    parser.add_argument("--dmp_key_path", required=True)
    parser.add_argument("--access_key_path", required=True)
    args = parser.parse_args()

    
    masterfile = generateResearchPaths(args.samples, args.bam_path_normal, args.bam_path_plasma_duplex, args.bam_path_plasma_simplex, args.bam_path_index_normal, args.bam_path_index_duplex, args.bam_path_index_simplex, args.maf_path, args.cna_path, args.sv_path)
    clincal_sample_sheet = generateClinicalPaths(args.samples, args.id_mapping, args.dmp_dir, args.mirror_bam_dir, args.mirror_access_bam_dir, args.dmp_key_path, args.access_key_path)
    #getPlasmaCNACalls(masterfile)
    #getPlasmaSVCalls(masterfile)
    #getDMPSVCalls(args.samples, args.id_mapping, args.dmp_dir)
    #getMSIScores(args.samples, args.id_mapping)
    #getDMPCNACalls(args.samples, args.id_mapping, args.dmp_dir)

    plasma_calls = getPlasmaCalls(masterfile)
    DMP_calls = getDMPCalls(args.samples, args.id_mapping, args.dmp_dir)
    maf_paths = generateMAFs(plasma_calls, DMP_calls)
    sample_sheet = aggregateSampleSheet(masterfile, clincal_sample_sheet, maf_paths)