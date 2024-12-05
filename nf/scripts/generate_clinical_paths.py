import pandas as pd
import argparse
import re

def generateClinicalPaths(samples, id_mapping, dmp_dir, mirror_bam_dir, mirror_access_bam_dir, dmp_key_path, access_key_path, output):

    # get samples and patient ids for samples
    samples_df = pd.read_csv(samples)
    samples_df['cmo_patient_id'] = samples_df['cmo_sample_id'].str.split('-').str[:2].str.join('-')
    samples_df = samples_df[(samples_df['cmo_sample_id'].str[0] == "C")].copy()
    unique_patient_ids = samples_df["cmo_patient_id"].unique()
    ids_df = pd.read_csv(id_mapping)

    sample_sheet = []
    impact_calls = []

    for cmo_patient in unique_patient_ids:
        matching_rows = ids_df.loc[ids_df['cmo_patient_id'] == cmo_patient, 'dmp_patient_id']
        dmp_id = matching_rows.iloc[0]
        if matching_rows.isna().any():
            print(f'Patient {cmo_patient} has no DMP ID. Skipping.')
            continue

        with open(dmp_key_path, 'r') as impact_key:
            for line in impact_key:
                if re.search(rf"{dmp_id}-T.*-IM.*", line):
                    sample_id, standard_bam = getBamPath(line, mirror_bam_dir)
                    sample_sheet.append({'cmo_patient_id': cmo_patient, 'dmp_patient_id': dmp_id, 'sample_id': sample_id, 'standard_bam': standard_bam, 'duplex_bam': 'NA', 'simplex_bam': 'NA', 'maf': 'NA'})

                if re.search(rf"{dmp_id}-N.*-IM.*", line):
                    sample_id, standard_bam = getBamPath(line, mirror_bam_dir)
                    sample_sheet.append({'cmo_patient_id': cmo_patient, 'dmp_patient_id': dmp_id, 'sample_id': sample_id, 'standard_bam': standard_bam, 'duplex_bam': 'NA', 'simplex_bam': 'NA', 'maf': 'NA'})

        with open(access_key_path, 'r') as access_key:
            for line in access_key:
                if re.search(rf"{dmp_id}.*-T.*-XS.*simplex.*", line):
                    sample_id, simplex_bam = getBamPath(line, mirror_access_bam_dir)
                    sample_sheet.append({'cmo_patient_id': cmo_patient, 'dmp_patient_id': dmp_id, 'sample_id': sample_id, 'standard_bam': 'NA', 'duplex_bam': 'NA', 'simplex_bam': simplex_bam, 'maf': 'NA'})

                if re.search(rf"{dmp_id}.*-T.*-XS.*duplex.*", line):
                    sample_id, duplex_bam = getBamPath(line, mirror_access_bam_dir)
                    sample_sheet.append({'cmo_patient_id': cmo_patient, 'dmp_patient_id': dmp_id, 'sample_id': sample_id, 'standard_bam': 'NA', 'duplex_bam': duplex_bam, 'simplex_bam': 'NA', 'maf': 'NA'})

                if re.search(rf"{dmp_id}.*-N.*-XS.*unfilter.*", line):
                    sample_id, standard_bam = getBamPath(line, mirror_access_bam_dir)
                    sample_sheet.append({'cmo_patient_id': cmo_patient, 'dmp_patient_id': dmp_id, 'sample_id': sample_id, 'standard_bam': standard_bam, 'duplex_bam': 'NA', 'simplex_bam': 'NA', 'maf': 'NA'})
        
        sample_sheet_df = pd.DataFrame(sample_sheet)
        sample_sheet_df = (sample_sheet_df.groupby('sample_id', as_index=False)
                     .agg(lambda x: next((val for val in x if val != 'NA'), 'NA'))
)
        sample_sheet_df.to_csv(output, index=False)

        # combine with the research bam paths to make complete sample sheet
        # split by patient

        with open(f'{dmp_dir}/data_mutations_extended.txt', 'r') as impact:
            print(dmp_id)
            for line in impact:
                if re.search(rf"{dmp_id}", line) and not re.search("GERMINLINE", line) and not re.search('sequenced_samples:', line):
                    cols = line.split(sep='\t')

                    impact_calls.append({'Hugo_Symbol': cols[0], 
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
        
        # combine impact MAF with duplex MAF to make all_unique_calls.maf (split by patient)
        
    impact_calls_df = pd.DataFrame(impact_calls)
    impact_calls_df.to_csv("impact_calls.csv", index=False)

                    
                

    # for each patient
        # Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2, Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, t_ref_count, t_alt_count, n_ref_count, n_alt_count, Variant_Classification
        # generate a maf file (union of all calls)
        # add row in sample sheet with maf file

def getBamPath(line, mirror):

    cols = line.split(sep=',')
    sample_id = cols[0]
    base_bam = cols[1]
    fl = base_bam[0]
    sl = base_bam[1]
    bam_path = f'{mirror}/{fl}/{sl}/{base_bam}.bam'

    return sample_id, bam_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate clinical ACCESS and clinical IMPACT BAM paths.")
    parser.add_argument("--samples", required=True, help="Path to samples CSV file.")
    parser.add_argument("--id_mapping", required=True, help="Path to CMO ID - DMP ID mapping file.")
    parser.add_argument("--dmp_dir", required=True)
    parser.add_argument("--mirror_bam_dir", required=True)
    parser.add_argument("--mirror_access_bam_dir", required=True)
    parser.add_argument("--dmp_key_path", required=True)
    parser.add_argument("--access_key_path", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    generateClinicalPaths(args.samples, args.id_mapping, args.dmp_dir, args.mirror_bam_dir, args.mirror_access_bam_dir, args.dmp_key_path, args.access_key_path, args.output)

