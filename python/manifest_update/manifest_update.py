"""
Author: Carmelina Charalambous
Date: June 21 2024
Description: Fills in ACCESS manifest file with specific paths ready for access_data_analysis
"""

import pandas as pd
import argparse

def update_manifest(input_file_path, output_file_name):

    # Load input  manifest
    df = pd.read_excel(input_file_path)

    # add base paths for normal/duplex/simplex bams and snv/cnv/sv files
    base_path_bams = "/work/access/production/data/bams/"
    base_path_small_variants = "/work/access/production/data/small_variants/"
    base_path_copy_number_variants = "/work/access/production/data/copy_number_variants/"
    base_path_structural_variants = "/work/access/production/data/structural_variants/"

    # Fill in bam normal
    df['bam_path_normal'] = base_path_bams + df['cmo_patient_id'] + '/' + df['cmo_sample_id_normal'] + '/current/' + df['cmo_sample_id_normal'] + '_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX.bam'

    # Fill in bam duplex
    df['bam_path_plasma_duplex'] = base_path_bams + df['cmo_patient_id']  + '/' + df['cmo_sample_id_plasma'] + '/current/' + df['cmo_sample_id_plasma'] + '_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.bam'

    #  Fill in bam  simplex
    df['bam_path_plasma_simplex'] = base_path_bams + df['cmo_patient_id'] + '/'  + df['cmo_sample_id_plasma'] + '/current/' + df['cmo_sample_id_plasma'] + '_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-simplex.bam'

    # Fill maf
    df['maf_path'] = base_path_small_variants + df['cmo_patient_id'] +  '/' + df['cmo_sample_id_plasma'] + '/current/' + df['cmo_sample_id_plasma'] + '.DONOR22-TP.combined-variants.vep_keptrmv_taggedHotspots_fillout_filtered.maf'

    # Fill in cna
    df['cna_path'] = base_path_copy_number_variants + df['cmo_patient_id'] + '/' + df['cmo_sample_id_plasma'] + '/current/' + df['cmo_sample_id_plasma'] + '_copynumber_segclusp.genes.txt'

    # Fill in sv
    df['sv_path'] = base_path_structural_variants + df['cmo_patient_id'] + '/'  + df['cmo_sample_id_plasma'] + '/current/'  + df['cmo_sample_id_plasma'] + '_AllAnnotatedSVs.txt'

    # Construct output file paths
    output_file_path_xlsx = f"{output_file_name}.xlsx"
    output_file_path_csv = f"{output_file_name}.csv"

    # Save the updated manifest to both Excel and CSV formats
    df.to_excel(output_file_path_xlsx, index=False)
    df.to_csv(output_file_path_csv, index=False)

    print(f'Saved updated manifest as {output_file_path_xlsx} and {output_file_path_csv}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Update manifest file with specified paths.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to the input manifest file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Prefix name for the output files (without extension)')

    args = parser.parse_args()

    update_manifest(args.input, args.output)