import pandas as pd
import argparse

def generatePaths(samples, output):
    # open samples.csv
    samples_df = pd.read_csv(samples)
    print(output)
    samples_df['cmo_patient_id'] = samples_df['cmo_sample_id'].str.split('-').str[:2].str.join('-')
    # create a new dataframe for storing bam paths
    bam_paths_df = []

    # create df with rows that are tumor only and research access only
    normal_df = samples_df[samples_df['tumor_normal'] == "normal"].copy()
    tumor_df = samples_df[(samples_df['tumor_normal'] == "tumor") & (samples_df['cmo_sample_id'].str[0] == "C")].copy()

    # Merge the tumor_df with normal_df on the cmo_patient_id column
    normal_df.rename(columns={'cmo_sample_id': 'cmo_sample_id_normal'}, inplace=True)
    tumor_df = tumor_df.merge(normal_df[['cmo_patient_id', 'cmo_sample_id_normal']], 
                           on='cmo_patient_id', 
                           how='left')

    # normal, duplex, simplex, maf, cna, sv
    base_path_bams = "/juno/work/access/production/data/bams/"
    base_path_small_variants = "/juno/work/access/production/data/small_variants/"
    base_path_copy_number_variants = "/juno/work/access/production/data/copy_number_variants/"
    base_path_structural_variants = "/juno/work/access/production/data/structural_variants/"

    tumor_df['bam_path_normal'] = base_path_bams + tumor_df['cmo_patient_id'] + '/' + tumor_df['cmo_sample_id_normal'] + '/current/' + tumor_df['cmo_sample_id_normal'] + '_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX.bam'
    tumor_df['bam_path_plasma_duplex'] = base_path_bams + tumor_df['cmo_patient_id']  + '/' + tumor_df['cmo_sample_id'] + '/current/' + tumor_df['cmo_sample_id'] + '_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.bam'
    tumor_df['bam_path_plasma_simplex'] = base_path_bams + tumor_df['cmo_patient_id'] + '/'  + tumor_df['cmo_sample_id'] + '/current/' + tumor_df['cmo_sample_id'] + '_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-simplex.bam'
    
    tumor_df['bam_index_normal'] = base_path_bams + tumor_df['cmo_patient_id'] + '/' + tumor_df['cmo_sample_id_normal'] + '/current/' + tumor_df['cmo_sample_id_normal'] + '_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX.bai'
    tumor_df['bam_index_plasma_duplex'] = base_path_bams + tumor_df['cmo_patient_id']  + '/' + tumor_df['cmo_sample_id'] + '/current/' + tumor_df['cmo_sample_id'] + '_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.bai'
    tumor_df['bam_index_plasma_simplex'] = base_path_bams + tumor_df['cmo_patient_id'] + '/'  + tumor_df['cmo_sample_id'] + '/current/' + tumor_df['cmo_sample_id'] + '_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-simplex.bai'

    tumor_df['maf_path'] = base_path_small_variants + tumor_df['cmo_patient_id'] +  '/' + tumor_df['cmo_sample_id'] + '/current/' + tumor_df['cmo_sample_id'] + '.DONOR22-TP.combined-variants.vep_keptrmv_taggedHotspots_fillout_filtered.maf'
    tumor_df['cna_path'] = base_path_copy_number_variants + tumor_df['cmo_patient_id'] + '/' + tumor_df['cmo_sample_id'] + '/current/' + tumor_df['cmo_sample_id'] + '_copynumber_segclusp.genes.txt'
    tumor_df['sv_path'] = base_path_structural_variants + tumor_df['cmo_patient_id'] + '/'  + tumor_df['cmo_sample_id'] + '/current/'  + tumor_df['cmo_sample_id'] + '_AllAnnotatedSVs.txt'

    # get absolute paths and not softlinked paths

    # output to csv bam paths
    tumor_df.to_csv(output, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate BAM paths.")
    parser.add_argument("--samples", required=True, help="Path to samples CSV file.")
    parser.add_argument("--output", required=True, help="Path for output BAM CSV file.")
    args = parser.parse_args()

    generatePaths(args.samples, args.output)