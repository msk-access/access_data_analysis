import pandas as pd
import os

patients_df = pd.read_csv("/juno/cmo/bergerlab/access_projects/EDD_R01/all_patients_deid.tsv", sep='\t')
samples_df = pd.read_csv("/juno/cmo/bergerlab/access_projects/EDD_R01/all_samples_deid.tsv", sep='\t')
cohorts = []

# list of patients in each cohort 
for cohort, group in patients_df.groupby(['cohort']):
     cohort_string = str(cohort).strip("()").replace(",", "").replace("'", "")
     cohorts.append(cohort_string)
     group.to_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort_string}_Mutations/masterfiles/patient_list.tsv', index=False, sep='\t')

# list of samples in each cohort
for cohort, group in samples_df.groupby(['cohort']):
     cohort_string = str(cohort).strip("()").replace(",", "").replace("'", "")
     group.to_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort_string}_Mutations/masterfiles/sample_list.tsv', index=False, sep='\t')

# generate bam paths

for cohort in cohorts:
     df = pd.read_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/sample_list.tsv', sep='\t')
     df = df[df['project'] != 'clinical_access']
     df = df[df['project'] != 'clinical_impact']
     base_path_bams = "/juno/work/access/production/data/bams/"
     base_path_small_variants = "/juno/work/access/production/data/small_variants/"
     base_path_copy_number_variants = "/juno/work/access/production/data/copy_number_variants/"
     base_path_structural_variants = "/juno/work/access/production/data/structural_variants/"
     df['bam_path_normal'] = base_path_bams + df['cmo_patient_id'] + '/' + df['cmo_sample_id_normal'] + '/current/' + df['cmo_sample_id_normal'] + '_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX.bam'
     df['bam_path_plasma_duplex'] = base_path_bams + df['cmo_patient_id']  + '/' + df['cmo_sample_id_plasma'] + '/current/' + df['cmo_sample_id_plasma'] + '_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.bam'
     df['bam_path_plasma_simplex'] = base_path_bams + df['cmo_patient_id'] + '/'  + df['cmo_sample_id_plasma'] + '/current/' + df['cmo_sample_id_plasma'] + '_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-simplex.bam'
     df['maf_path'] = base_path_small_variants + df['cmo_patient_id'] +  '/' + df['cmo_sample_id_plasma'] + '/current/' + df['cmo_sample_id_plasma'] + '.DONOR22-TP.combined-variants.vep_keptrmv_taggedHotspots_fillout_filtered.maf'
     df['cna_path'] = base_path_copy_number_variants + df['cmo_patient_id'] + '/' + df['cmo_sample_id_plasma'] + '/current/' + df['cmo_sample_id_plasma'] + '_copynumber_segclusp.genes.txt'
     df['sv_path'] = base_path_structural_variants + df['cmo_patient_id'] + '/'  + df['cmo_sample_id_plasma'] + '/current/'  + df['cmo_sample_id_plasma'] + '_AllAnnotatedSVs.txt'
     df.to_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/masterfile.tsv', index=False, sep='\t')
     
     paths = pd.read_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/masterfile.tsv', sep='\t')
     path_columns = ['bam_path_normal', 'bam_path_plasma_duplex', 'bam_path_plasma_simplex', 'maf_path', 'cna_path', 'sv_path']
     
     with open(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/check_paths.txt', 'w') as output:
        for index,row in paths.iterrows():
             sample = row['cmo_sample_id_plasma']
             for col in path_columns:
                 if col == 'bam_path_normal' and pd.isna(row[col]):
                     continue
                 path = row[col].strip()
                 if not os.path.isfile(path):
                     output.write(f"{sample} {col} does not exist\n")
 

for cohort in cohorts:

    # list of samples per patient with clinical 
    c = pd.read_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/sample_list.tsv', sep='\t')
    # change column name from 'cmo_sample_id_plasma' to 'cmo_sample_id' for report generation
    c.rename(columns={'cmo_sample_id_plasma':'cmo_sample_id'}, inplace=True)
    for patient, group in c.groupby(['cmo_patient_id']):
        patient_string = str(patient).strip("()").replace(",", "").replace("'", "")
        group.to_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/patients/clinical/{patient_string}_masterfile.tsv', index=False, sep='\t')


    # masterfiles with paths per patient
    p = pd.read_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/masterfile.tsv', sep='\t')
    for patient, group in p.groupby(['cmo_patient_id']):
        patient_string = str(patient).strip("()").replace(",", "").replace("'", "")
        group.to_csv(f'/juno/cmo/bergerlab/access_projects/EDD_R01/{cohort}_Mutations/masterfiles/patients/{patient_string}_masterfile.tsv', index=False, sep='\t')