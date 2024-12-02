import sys
import pandas as pd
import csv
import argparse
import os

def checkPaths(bam_paths, output):

    paths_missing = False
    # open bam paths
    bam_paths = pd.read_csv(bam_paths)
    path_columns = ['bam_path_normal', 'bam_path_plasma_duplex', 'bam_path_plasma_simplex', 'maf_path', 'cna_path', 'sv_path', 'bam_index_normal', 'bam_index_plasma_duplex', 'bam_index_plasma_simplex']

    with open(output, 'w') as output:
        writer = csv.writer(output)
        for index,row in bam_paths.iterrows():
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check if BAM paths exist.")
    parser.add_argument("--bam_paths", required=True, help="Path to BAM paths CSV file.")
    parser.add_argument("--output", required=True, help="Path for output missing BAMS CSV file.")
    args = parser.parse_args()

    checkPaths(args.bam_paths, args.output)