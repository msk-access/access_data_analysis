import pandas as pd
import os
import argparse

def unlinkPaths(bam_paths, output):

    bam_paths = pd.read_csv(bam_paths)
    path_columns = ['bam_path_normal', 'bam_path_plasma_duplex', 'bam_path_plasma_simplex', 'maf_path', 'cna_path', 'sv_path', 'bam_index_normal', 'bam_index_plasma_duplex', 'bam_index_plasma_simplex']


    # Loop through each row and update the path columns with absolute paths
    for index, row in bam_paths.iterrows():
        for col in path_columns:
            # skip NA paths if normal is missing
            if pd.notna(row[col]):
                # Replace relative or softlinked path with absolute path
                absolute_path = os.path.realpath(row[col]) 
                bam_paths.at[index, col] = absolute_path

    # Write the updated DataFrame back to a new CSV file
    bam_paths.to_csv(output, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Replace BAM paths with absolute paths.")
    parser.add_argument("--bam_paths", required=True, help="Path to BAM paths CSV file.")
    parser.add_argument("--output", required=True, help="Path for output BAMS CSV file with absolute paths.")
    args = parser.parse_args()

    unlinkPaths(args.bam_paths, args.output)