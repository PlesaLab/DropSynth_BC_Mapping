import pandas as pd
import sys
import glob

def combine_starcode_files(input_pattern, output_file):
    # Get list of input files
    file_names = glob.glob(input_pattern)

    if not file_names:
        print(f"No files found matching pattern: {input_pattern}")
        sys.exit(1)

    # Initialize an empty dataframe
    merged_df = pd.DataFrame()

    # Iterate over each file
    for file_name in file_names:
        df = pd.read_csv(file_name, sep='\t', header=None, names=['barcode', 'count'])
        merged_df = pd.concat([merged_df, df])

    # Group by barcode and sum the counts
    merged_df = merged_df.groupby('barcode').sum().reset_index()

    # Save the merged dataframe to a new TSV file
    merged_df.to_csv(output_file, sep='\t', index=False, header=False)

    print(f"Files merged successfully! Output saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python starcode_combine.py <input_file_pattern> <output_file>")
        sys.exit(1)

    input_pattern = sys.argv[1]
    output_file = sys.argv[2]
    combine_starcode_files(input_pattern, output_file)
