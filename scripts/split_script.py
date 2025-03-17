import os
import sys
import subprocess
import argparse

def split_fastq_gz(input_file, output_prefix, target_size_mb=1000):
    print("Starting split_fastq_gz function...")  # DEBUG: Function start

    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_prefix)
    os.makedirs(output_dir, exist_ok=True)
    print(f"Ensured output directory exists: {output_dir}")  # DEBUG: Directory check

    # Step 1: Count the number of lines in the original .fastq.gz file
    print("Counting lines in the input file...")  # DEBUG: Step 1 start
    line_count_cmd = f"zcat {input_file} | wc -l"
    print(f"Executing command: {line_count_cmd}")  # DEBUG: Command
    line_count_result = subprocess.run(line_count_cmd, shell=True, capture_output=True, text=True)
    print("Command executed successfully.")  # DEBUG: Command success
    total_lines = int(line_count_result.stdout.strip())
    print(f"Total line count: {total_lines}")  # DEBUG: Line count
    print("Counting lines complete.")  # DEBUG: Step 1 end

    # Step 2: Calculate the number of lines per part
    print("Calculating lines per part...")  # DEBUG: Step 2 start
    file_size_bytes = os.path.getsize(input_file)
    bytes_per_line = file_size_bytes / total_lines
    target_lines = int((target_size_mb * 1024 * 1024) / bytes_per_line)
    lines_per_part = (target_lines // 4) * 4  # Ensure it's divisible by 4
    print(f"Lines per part: {lines_per_part}")  # DEBUG: Lines per part
    print("Calculation complete.")  # DEBUG: Step 2 end

    # Step 3: Split the file
    print("Splitting the file...")  # DEBUG: Step 3 start
    split_cmd = f"zcat {input_file} | split -d -l {lines_per_part} - {output_prefix}"
    print(f"Executing command: {split_cmd}")  # DEBUG: Command
    try:
        subprocess.run(split_cmd, shell=True, check=True)
        print("Splitting command executed successfully.")  # DEBUG: Command success
    except subprocess.CalledProcessError as e:
        print(f"Error during splitting: {e}")  # DEBUG: Error
        sys.exit(1)
    print("Splitting complete.")  # DEBUG: Step 3 end

    # Count the actual number of parts created
    actual_parts = len([f for f in os.listdir(output_dir) if f.startswith(os.path.basename(output_prefix))])
    print(f"Actual number of parts created: {actual_parts}")  # DEBUG: Actual parts

    # Step 4: Compress the resulting files
    print("Compressing the resulting files...")  # DEBUG: Step 4 start
    for i in range(actual_parts):
        part = f"{output_prefix}{str(i).zfill(2)}"
        if os.path.exists(part):  # Check if the file exists before compressing
            compress_cmd = f"gzip {part}"
            print(f"Executing command: {compress_cmd}")  # DEBUG: Command
            try:
                subprocess.run(compress_cmd, shell=True, check=True)
                print(f"Compression of {part} successful.")  # DEBUG: Compression success
            except subprocess.CalledProcessError as e:
                print(f"Error during compression of {part}: {e}")  # DEBUG: Error
                sys.exit(1)
    print("Compression complete.")  # DEBUG: Step 4 end

    print(f"Splitting complete. Files: {[f'{output_prefix}{str(i).zfill(2)}.gz' for i in range(actual_parts)]}")
    print("split_fastq_gz function complete.")  # DEBUG: Function end

if __name__ == "__main__":
    print("Starting main execution block...")  # DEBUG: Main start

    parser = argparse.ArgumentParser(description="Split a FASTQ.gz file into multiple parts.")
    parser.add_argument("input_file", help="Path to the input FASTQ.gz file")
    parser.add_argument("output_prefix", help="Prefix for the output files")
    parser.add_argument("--target_size_mb", type=int, default=1000, help="Target size of each split in megabytes (default: 1000)")
    
    args = parser.parse_args()
    print(f"Parsed arguments: {args}")  # DEBUG: Arguments

    print("Calling split_fastq_gz function with parsed arguments...")  # DEBUG: Function call
    split_fastq_gz(args.input_file, args.output_prefix, args.target_size_mb)
    print("split_fastq_gz function call complete.")  # DEBUG: Function call complete

    print("Main execution block complete.")  # DEBUG: Main end
