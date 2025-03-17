#!/usr/bin/env python3
import argparse
from Bio import SeqIO

def read_fasta(filename):
    """
    Reads a FASTA file and returns a list of tuples (header, sequence).
    """
    records = []
    for record in SeqIO.parse(filename, "fasta"):
        header = f">{record.id}"  # Ensure correct FASTA format
        sequence = str(record.seq)
        records.append((header, sequence))
    return records

def write_fasta(records, output_file):
    """
    Writes sorted sequences to a FASTA file.
    """
    with open(output_file, 'w') as out:
        for header, seq in records:
            out.write(f"{header}\n{seq}\n")

def main():
    parser = argparse.ArgumentParser(description="Sort a FASTA file by header.")
    parser.add_argument('-i', '--input', required=True, help="Input FASTA file")
    parser.add_argument('-o', '--output', required=True, help="Output sorted FASTA file")
    args = parser.parse_args()

    # Read and sort sequences
    records = read_fasta(args.input)
    records.sort(key=lambda x: x[0])  # Sort by header (alphabetical order)

    # Save sorted FASTA file
    write_fasta(records, args.output)

    print(f"✅ Sorted FASTA file saved to {args.output}")

if __name__ == '__main__':
    main()
