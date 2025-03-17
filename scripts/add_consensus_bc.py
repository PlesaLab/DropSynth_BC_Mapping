#!/usr/bin/env python3
import argparse
from Bio import SeqIO

def read_starcode_clusters(clustered_fp):
    """
    Reads the Starcode output file and maps each original barcode to its consensus barcode.
    """
    bc_map = {}  # Dictionary mapping original barcode -> consensus barcode
    with open(clustered_fp, 'r') as f:
        for line in f:
            consensus_bc, count, collapsed_bcs = line.strip().split('\t')
            collapsed_bcs = collapsed_bcs.split(',')
            for bc in collapsed_bcs:
                bc_map[bc] = consensus_bc
    return bc_map

def assign_consensus_barcodes(input_fasta, bc_map):
    """
    Reads a FASTA file, replaces barcodes with their consensus barcode, 
    and outputs a new FASTA file with updated headers.
    """
    updated_records = []
    
    for record in SeqIO.parse(input_fasta, "fasta"):
        original_bc = record.id.split('.')[0]  # Extract barcode from header
        consensus_bc = bc_map.get(original_bc, original_bc)  # Assign consensus barcode if available
        record.id = consensus_bc
        record.description = ""  # Ensure no extra description
        updated_records.append(record)

    return updated_records

def save_fasta(records, output_fasta):
    """
    Saves updated barcode-sequence records to a FASTA file.
    """
    with open(output_fasta, "w") as out_f:
        SeqIO.write(records, out_f, "fasta")

    print(f"✅ Consensus barcodes assigned. Output saved to {output_fasta}")

def main():
    parser = argparse.ArgumentParser(description="Assign consensus barcodes to sequences and save as FASTA.")
    parser.add_argument('-i', dest='input_fasta', help="Input FASTA file", required=True)
    parser.add_argument('-c', dest='clustered_fp', help="Starcode output file of clustered barcodes", required=True)
    parser.add_argument('-o', dest='out_fasta', help="Output FASTA file", required=True)
    args = parser.parse_args()

    # Read Starcode clusters and process sequences
    bc_map = read_starcode_clusters(args.clustered_fp)
    updated_records = assign_consensus_barcodes(args.input_fasta, bc_map)

    if updated_records:
        save_fasta(updated_records, args.out_fasta)

if __name__ == "__main__":
    main()
