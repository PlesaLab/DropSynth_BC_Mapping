#!/usr/bin/env python3
import argparse
import multiprocessing as mp
import pyabpoa as pa
import random

def calculate_consensus_score(aligned_result):
    """
    Calculate a consensus score from pyabpoa alignment result.

    Parameters:
    - aligned_result: msa_result object from pyabpoa.msa_aligner.msa()

    Returns:
    - float: Normalized consensus score between 0 and 1.
    """
    consensus_coverage = aligned_result.cons_cov
    n_seq = aligned_result.n_seq

    # flatten coverage list if nested
    if isinstance(consensus_coverage[0], list):
        consensus_coverage = [item for sublist in consensus_coverage for item in sublist]

    normalized_coverages = [cov / n_seq for cov in consensus_coverage]
    score = sum(normalized_coverages) / len(normalized_coverages)
    return score

def align_seqs(chunk):
    """
    Align a chunk of sequences using the abPOA aligner.
    Limits the number of sequences to 200.
    Returns a tuple with the barcode, number of sequences, consensus sequence, and score.
    """
    barcode, seqs = chunk
    seqs_len = len(seqs)
    if seqs_len == 1:
        return barcode, 1, seqs[0], -1  # score -1 for a single sequence
    elif seqs_len > 200:
        seqs = random.sample(seqs, 200)
        seqs_len = 200
    alignment = pa.msa_aligner().msa(seqs, out_cons=True, out_msa=False)
    score = calculate_consensus_score(alignment)
    consensus = alignment.cons_seq[0]
    return barcode, seqs_len, consensus, score

def fasta_records(file_handle):
    """
    Generator that yields FASTA records as tuples: (header, sequence)
    Assumes that header lines start with '>'.
    """
    header = None
    seq_lines = []
    for line in file_handle:
        line = line.rstrip('\n')
        if line.startswith('>'):
            if header is not None:
                yield header, "".join(seq_lines)
            header = line[1:].strip()  # remove '>' and extra whitespace
            seq_lines = []
        else:
            seq_lines.append(line.strip())
    if header is not None:
        yield header, "".join(seq_lines)

def get_seqs(file_handle):
    """
    Generator that groups FASTA records by consensus barcode.
    Assumes the FASTA file is sorted by header.
    
    Yields:
      tuple: (barcode, [list of sequences])
    """
    current_barcode = None
    current_chunk = []
    
    for header, seq in fasta_records(file_handle):
        # Use the first token from the header as the barcode.
        barcode = header.split()[0]
        if current_barcode is None:
            current_barcode = barcode
        if barcode != current_barcode:
            yield current_barcode, current_chunk
            current_barcode = barcode
            current_chunk = []
        current_chunk.append(seq)
    
    if current_chunk:
        yield current_barcode, current_chunk

def main():
    parser = argparse.ArgumentParser(
        description="Group sequences by consensus barcode and compute consensus gene using pyabpoa"
    )
    parser.add_argument('-i', dest='input_fp', help="Input sorted consensus barcode FASTA file", required=True)
    parser.add_argument('-o', dest='out_fp', help="Output FASTA file with consensus gene", required=True)
    parser.add_argument('-s', dest='score_fp', help="Output text file with consensus scores", required=True)
    parser.add_argument('-p', dest='threads', default=12, help="Threads", type=int)
    args = parser.parse_args()

    with open(args.input_fp, 'r') as f_in, \
         open(args.out_fp, 'w') as f_out, \
         open(args.score_fp, 'w') as f_score, \
         mp.Pool(processes=args.threads) as pool:

        for result in pool.imap_unordered(align_seqs, get_seqs(f_in), chunksize=5000):
            barcode, num_seqs, consensus, score = result
            # Write the consensus gene record with header as the consensus barcode only.
            f_out.write(f">{barcode}\n{consensus}\n")
            f_score.write(f"{barcode}\t{score}\n")

if __name__ == "__main__":
    main()
