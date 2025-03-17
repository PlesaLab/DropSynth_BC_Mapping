import regex
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import collections
import csv
import gzip
from concurrent.futures import ThreadPoolExecutor
import argparse

def process_record(record, start_motif, end_motif, motif):
    global count, match_count, count_no_startseq, count_no_endseq_found, rc_count, d, dnum, outputseqs_noN
    count += 1  # count read
    rc_flag = 0
    seq = str(record.seq).upper()
    match = regex.search(motif, seq, regex.BESTMATCH)

    if match is None:  # check other strand, maybe it's reversed?
        seq = str(record.seq.reverse_complement()).upper()
        match = regex.search(motif, seq, regex.BESTMATCH)
        rc_flag = 1
        rc_count += 1

    if match is not None:
        match_count += 1
        barcode = match.group(3)
        bc_index = seq.find(barcode)
        sequence = seq[0:bc_index]

        d[barcode].append(sequence)  # store as a dictionary key = barcode, value = list of sequences

        if barcode not in dnum:  # new barcode not seen before
            dnum[barcode] = 1
        else:  # possible collision
            dnum[barcode] = dnum[barcode] + 1

        startseq_index = sequence.find(start_motif)
        if startseq_index != -1:
            endseq_index = sequence.find(end_motif)
            if endseq_index != -1:
                seq_between_start_end_motif = sequence[startseq_index + len(start_motif):endseq_index]
                if 'N' not in seq_between_start_end_motif:
                    if seq_between_start_end_motif != "":
                        outputseqs_noN.append(SeqRecord(id=barcode, seq=Seq(seq_between_start_end_motif), description=""))
            else:
                count_no_endseq_found += 1
        else:
            count_no_startseq += 1

def process_fastq(input_file, output_prefix, max_threads=80, start_motif="CATATG", end_motif="TAAGGTACCTAA"):
    global count, match_count, count_no_startseq, count_no_endseq_found, rc_count, d, dnum, outputseqs_noN
    
    count = 0
    match_count = 0
    count_no_startseq = 0
    count_no_endseq_found = 0
    rc_count = 0
    d = collections.defaultdict(list)
    dnum = dict()
    outputseqs_noN = []

    motif = r'((CTGCCGAACAGC)(....................)(AGGTGAAGAGCC)){e<3}'  # e<4 = less than 4 errors #sequences flanking the barcode

    with gzip.open(input_file, "rt") as handle:
        records = list(SeqIO.parse(handle, "fastq"))

    print(f"{len(records)} total records")

    with ThreadPoolExecutor(max_threads) as executor:
        executor.map(lambda record: process_record(record, start_motif, end_motif, motif), records)

    print(f"{rc_count} times tried RC {len(records)}")
    print(f"{len(d)} BCs found out of {len(records)}")
    print(f"{count_no_startseq} no start (NdeI) site out of {match_count}")
    print(f"{count_no_endseq_found} no End site found (but has start) out of {match_count - count_no_startseq}")
    print(f"{len(outputseqs_noN)} no ambiguous sites out of {match_count - count_no_startseq - count_no_endseq_found}")

    # Write output files
    with open(f"{output_prefix}.bc_stats.csv", "w") as outfile:
        writer = csv.writer(outfile)
        for k, v in dnum.items():
            writer.writerow([k, str(v)])

    with open(f"{output_prefix}.bc_stats_for_starcode.tsv", "w") as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        for k, v in dnum.items():
            writer.writerow([k, str(v)])

    with open(f"{output_prefix}.bc_list.csv", "w") as outfile:
        writer = csv.writer(outfile)
        for k, v in d.items():
            for ss in v:
                writer.writerow([k, str(ss)])

    with open(f"{output_prefix}.sorted.fasta", 'w') as output_file_x:
        SeqIO.write(outputseqs_noN, output_file_x, "fasta")

    print(f"{len(dnum)} write out bc_stats.csv")
    print(f"{sum(len(v) for v in d.values())} write out bc_list.csv")
    print(f"{len(outputseqs_noN)} write out fasta no N")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process FASTQ file for barcode analysis.")
    parser.add_argument("input_file", help="Path to the input FASTQ.gz file")
    parser.add_argument("output_prefix", help="Prefix for the output files")
    parser.add_argument("--max_threads", type=int, default=80, help="Maximum number of threads to use")
    parser.add_argument("--start_motif", default="CATATG", help="Start motif sequence")
    parser.add_argument("--end_motif", default="TAAGGTACCTAA", help="End motif sequence")
    
    args = parser.parse_args()
    
    process_fastq(args.input_file, args.output_prefix, args.max_threads, args.start_motif, args.end_motif)
