# DropSynth Barcode Mapping Pipeline

This repository implements an automated pipeline to process high-throughput Nanopore sequencing data associated with barcoded DropSynth gene libraries. The workflow covers steps from splitting large FASTQ files and extracting barcodes to generating consensus gene sequences and determining mutant sequences. The Makefile ties together several Python scripts that each perform a dedicated task in the pipeline.

## Workflow Overview

The `Makefile` pipeline consists of the following main steps:

1. **FASTQ File Splitting** (*See split_script.py*)   
   Large FASTQ.gz files are split into smaller, more manageable parts using `split_script.py`.

2. **Barcode Extraction and Processing**  
   The `barcode_processing.py` script reads input FASTQ files, identifies specific barcode regions based on defined motifs, and extracts associated sequences. It outputs:
   - Barcode statistics (`.bc_stats.csv` and `.bc_stats_for_starcode.tsv`)
   - A list of barcode-associated sequences (`.bc_list.csv`)
   - A sorted FASTA file of sequences without ambiguous sites  
   *(See [barcode_processing.py](&#8203;:contentReference[oaicite:1]{index=1}))*

3. **Barcode Clustering (Starcode)**  
   An external tool (Starcode) is used to cluster similar barcodes. If multiple clustering output files are generated, they can be merged using `starcode_combine.py` to create a unified consensus barcode count file.  
   *(See [starcode_combine.py](&#8203;:contentReference[oaicite:2]{index=2}))*

4. **Assigning Consensus Barcodes**  
   The `add_consensus_bc.py` script takes a FASTA file and the Starcode clustering output. It assigns consensus barcodes to each record (updating the headers) and outputs a new FASTA file.  
   *(See [add_consensus_bc.py](&#8203;:contentReference[oaicite:3]{index=3}))*

5. **Sorting Consensus Sequences**  
   The resulting consensus FASTA file is sorted alphabetically by header using `sort_consensus_bc.py` to ensure proper ordering for downstream analysis.  
   *(See [sort_consensus_bc.py](&#8203;:contentReference[oaicite:4]{index=4}))*

6. **Consensus Gene Generation**  
   With a sorted FASTA file of consensus barcodes, `get_consensus_gene.py` groups sequences by barcode and aligns them using pyabpoa. For each group, a consensus gene is generated and scored, with results saved as:
   - A consensus gene FASTA file (headers correspond to consensus barcodes)
   - A text file listing consensus scores  
   *(See [get_consensus_gene.py](&#8203;:contentReference[oaicite:5]{index=5}))*

7. **Alternative Translation Processing**  
   Optionally, nucleotide sequences can be translated into protein sequences using `process_alt_trans.py`, which outputs the results as a CSV file.  
   *(See [process_alt_trans.py](&#8203;:contentReference[oaicite:6]{index=6}))*

8. **SAM File Parsing and Mutation Analysis**  
   Finally, the `parse_sam_script.py` script parses a SAM file (generated, for example, by BBMap) alongside reference protein sequences and barcode information. This step performs pairwise alignments to identify mutations, generate mutant IDs, and output two CSV reports:
   - Barcode-to-mutant mapping (`C5seqs_mutID_all.csv`)
   - Aggregated mutant information (`C5seqs_mutID_info_all.csv`)  
   *(See [parse_sam_script.py](&#8203;:contentReference[oaicite:7]{index=7}))*

## Dependencies

This pipeline uses the following dependencies, which can be installed via Conda using the provided `environment.yml` file:

- **Python 3.9**
- **biopython** – for sequence handling and FASTA/FASTQ parsing
- **pandas** – for data manipulation
- **numpy** – for numerical operations
- **regex** – for flexible motif matching
- **starcode 1.4** – for clustering similar barcodes
- **minimap2 2.28** – for read mapping (alternative to BBMap)
- **samtools 1.21** – for SAM/BAM file processing
- **bbmap 39.17** – for read mapping (default mapping tool)

Additionally, the following pip packages are installed:
- **edlib** – for sequence alignment
- **pyabpoa** – for multiple sequence alignment in consensus gene generation

To set up the environment, run:
```bash
conda env create -f environment.yml
conda activate newenv
```

## Makefile Targets and Running the Pipeline

The workflow is orchestrated via a Makefile. To facilitate running the entire pipeline in a high-performance computing environment, a shell script named `consensus_makefile.sh` is provided. This script performs the following steps:

1. Dynamically sets the project directory.

2. Loads necessary modules (e.g., `miniconda3`).

3. Activates the Conda environment (`newenv`).

4. Exports the thread count from the `SLURM_CPUS_PER_TASK` environment variable for OpenMP-based programs.

5. Runs the Makefile with the default mapping tool (BBMap).

The relevant portion of the `consensus_makefile.sh` script is:
```bash
# Set the project directory dynamically to the script's location
PROJ_DIR=$(pwd)
cd "$PROJ_DIR"

# Load necessary modules
module purge
module load miniconda3

# Activate the Conda environment
conda activate newenv

# Set number of threads for OpenMP-based programs
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Export thread count for use in Makefile
export THREAD_COUNT=$SLURM_CPUS_PER_TASK

###---- Leave as is if you want to use BBMap ---------------------
# Run the Makefile with MAPPING_TOOL=bbmap (default)
make -j $SLURM_CPUS_PER_TASK

###---- Uncomment if you want to use MiniMap instead of BBMap-----
# Run the Makefile with MAPPING_TOOL=minimap
#make -j $SLURM_CPUS_PER_TASK MAPPING_TOOL=minimap

# Run only the SAM parser using the MiniMap input
#make parse_sam_only_minimap
```

The `Makefile` itself includes targets for:

- **Splitting FASTQ files** (`split`)
- **Processing barcodes** (`process_barcodes`)
- **Assigning consensus barcodes** (`assign_consensus`)
- **Sorting consensus sequences** (`sort_consensus`)
- **Generating consensus gene sequences** (`consensus_gene`)
- **Parsing SAM files for mutation analysis** (`parse_sam`)
- **Alternative translation** (`alt_translation`)

To run the full pipeline, simply execute the `consensus_makefile.sh` script. In an HPC environment, submit this script as a batch job where `SLURM_CPUS_PER_TASK` defines the number of threads available.


