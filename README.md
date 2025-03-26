# DropSynth Barcode Mapping Pipeline

This repository implements an automated pipeline to process high-throughput Nanopore sequencing data associated with barcoded DropSynth gene libraries. The workflow covers steps from splitting large FASTQ files and extracting barcodes to generating consensus gene sequences and determining mutant sequences. The Makefile ties together several Python scripts that each perform a dedicated task in the pipeline. The default mapping 

## Workflow Overview

> [!NOTE]
> The `Makefile` pipeline consists of the following main steps:

1. **FASTQ File Splitting** (*split_script.py*)   
   Large FASTQ.gz files are split into smaller, more manageable parts using `split_script.py`.

2. **Barcode Extraction and Processing** (*barcode_processing.py*)   
   The `barcode_processing.py` script reads input FASTQ files, identifies specific barcode regions based on defined motifs, and extracts associated sequences. It outputs:
   - Barcode statistics (`.bc_stats.csv` and `.bc_stats_for_starcode.tsv`)
   - A list of barcode-associated sequences (`.bc_list.csv`)
   - A sorted FASTA file of sequences without ambiguous sites

3. **Barcode Clustering (Starcode)** (*starcode_combine.py*)   
   An external tool (Starcode) is used to cluster similar barcodes. If multiple clustering output files are generated, they can be merged using `starcode_combine.py` to create a unified consensus barcode count file.

4. **Assigning Consensus Barcodes** (*add_consensus_bc.py*)   
   The `add_consensus_bc.py` script takes a FASTA file and the Starcode clustering output. It assigns consensus barcodes to each record (updating the headers) and outputs a new FASTA file.  

5. **Sorting Consensus Sequences** (*sort_consensus_bc.py*)   
   The resulting consensus FASTA file is sorted alphabetically by header using `sort_consensus_bc.py` to ensure proper ordering for downstream analysis.

6. **Consensus Gene Generation** (*get_consensus_gene.py*)   
   With a sorted FASTA file of consensus barcodes, `get_consensus_gene.py` groups sequences by barcode and aligns them using pyabpoa. For each group, a consensus gene is generated and scored, with results saved as:
   - A consensus gene FASTA file (headers correspond to consensus barcodes)
   - A text file listing consensus scores

7. **Consensus Gene Translation** (*process_alt_trans.py*)   
   Translate the consensus gene sequences (nucleotides) into protein sequences (amino acids) using `process_alt_trans.py`, which outputs the results as a CSV file for downstream parsing.  

8. **Map Consensus Genes to Reference Genes** (*BBMap*)   
   Matches (maps) the consensus gene FASTA file (with headers corresponding to consensus barcodes) to the reference genes file (*.genes*) and generates a SAM file for downstream parsing. 
   - Default Mapping: **BBMAP**
   - Alternative Mapping: **MINIMAP**

9. **SAM File Parsing and Mutation Analysis** (*parse_sam_script.py*)   
   Finally, the `parse_sam_script.py` script parses a SAM file (default: **BBMap**) using the reference proteins file (*.proteins*) and barcode information. This step performs pairwise alignments to identify mutations, generate mutant IDs, and output two CSV reports:
   - Barcode-to-mutant mapping csv: `C5seqs_mutID_all.csv`
   - Aggregated mutant information csv: `C5seqs_mutID_info_all.csv`

## Dependencies

> [!NOTE]
> This pipeline uses the following dependencies, which can be installed via Conda using the provided `environment.yml` file:

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

> [!NOTE]
> The workflow is orchestrated via a `Makefile`. To facilitate running the entire pipeline in a high-performance computing environment, a shell script named `consensus_makefile.sh` is provided. This script performs the following steps:

1. Dynamically sets the project directory.

2. Loads necessary modules (e.g., `miniconda3`).

3. Activates the Conda environment (`newenv`).

4. Exports the thread count from the `SLURM_CPUS_PER_TASK` environment variable for OpenMP-based programs.

5. Runs the `Makefile` with the default mapping tool (BBMap).

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
> [!TIP]
> The `Makefile` itself includes targets for:

- **Splitting FASTQ files** (`split_script.py`)
- **Processing barcodes** (`barcode_processing.py`)
- **Clustering barcodes with Starcode** ('starcode_combine.py')
- **Assigning consensus barcodes** (`add_consensus_bc.py`)
- **Sorting consensus barcodes** (`sort_consensus_bc.py`)
- **Generating consensus gene sequences** (`get_consensus_gene.py`)
- **Translating consensus gene sequences** (`process_alt_trans.py`)
- **Parsing SAM files for perfect and mutant variants** (`parse_sam_script.py`)

To run the full pipeline, simply execute the `consensus_makefile.sh` script. In an HPC environment, submit this script as a batch job where `SLURM_CPUS_PER_TASK` defines the number of threads available.
```bash
sbatch consensus_makefile.sh
```
