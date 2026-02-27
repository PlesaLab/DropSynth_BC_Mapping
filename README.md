# DropSynth Barcode Mapping Pipeline

This repository implements an automated analysis pipeline for high-throughput Nanopore sequencing data generated from barcoded DropSynth gene libraries. The workflow spans preprocessing of large FASTQ files, barcode extraction and clustering, consensus gene reconstruction, translation, reference mapping, and mutation analysis.

The pipeline is orchestrated through a `Makefile` that links modular Python scripts, each performing a dedicated stage in the analysis. The default mapping tool is **BBMap**, with optional support for **Minimap2**.

This repository is intended for research use only.

---

## 🔍 Workflow Overview

> **NOTE**  
> The `Makefile` pipeline consists of the following main steps:

### 1. FASTQ File Splitting (`split_script.py`)
Large `FASTQ.gz` files are split into smaller, more manageable chunks for parallel processing.

### 2. Barcode Extraction and Processing (`barcode_processing.py`)
Reads input FASTQ files, identifies barcode regions based on defined motifs, and extracts associated sequences.

Outputs:
- Barcode statistics (`.bc_stats.csv`, `.bc_stats_for_starcode.tsv`)
- Barcode-associated read list (`.bc_list.csv`)
- Sorted FASTA file of sequences without ambiguous bases

### 3. Barcode Clustering (Starcode) (`starcode_combine.py`)
Barcodes are clustered using **Starcode**. If multiple clustering output files are generated, they can be merged using `starcode_combine.py` to create a unified consensus barcode count file.

### 4. Assigning Consensus Barcodes (`add_consensus_bc.py`)
Assigns consensus barcodes to each FASTA record based on Starcode clustering output and updates sequence headers accordingly.

### 5. Sorting Consensus Sequences (`sort_consensus_bc.py`)
Sorts consensus FASTA records alphabetically by header to ensure reproducible ordering for downstream analysis.

### 6. Consensus Gene Generation (`get_consensus_gene.py`)
Sequences are grouped by barcode and aligned using **pyabpoa**. For each barcode group:
- A nucleotide consensus gene sequence is generated
- Alignment scores are recorded

Outputs:
- Consensus gene FASTA file
- Consensus score report

### 7. Consensus Gene Translation (`process_alt_trans.py`)
Translates consensus nucleotide sequences into amino acid sequences and exports results as a CSV file for downstream parsing.

### 8. Map Consensus Genes to Reference Genes (BBMap / Minimap2)
Consensus gene sequences are mapped to the reference gene set (`*.genes`) to generate a SAM file for downstream mutation parsing.

- Default: **BBMap**
- Alternative: **Minimap2**

### 9. SAM File Parsing and Mutation Analysis (`parse_sam_script.py`)
Parses SAM files using:
- Reference protein file (`*.proteins`)
- Barcode assignments

Performs pairwise alignments to:
- Identify mutations
- Assign mutant IDs
- Generate aggregated mutation summaries

Outputs:
- Barcode-to-mutant mapping: `C5seqs_mutID_all.csv`
- Aggregated mutant information: `C5seqs_mutID_info_all.csv`

---

## 🚀 Installation & Setup

**Clone the Repo**

```bash
# SSH Clone
git clone git@github.com:PlesaLab/DropSynth_BC_Mapping.git
cd DropSynth_BC_Mapping

# HTML Clone
git clone https://github.com/PlesaLab/DropSynth_BC_Mapping.git
cd DropSynth_BC_Mapping
```

---

## 📥 Input File Requirements

The pipeline requires the following inputs:

### Raw Sequencing Data
- `*.fastq.gz` — Nanopore sequencing reads

### Reference Design Files
- `*.genes` — reference nucleotide design sequences
- `*.proteins` — reference amino acid sequences

The `.genes` file is used for nucleotide mapping, while `.proteins` is used during mutation classification and amino acid comparisons.

---

## 📤 Output Structure

The pipeline produces intermediate and final outputs organized by analysis stage.

### Barcode-Level Outputs
- `*.bc_stats.csv`
- `*.bc_list.csv`
- `*.bc_consensus.fasta`

### Consensus Gene Outputs
- `*.consensus_gene.fasta`
- `*.consensus_scores.txt`
- `*.translated.csv`

### Mapping Outputs
- `*.sam`

### Mutation Reports
- `C5seqs_mutID_all.csv`
- `C5seqs_mutID_info_all.csv`

These outputs serve as primary inputs for downstream QC analysis, mutation frequency estimation, and library performance evaluation.

---

## 💻 Dependencies

This pipeline relies on Conda-managed dependencies defined in `environment.yml`.

### Core Dependencies
- Python 3.9
- biopython
- pandas
- numpy
- regex

### External Tools
- starcode 1.4
- minimap2 2.28
- samtools 1.21
- bbmap 39.17

### Pip Packages
- edlib
- pyabpoa

### Environment Setup

```bash
conda env create -f environment.yml
conda activate newenv
```

---

## 🚀 Running the Pipeline

The workflow is executed via a `Makefile`, with HPC support provided through `consensus_makefile.sh`.

### consensus_makefile.sh Overview

This script:
1. Sets the project directory
2. Loads required modules
3. Activates the Conda environment
4. Sets thread counts from SLURM
5. Runs the pipeline

```bash
# Set project directory
PROJ_DIR=$(pwd)
cd "$PROJ_DIR"

# Load modules
module purge
module load miniconda3

# Activate environment
conda activate newenv

# Thread configuration
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export THREAD_COUNT=$SLURM_CPUS_PER_TASK

# Run with BBMap (default)
make -j $SLURM_CPUS_PER_TASK

# Alternative: MiniMap
# make -j $SLURM_CPUS_PER_TASK MAPPING_TOOL=minimap
```

### Submitting as SLURM job

```bash
sbatch consensus_makefile.sh
```

### Makefile Targets

The `Makefile` contains targets for each pipeline stage:
- FASTQ splitting
- Barcode processing
- Starcode clustering
- Consensus barcode assignment
- FASTA sorting
- Consensus gene generation
- Translation
- SAM parsing

Users may run individual stages if needed for debugging or iterative analysis.

---

## 🗂️ Configuration and Customization

### Mapping Tool Selection

Default:

```bash
make
```

Minimap:

```bash
make MAPPING_TOOL=minimap
```

### Thread Control

```bash
export THREAD_COUNT=<N>
```

In SLURM:

```bash
$SLURM_CPUS_PER_TASK
```

### Barcode Definitions and Filtering

Parameters for:
- motif detection
- trimming
- filtering thresholds

are defined in:
- `barcode_processing.py`
- `parse_sam_script.py`

These can be modified for:
- alternative barcode architectures
- new library formats
- updated filtering criteria

---

## 🏗️ Performance Notes

- **BBMAP** provides high sensitivity for synthetic constructs with low divergence
- **Minimap2** offers faster mapping for very large datasets (and longer constructs)
- **pyabpoa** improves consensus accuracy for noisy Nanopore reads

Runtime depends on:
- total read depth
- barcode diversity
- thread availability

Recommended resources for large libraries (>10k barcodes):
- ≥64 GB RAM
- ≥16 CPU threads

---

## 🛠️ Troubleshooting

### Conda Environment Issues

```bash
conda env remove -n newenv
conda env create -f environment.yml
conda activate newenv
```

### Starcode Not Found

```bash
which starcode
```

### Empty Consensus Outputs

Possible causes:
- incorrect barcode motif definition
- excessive filtering
- failed barcode clustering

### Low Mapping Rates

Check:
- reference gene file correctness
- mapping tool selection
- read orientation and trimming parameters

---

## 🔁 Reproducibility

This pipeline is designed for reproducible analysis of DropSynth barcoded libraries.

Key reproducibility features:
- deterministic FASTA sorting
- modular scripts with defined inputs/outputs
- Makefile-driven orchestration
- 1environment.yml1 for dependency locking

---

## 📄 License

This repository is released under an academic-use license. See `LICENSE` for details.

---

## ⚙️ Maintainers
 
- Karl Romanowicz (krom@uoregon.edu)
- Calin Plesa (calin@uoregon.edu)
