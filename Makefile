SHELL := /bin/bash

# Define Directories
PROJ_DIR := $(shell pwd)
INPUT_DIR := $(PROJ_DIR)/libs
REF_DIR := $(PROJ_DIR)/refs
SCRIPT_DIR := $(PROJ_DIR)/scripts
OUTPUT_DIR := $(PROJ_DIR)/out

# Define Variables
INPUT_FASTQ := $(wildcard $(INPUT_DIR)/*.merged.fastq.gz)
OUTPUT_PREFIX ?= $(OUTPUT_DIR)/part_
REF_GENOME := $(wildcard $(REF_DIR)/*.genes)
PROTEIN_FASTA := $(wildcard $(REF_DIR)/*.proteins)
BC_INFO_FILE := merged.sorted_noN_aa
TRANS_FILES := merged.sorted_noN_aa.csv
SAM_FILES := merged.sorted

# Define Software
PYTHON := python
BBMAP := bbmap.sh
SAMTOOLS := samtools
MINIMAP := minimap2
STARCODE := starcode

# Define Python scripts
SPLIT_SCRIPT := $(SCRIPT_DIR)/split_script.py
BC_PROCESSING := $(SCRIPT_DIR)/barcode_processing.py
STARCODE_COMBINE := $(SCRIPT_DIR)/starcode_combine.py
ADD_CONSENSUS_SCRIPT := $(SCRIPT_DIR)/add_consensus_bc.py
SORT_CONSENSUS_SCRIPT := $(SCRIPT_DIR)/sort_consensus_bc.py
CONSENSUS_GENE_SCRIPT := $(SCRIPT_DIR)/get_consensus_gene.py
ALT_TRANS_SCRIPT := $(SCRIPT_DIR)/process_alt_trans.py
SAM_PARSE_SCRIPT := $(SCRIPT_DIR)/parse_sam_script.py

# Define number of threads
THREADS ?= $(shell echo $$THREAD_COUNT)
THREADS := $(or $(THREADS),20)

# Choose mapping tool: bbmap or minimap
MAPPING_TOOL ?= bbmap

#========================================================================
# Primary targets
all: prepare split process_barcodes combine_fasta combine_starcode run_starcode add_consensus_bc sort_consensus_bc consensus_gene process_alt_trans run_mapping parse_sam

# Declare phony targets
.PHONY: all prepare split process_barcodes combine_fasta combine_starcode run_starcode add_consensus_bc sort_consensus_bc consensus_gene process_alt_trans run_mapping parse_sam

# Use secondary expansion for pattern rules
.SECONDEXPANSION:

# Generate lists of output files for all steps (number of files may vary based on splitting)
SPLIT_FILES := $(wildcard $(OUTPUT_PREFIX)*.gz)
PROCESSED_BC_FILES := $(patsubst %.gz, %.sorted.fasta, $(SPLIT_FILES))
STARCODE_FILES := $(patsubst %.gz, %.bc_stats_for_starcode.tsv, $(SPLIT_FILES))
MERGED_FASTA := $(OUTPUT_DIR)/merged.sorted_noN.fasta
STARCODE_INPUT := $(OUTPUT_DIR)/merged.bc_stats_for_starcode.tsv
STARCODE_OUTPUT := $(OUTPUT_DIR)/barcodes.collapse_d1.tsv
ADD_CONSENSUS_OUTPUT := $(OUTPUT_DIR)/merged.consensus_bc.fasta
SORTED_CONSENSUS_OUTPUT := $(OUTPUT_DIR)/merged.consensus_bc.sorted.fasta
CONSENSUS_GENE_FASTA := $(OUTPUT_DIR)/consensus_gene.fasta
CONSENSUS_SCORES := $(OUTPUT_DIR)/consensus_scores.txt
ALT_TRANS_OUTPUT := $(OUTPUT_DIR)/merged.sorted_noN_aa.csv
BBMAP_OUTPUT := $(OUTPUT_DIR)/merged.sorted.map.sam
BBMAP_UNALIGNED := $(OUTPUT_DIR)/merged.sorted.map.unaligned.sam
MINIMAP_OUTPUT := $(OUTPUT_DIR)/merged.sorted.map.sam
MINIMAP_UNALIGNED := $(OUTPUT_DIR)/merged.sorted.map.unaligned.fastq
PARSE_SAM_OUTPUT1 := $(OUTPUT_DIR)/C5seqs_mutID_all.csv
PARSE_SAM_OUTPUT2 := $(OUTPUT_DIR)/C5seqs_mutID_info_all.csv

# Ensure intermediate files are not deleted
.SECONDARY: $(SPLIT_FILES) $(PROCESSED_BC_FILES) $(STARCODE_FILES) $(MERGED_FASTA) $(STARCODE_INPUT) $(STARCODE_OUTPUT) $(ADD_CONSENSUS_OUTPUT) $(SORTED_CONSENSUS_OUTPUT) $(CONSENSUS_GENE_FASTA) $(CONSENSUS_SCORES) $(ALT_TRANS_OUTPUT) $(BBMAP_OUTPUT) $(BBMAP_UNALIGNED) $(MINIMAP_OUTPUT) $(MINIMAP_UNALIGNED) $(PARSE_SAM_OUTPUT1) $(PARSE_SAM_OUTPUT2)

# Retain only for interactive debugging to run steps independently
.PHONY: split_only process_barcodes_only combine_fasta_only combine_starcode_only run_starcode_only add_consensus_bc_only sort_consensus_bc_only consensus_gene_only process_alt_trans_only run_bbmap_only run_minimap_only parse_sam_only_bbmap parse_sam_only_minimap

split_only:
	echo "Splitting FASTQ file: $(INPUT_FASTQ)"
	$(PYTHON) $(SPLIT_SCRIPT) $(INPUT_FASTQ) $(OUTPUT_PREFIX)
	echo "Splitting complete."

process_barcodes_only:
	echo "Processing barcodes from split FASTQ files"
	for file in $(OUTPUT_PREFIX)*.gz; do \
		echo "Processing $${file}"; \
		$(PYTHON) $(BC_PROCESSING) "$${file}" "$${file%.gz}" --max_threads $(THREADS) --start_motif CATATG --end_motif TAAGGTACCTAA || exit 1; \
	done
	echo "Barcode processing complete."

combine_fasta_only:
	echo "Combining FASTA files"
	cat $(OUTPUT_PREFIX)*.sorted.fasta > $(MERGED_FASTA)
	echo "FASTA files combined."

combine_starcode_only:
	echo "Combining and processing Starcode files"
	$(PYTHON) $(STARCODE_COMBINE) "$(OUTPUT_PREFIX)*.bc_stats_for_starcode.tsv" $(STARCODE_INPUT)
	echo "Starcode files combined."

run_starcode_only:
	echo "Running Starcode for clustering barcodes"
	$(STARCODE) -d1 --sphere --print-clusters -i $(STARCODE_INPUT) --output $(STARCODE_OUTPUT)
	echo "Starcode clustering complete."

add_consensus_bc_only:
	echo "Adding consensus barcode using add_consensus_bc.py"
	$(PYTHON) $(ADD_CONSENSUS_SCRIPT) -i $(MERGED_FASTA) -c $(STARCODE_OUTPUT) -o $(ADD_CONSENSUS_OUTPUT)
	echo "Consensus barcode addition complete."

sort_consensus_bc_only:
	echo "Sorting consensus barcode file..."
	$(PYTHON) $(SORT_CONSENSUS_SCRIPT) -i $(ADD_CONSENSUS_OUTPUT) -o $(SORTED_CONSENSUS_OUTPUT)
	echo "Sorting complete."

consensus_gene_only:
	echo "Calling consensus gene using get_consensus_gene.py"
	$(PYTHON) $(CONSENSUS_GENE_SCRIPT) -i $(SORTED_CONSENSUS_OUTPUT) -o $(CONSENSUS_GENE_FASTA) -s $(CONSENSUS_SCORES) -p $(THREADS)
	echo "Consensus gene calling complete."

process_alt_trans_only:
	echo "Processing alternative transcripts"
	$(PYTHON) $(ALT_TRANS_SCRIPT) $(CONSENSUS_GENE_FASTA) $(ALT_TRANS_OUTPUT)
	echo "Alternative transcript processing completed."

run_bbmap_only:
	echo "Running BBMap"
	$(BBMAP) ref=$(REF_GENOME) in=$(CONSENSUS_GENE_FASTA) outm=$(BBMAP_OUTPUT) outu=$(BBMAP_UNALIGNED)
	echo "BBMap completed."

run_minimap_only:
	$(MINIMAP) -ax map-ont -o $(MINIMAP_OUTPUT) $(REF_GENOME) $(CONSENSUS_GENE_FASTA)
	echo "Extracting unmapped reads"
	$(SAMTOOLS) fastq -f 4 $(MINIMAP_OUTPUT) > $(MINIMAP_UNALIGNED)
	echo "minimap2 completed."

parse_sam_only_bbmap:
	echo "Parsing SAM alignments"
	$(PYTHON) $(SAM_PARSE_SCRIPT) "$(PROTEIN_FASTA)" "$(ALT_TRANS_OUTPUT)" "$(BBMAP_OUTPUT)" "$(OUTPUT_DIR)/"
	echo "SAM parsing completed."

parse_sam_only_minimap:
	echo "Parsing SAM alignments"
	$(PYTHON) $(SAM_PARSE_SCRIPT) "$(PROTEIN_FASTA)" "$(ALT_TRANS_OUTPUT)" "$(MINIMAP_OUTPUT)" "$(OUTPUT_DIR)/"
	echo "SAM parsing completed."

#========================================================================
# STEP 1: Create Output Directories (as needed)
prepare:
	mkdir -p $(OUTPUT_DIR)

#========================================================================
# STEP 2: Target to split the fastq.gz file
split: prepare
	echo "Splitting FASTQ file: $(INPUT_FASTQ)"
	$(PYTHON) $(SPLIT_SCRIPT) $(INPUT_FASTQ) $(OUTPUT_PREFIX)
	echo "Splitting complete."

#========================================================================
# STEP 3: Process barcodes from split FASTQ files
process_barcodes: split
	echo "Processing barcodes from split FASTQ files"
	for file in $(OUTPUT_PREFIX)*.gz; do \
	  echo "Processing $$file"; \
	  $(PYTHON) $(BC_PROCESSING) "$$file" "$${file%.gz}" --max_threads $(THREADS) --start_motif CATATG --end_motif TAAGGTACCTAA || exit 1; \
	done
	echo "Barcode processing complete."

#========================================================================
# STEP 4: Combining FASTA files
combine_fasta: process_barcodes
	echo "Combining FASTA files"
	cat $(OUTPUT_PREFIX)*.sorted.fasta > $(MERGED_FASTA)
	echo "FASTA files combined."

#========================================================================
# STEP 5: Combining and processing Starcode files
combine_starcode: process_barcodes
	echo "Combining and processing Starcode files"
	$(PYTHON) $(STARCODE_COMBINE) "$(OUTPUT_PREFIX)*.bc_stats_for_starcode.tsv" $(STARCODE_INPUT)
	echo "Starcode files combined."

#========================================================================
# STEP 6: Running Starcode for clustering barcodes
run_starcode: combine_starcode
	echo "Running Starcode for clustering barcodes"
	$(STARCODE) -d1 --sphere --print-clusters -i $(STARCODE_INPUT) --output $(STARCODE_OUTPUT)
	echo "Starcode clustering complete."

#========================================================================
# STEP 7 (new): Add consensus barcodes using add_consensus_bc.py
add_consensus_bc: run_starcode
	echo "Adding consensus barcode using add_consensus_bc.py"
	$(PYTHON) $(ADD_CONSENSUS_SCRIPT) -i $(MERGED_FASTA) -c $(STARCODE_OUTPUT) -o $(ADD_CONSENSUS_OUTPUT)
	echo "Consensus barcode addition complete."

#========================================================================
# STEP 8 (new): Sort the consensus barcode file
sort_consensus_bc: add_consensus_bc
	echo "Sorting consensus barcode file..."
	$(PYTHON) $(SORT_CONSENSUS_SCRIPT) -i $(ADD_CONSENSUS_OUTPUT) -o $(SORTED_CONSENSUS_OUTPUT)
	echo "Sorting complete."

#========================================================================
# STEP 9 (new): Generate consensus gene sequences using get_consensus_gene.py
consensus_gene: sort_consensus_bc
	echo "Calling consensus gene using get_consensus_gene.py"
	$(PYTHON) $(CONSENSUS_GENE_SCRIPT) -i $(SORTED_CONSENSUS_OUTPUT) -o $(CONSENSUS_GENE_FASTA) -s $(CONSENSUS_SCORES) -p $(THREADS)
	echo "Consensus gene calling complete."

#========================================================================
# STEP 10: Processing Alternative Transcripts
process_alt_trans: consensus_gene
	echo "Processing alternative transcripts"
	$(PYTHON) $(ALT_TRANS_SCRIPT) $(CONSENSUS_GENE_FASTA) $(ALT_TRANS_OUTPUT)
	echo "Alternative transcript processing completed."

#========================================================================
# STEP 11: Running BBMap or MiniMap (modified to use consensus gene FASTA)
run_mapping: process_alt_trans
ifeq ($(MAPPING_TOOL),bbmap)
	echo "Running BBMap"
	$(BBMAP) ref=$(REF_GENOME) in=$(CONSENSUS_GENE_FASTA) outm=$(BBMAP_OUTPUT) outu=$(BBMAP_UNALIGNED)
	echo "BBMap completed."
else ifeq ($(MAPPING_TOOL),minimap)
	echo "Running MiniMap"
	$(MINIMAP) -ax map-ont -o $(MINIMAP_OUTPUT) $(REF_GENOME) $(CONSENSUS_GENE_FASTA)
	echo "Extracting unmapped reads"
	$(SAMTOOLS) fastq -f 4 $(MINIMAP_OUTPUT) > $(MINIMAP_UNALIGNED)
	echo "minimap2 completed."
else
	$(error Invalid MAPPING_TOOL specified. Use "bbmap" or "minimap")
endif

#========================================================================
# STEP 12: Parsing SAM files with new consensus gene FASTA
parse_sam: run_mapping
ifeq ($(MAPPING_TOOL),bbmap)
	echo "Parsing SAM alignments from BBMap"
	$(PYTHON) $(SAM_PARSE_SCRIPT) "$(PROTEIN_FASTA)" "$(ALT_TRANS_OUTPUT)" "$(BBMAP_OUTPUT)" "$(OUTPUT_DIR)/"
	echo "SAM parsing with BBMap completed."
else ifeq ($(MAPPING_TOOL),minimap)
	echo "Parsing SAM alignments from MiniMap"
	$(PYTHON) $(SAM_PARSE_SCRIPT) "$(PROTEIN_FASTA)" "$(ALT_TRANS_OUTPUT)" "$(MINIMAP_OUTPUT)" "$(OUTPUT_DIR)/"
	echo "SAM parsing with MiniMap completed."
else
	$(error Invalid MAPPING_TOOL specified. Use "bbmap" or "minimap")
endif
