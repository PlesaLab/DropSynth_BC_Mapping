#!/bin/bash
#SBATCH --account=plesalab
#SBATCH --partition=compute
#SBATCH --job-name=lib1_BCmapping_consensus
#SBATCH --output=/home/%u/%x-%j.log
#SBATCH --error=/home/%u/%x-%j.log
#SBATCH --time=24:00:00
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=krom@uoregon.edu
#SBATCH --mail-type=BEGIN,END,FAIL

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