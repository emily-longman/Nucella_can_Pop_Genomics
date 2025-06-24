#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=mitogenome_analysis 

# Specify partition
#SBATCH --partition=general

# Request nodes

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=30:00:00

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=20G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu

#--------------------------------------------------------------------------------

# Determine mitochondrial data for NCBI
# Note: Using a N. lapillus mitogenome (https://www.ncbi.nlm.nih.gov/search/all/?term=PP147811) 
# assembled as part of Garcia-Souto et al. 2024 paper (https://link.springer.com/article/10.1007/s00227-024-04424-3)

#--------------------------------------------------------------------------------

# Load modules 
# Call package (installed with conda)
module load python3.11-anaconda/2024.02-1
source ${ANACONDA_ROOT}/etc/profile.d/conda.sh
#conda create --name exonerate #If you haven't already done so, create and name the environment
source activate exonerate #activate the environment
#conda install -c bioconda exonerate # If you haven't already done so, install the program

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

# Mitogenome - downloaded from https://www.ncbi.nlm.nih.gov/search/all/?term=PP147811
mitogenome=$WORKING_FOLDER/data/processed/mitochondrial_genome/PP147811.1.fna

# N. canaliculata draft genome
genome=/gpfs3/scratch/elongman/NCBI_genome/N.canaliculata_assembly.fasta

#--------------------------------------------------------------------------------

# Change directory
cd $WORKING_FOLDER/data/processed/mitochondrial_genome_2

# Align N. lapillus mitogenome to draft N. canaliculata genome
exonerate --model est2genome $mitogenome $genome

#--------------------------------------------------------------------------------

# Deactivate conda
conda deactivate