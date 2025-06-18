#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=guide_file_pop_names

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss
#SBATCH --time=10:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=1G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will generate a guide file with the population names in the same order as the vcf.

#--------------------------------------------------------------------------------

# Load modules 
module load gcc/13.3.0-xp3epyt
module load bcftools/1.19-iq5mwek

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

#--------------------------------------------------------------------------------

# Make file with population names
bcftools query -l $WORKING_FOLDER/data/processed/fastq_to_vcf/vcf_freebayes/N.canaliculata_pops.vcf.gz > $WORKING_FOLDER/guide_files/N.canaliculata_pops.vcf_pop_names.txt
