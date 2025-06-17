#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=convert_ONT_FQtoFA

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1  

# Request CPUs per task
#SBATCH -c 1

# Reserve walltime -- hh:mm:ss --30 hour limit
#SBATCH --time=10:00:00 

# Request memory for the entire job 
#SBATCH --mem=60G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out 

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will convert the filtered ONT data from fastq to fasta format. 
# NOTE: make sure you use the same filter length as what was used for SparseAssembler and DBG2OLC.

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/genome_assembly

# This part of the script will check and generate, if necessary, all of the output folders used in the script

# Make consensus directory
if [ -d "consensus" ]
then echo "Working consensus folder exist"; echo "Let's move on."; date
else echo "Working consensus folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/genome_assembly/consensus; date
fi

#--------------------------------------------------------------------------------

# Convert ONT from fastq to fasta
cat $WORKING_FOLDER/data/processed/genome_assembly/ONT_fltlong/Nuc.2000.fltlong.fastq | sed -n '1~4s/^@/>/p;2~4p' > $WORKING_FOLDER/data/processed/genome_assembly/consensus/Nuc.2000.fltlong.FQtoFA.fasta