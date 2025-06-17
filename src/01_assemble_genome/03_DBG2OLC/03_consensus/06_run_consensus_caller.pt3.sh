#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=consensus_pt3

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss 
#SBATCH --time=28:00:00 

# Request memory for the entire job
#SBATCH --mem=60G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out 

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run the third step in the consensus script.

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

#--------------------------------------------------------------------------------

final_dir=$WORKING_FOLDER/data/processed/genome_assembly/consensus/final_assembly_array

# Cat files together to produce a final consensus
for file in $(find ${final_dir} -name "*.consensus.fasta"); do
cmd="cat ${file};"
eval $cmd
done > $WORKING_FOLDER/data/processed/genome_assembly/consensus/final_assembly.fasta