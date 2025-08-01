#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=index_genome

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=24:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=60G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script with index the reference genome.

#--------------------------------------------------------------------------------

# Load modules 
module load gcc/13.3.0-xp3epyt
module load samtools/1.19.2-pfmpoam
bwa=/netfiles/nunezlab/Shared_Resources/Software/bwa-mem2-2.2.1_x64-linux/bwa-mem2.avx2
PICARD=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

# Genome from first round of pilon
REFERENCE=$WORKING_FOLDER/data/processed/genome_assembly/pilon/polished_genome_round_1/polished_assembly.fasta

#--------------------------------------------------------------------------------

# Move to the directory that the genome is currently stored
cd $WORKING_FOLDER/data/processed/genome_assembly/pilon/polished_genome_round_1

# Index database sequences in the FASTA format 
$bwa index $REFERENCE 

# Generate the FASTA sequence dictionary file
java -jar $PICARD CreateSequenceDictionary \
R=$REFERENCE

# Generate the fasta index file for the reference
samtools faidx $REFERENCE