#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=trim_genome_short_reads 

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --cpus-per-task=40 
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=24:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=64G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out 

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script cleans the raw data using fastp.

#--------------------------------------------------------------------------------

# Load modules  
fastp=/gpfs1/home/e/l/elongman/software/fastp

#--------------------------------------------------------------------------------

# Define important file locations

# RAW_DATA is the core folder where all of the raw data is located.
RAW_DATA=/netfiles/pespenilab_share/Nucella/raw

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/genome_assembly

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "trim_AVITI" ]
then echo "Working trim_AVITI folder exist"; echo "Let's move on."; date
else echo "Working trim_AVITI folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/genome_assembly/trim_AVITI; date
fi

#--------------------------------------------------------------------------------

# Define the reads to be processed
READ1=$RAW_DATA/Project_ESEL_NC3/NC3_R1.fastq.gz
READ2=$RAW_DATA/Project_ESEL_NC3/NC3_R2.fastq.gz

echo "Read1:" ${READ1} "Read2:" ${READ2}

#--------------------------------------------------------------------------------

# Use fastp to trim the AVITI reads
$fastp \
-i ${READ1} -I ${READ2} \
-o $WORKING_FOLDER/data/processed/genome_assembly/trim_AVITI/NC3_R1_clean.fq.gz \
-O $WORKING_FOLDER/data/processed/genome_assembly/trim_AVITI/NC3_R2_clean.fq.gz \
--detect_adapter_for_pe \
--trim_front1 8 \
--trim_poly_g \
--thread 16 \
--cut_right \
--cut_window_size 6 \
--qualified_quality_phred 20 \
--length_required 35 \
--html $WORKING_FOLDER/data/processed/genome_assembly/trim_AVITI/NC3_fastp.html \
--json $WORKING_FOLDER/data/processed/genome_assembly/trim_AVITI/fastp/NC3_fastp.json
