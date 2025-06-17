#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=fltlong 

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --cpus-per-task=40 
#SBATCH --nodes=1

# Reserve walltime -- hh:mm:ss
#SBATCH --time=24:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=80G 

# Submit job array
#SBATCH --array=0-5

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/fltlong.%A_%a.out 

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Filter ONT data by length using Filtlong (https://github.com/rrwick/Filtlong)

#--------------------------------------------------------------------------------

# Load modules  
filtlong=/gpfs1/home/e/l/elongman/software/Filtlong/bin/filtlong

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "ONT_fltlong" ]
then echo "Working ONT_fltlong folder exist"; echo "Let's move on."; date
else echo "Working ONT_fltlong folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/ONT_fltlong; date
fi

#--------------------------------------------------------------------------------

# Specify array filter lengths
arr=(1000 2000 3500 5000 7500 10000)
L="${arr[$SLURM_ARRAY_TASK_ID]}"
echo ${L}

#--------------------------------------------------------------------------------

# Filter ONT using Filtlong
$filtlong \
--min_length $L \
$WORKING_FOLDER/data/raw/ONT/FC_all.ONT.nuc.fastq.gz | gzip > $WORKING_FOLDER/data/processed/ONT_fltlong/Nuc.$L.fltlong.fastq.gz