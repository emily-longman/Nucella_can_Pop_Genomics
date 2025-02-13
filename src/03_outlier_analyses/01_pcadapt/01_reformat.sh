#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=reformat_VCF

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=10:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=8G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will use plink to reformat the vcf so that it is in bed format.

# Load software
plink=/gpfs1/home/e/l/elongman/software/plink_linux_x86_64_20241022/plink

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/outlier_analyses

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "pcadapt" ]
then echo "Working pcadapt folder exist"; echo "Let's move on."; date
else echo "Working pcadapt folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/outlier_analyses/pcadapt; date
fi

#--------------------------------------------------------------------------------