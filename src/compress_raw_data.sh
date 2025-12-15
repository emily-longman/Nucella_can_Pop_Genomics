#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=compress_raw_data

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss
#SBATCH --time=30:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=100G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will move the raw data to gpfs3 and compress it. 

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

# NETFILES is for long term storage
NETFILES=/netfiles/pespenilab_share/Nucella

# TMP_FOLDER is the path to gpfs3 where there is lots of storage
TMP_FOLDER=/gpfs3/scratch/elongman

#--------------------------------------------------------------------------------

# Copy Population genomics raw data from Netfiles to TMP_FOLDER

# Change directory
cd $TMP_FOLDER

# Copy folder
scp -r $NETFILES/raw/Population_genomics $TMP_FOLDER

#--------------------------------------------------------------------------------

# Tar and gzip file
#tar -cvzf Population_genomics.tar.gz Population_genomics

echo "done"