#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=download_short_reads

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss
#SBATCH --time=07-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=100G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out 

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# Download AVITI short reads from UC Davis DNA Technologies Core.

#--------------------------------------------------------------------------------

# Define important file locations

# RAW_DATA is the core folder where all of the raw data is located.
RAW_DATA=/netfiles/pespenilab_share/Nucella/raw

#--------------------------------------------------------------------------------

# Move to the directory where you want the output files to be saved
cd $RAW_DATA

# Use rsync to download the data
rsync -avL slimsdata.genomecenter.ucdavis.edu::slims/v9czo36oz7/Unaligned/Project_ESEL_NC3 .