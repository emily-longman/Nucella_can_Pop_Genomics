#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=download_short_reads

# Specify partition
#SBATCH --partition=bigmemwk

# Request nodes
#SBATCH --ntasks-per-node=1 # this is number of CPUs you want to use for parallel computing [also referred to as threads] 

# Reserve walltime -- hh:mm:ss --7 day limit
#SBATCH --time=07-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=100G #<= this may depend on your resources

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# Download short reads

# Move to the directory where you want the output files to be saved (this is one folder higher than where first downloaded the data)
cd /netfiles/pespenilab_share/Nucella/raw/

# Use rsync to download the data
rsync -avL slimsdata.genomecenter.ucdavis.edu::slims/v9czo36oz7/Unaligned/Project_ESEL_NC3 .