#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=fastqc_Nucella_ONT_filt_2000

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --cpus-per-task=1 
#SBATCH --nodes=1 

#SBATCH -c 5 #number of cores

# Reserve walltime --time=<dd-hh:mm:ss>
#SBATCH --time=30:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=1000G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# Load modules 
module load gcc/13.3.0-xp3epyt
module load fastqc/0.12.1-qxseug5

#--------------------------------------------------------------------------------

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed

#--------------------------------------------------------------------------------

# Generate Folders and Files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

# Move to working directory
cd $WORKING_FOLDER

# Generating new folder
if [ -d "qc_ONT" ]
then echo "Working qc_ONT folder exist"; echo "Let's move on"; date
else echo "Working qc_ONT folder doesnt exist. Let's fix that"; mkdir $WORKING_FOLDER/qc_ONT; date
fi

#--------------------------------------------------------------------------------

# Lets do some QC on the reads
fastqc $WORKING_FOLDER/Nuc.2000.fltlong.fastq \
--outdir $WORKING_FOLDER/qc_ONT

#--------------------------------------------------------------------------------

echo "done"