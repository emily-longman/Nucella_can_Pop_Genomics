#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=baypass_omega

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=30:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=40G 

# Request CPU
#SBATCH --cpus-per-task=5

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will perform the first step in running baypass. 

# Load modules 
module load gcc/10.5.0
baypass=/gpfs1/home/e/l/elongman/software/baypass_public/sources/g_baypass

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed

#--------------------------------------------------------------------------------

# Change directory 
cd $WORKING_FOLDER/GEA/baypass/omega_file

# Run baypass - this will generate the omega file which will be used in the subsequent scripts
$baypass -npop 19 \
-gfile $WORKING_FOLDER/GEA/baypass/genobaypass \
-poolsizefile $WORKING_FOLDER/GEA/baypass/poolsize \
-d0yij 8 \
-outprefix NC_baypass \
-npilot 100 -nthreads 5

#-npop: number of pools
#-gfile: gfile input
#-poolsizefile: poolsize file
#-d0yij = is something to do with the initial read counts and their distribution seeded into the model, and is specifically for Pool-seq data. I followed the tip in the manual that says to set it to 1/5th of the minimum pool size
#-outprefix = the header of the output files
#-npilot = number of pilot runs