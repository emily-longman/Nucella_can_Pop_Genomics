#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=breaker_sing

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=8

# Reserve walltime -- hh:mm:ss --7 day limit 
#SBATCH --time=6-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=30G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Load modules
module load singularity/3.7.1

#Define important file locations

#Working folder is core folder where this pipeline is being run.
SCRIPTS_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/00_Genome_short_read

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $SCRIPTS_FOLDER

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "23_braker_singularity" ]
then echo "Working 23_braker_singularity folder exist"; echo "Let's move on."; date
else echo "Working 23_braker_singularity folder doesnt exist. Let's fix that."; mkdir $SCRIPTS_FOLDER/23_braker_singularity; date
fi

#--------------------------------------------------------------------------------

# Change directory
cd $SCRIPTS_FOLDER/23_braker_singularity

# Build the sif-file for the braker singularity (https://hub.docker.com/r/teambraker/braker3)

#singularity build braker3.sif docker://teambraker/braker3:latest

# Test with: 
# singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test1.sh . # Get test script
# export BRAKER_SIF=$SCRIPTS_FOLDER/23_breaker_singularity/braker3.sif
# bash test1.sh

#--------------------------------------------------------------------------------

# Change directory
cd $SCRIPTS_FOLDER

# Export path to braker sif
export BRAKER_SIF=$SCRIPTS_FOLDER/23_braker_singularity/braker3.sif

# Change permissions for associated script
chmod +x 23_G_braker.sh

# Use the singularity to run angsd
singularity exec $SCRIPTS_FOLDER/23_braker_singularity/braker3.sif $SCRIPTS_FOLDER/23_G_braker.sh