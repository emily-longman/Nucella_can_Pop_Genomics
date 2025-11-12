#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=moments_admix

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=5:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=20G 

# Submit job array
#SBATCH --array=0-2

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Set up environemnt
#load dependencies
#module load python3.12-anaconda/2024.06-1

#Step 1. Create conda environemt
#This needs to be done  only once

#conda create \
#-n moments_ekl \
#python=3.12 \
#dadi \
#ipykernel \
#-c conda-forge

# Activates conda environemnt
#source activate moments_ekl

#Installs moments
#conda install moments -c bioconda

# install pandas
#conda install pandas

#--------------------------------------------------------------------------------

# Load modules 
module load python3.12-anaconda/2024.06-1
source ${ANACONDA_ROOT}/etc/profile.d/conda.sh
#conda activate moments_jcbn
conda activate moments_ekl

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/pop_structure

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "admix_relaxed_filt" ]
then echo "Working admix_relaxed_filt folder exist"; echo "Let's move on."; date
else echo "Working admix_relaxed_filt folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/pop_structure/admix_relaxed_filt; date
fi

#--------------------------------------------------------------------------------

# Change directory
cd $WORKING_FOLDER/data/processed/pop_structure

# Move guide file
scp $WORKING_FOLDER/guide_files/trios_guide.txt .

moment_script="$WORKING_FOLDER/src/03_population_structure/06_simulations_and_inference/01_moments/02_nucella.run_moments_admix_relaxed_filt.py"

  python $moment_script \
	${SLURM_ARRAY_TASK_ID}

#--------------------------------------------------------------------------------

conda deactivate

#--------------------------------------------------------------------------------

### Print the time
  echo "ended at"  `date`
