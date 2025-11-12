#!/bin/sh
#
#SBATCH -J moments_admix
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=5:00:00
#SBATCH -o ./slurmOutput/admix.%A_%a.out # Standard output
#SBATCH -p general
#SBATCH --array=0-2

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
WORKING_FOLDER=.

#--------------------------------------------------------------------------------

# Change directory
cd $WORKING_FOLDER/data/processed/pop_structure

# Move guide file
scp $WORKING_FOLDER/guide_files/trios_guide.txt .

moment_script="$WORKING_FOLDER/src/03_population_structure/06_simulations_and_inference/01_moments/02_nucella.run_moments_admix.py"

  python $moment_script \
	${SLURM_ARRAY_TASK_ID}

#--------------------------------------------------------------------------------

conda deactivate

#--------------------------------------------------------------------------------

### Print the time
  echo "ended at"  `date`
