#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=moments_3st_alt_mutation

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=30:00:00 

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

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

#--------------------------------------------------------------------------------

# Load modules 
module load python3.12-anaconda/2024.06-1
source ${ANACONDA_ROOT}/etc/profile.d/conda.sh
#conda activate moments_jcbn
#conda activate moments_ekl

# Use JCBN conda environment
# Move to working directory
#cd $WORKING_FOLDER/guide_files
#conda env create -f moments_jcbn.environment.yml
conda activate moments_jcbn

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/pop_structure

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "o3step_alt_mutation" ]
then echo "Working o3step_alt_mutation folder exist"; echo "Let's move on."; date
else echo "Working o3step_alt_mutation folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/pop_structure/o3step_alt_mutation; date
fi

#--------------------------------------------------------------------------------

# Change directory
cd $WORKING_FOLDER/data/processed/pop_structure

# Move guide file
scp $WORKING_FOLDER/guide_files/trios_guide.txt .

moment_script="$WORKING_FOLDER/src/03_population_structure/06_simulations_and_inference/01_moments/02_nucella.run_moments_3stepstone_alt_mutation.py"

  python $moment_script \
	${SLURM_ARRAY_TASK_ID}

#--------------------------------------------------------------------------------

conda deactivate

#--------------------------------------------------------------------------------

### Print the time
  echo "ended at"  `date`
