#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=1:00:00
#SBATCH -o ./SlurmOut/slim.%A_%a.out # Standard output
#SBATCH -p general
#SBATCH --array=1-500

#--------------------------------------------------------------------------------

# Load modules 
module load slim

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

#--------------------------------------------------------------------------------

repid=${SLURM_ARRAY_TASK_ID}
root=$WORKING_FOLDER/data/processed/pop_structure/slim_glacial_refugiaNorth

cd $WORKING_FOLDER/data/processed/pop_structure
mkdir slim_glacial_refugiaNorth

#--------------------------------------------------------------------------------
    
    slim \
	-d "repId=${repid}" \
	-d "root='${root}'" \
     $WORKING_FOLDER/src/03_population_structure/06_simulations_and_inference/02_SLIM_sims/01_Glaciation_refugiaNorth.slim
     
