#!/bin/sh
#
#SBATCH -J genoTrios
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=5:00:00
#SBATCH -o ./slurmOutput/genom.%A_%a.out # Standard output
#SBATCH -p general
#SBATCH --array=1-3

#--------------------------------------------------------------------------------

# Load modules 
module load Rgeospatial

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

#--------------------------------------------------------------------------------

cd $WORKING_FOLDER

guide=$WORKING_FOLDER/guide_files/trios_guide.txt

#--------------------------------------------------------------------------------

# Run accompanying R script
Rscript --vanilla $WORKING_FOLDER/src/03_population_structure/06_simulations_and_inference/01_moments/01_Genomalicious.Trios.Nucella.R \
${SLURM_ARRAY_TASK_ID} \
$guide
