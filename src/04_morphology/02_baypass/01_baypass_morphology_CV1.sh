#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=baypass_morphology_CV_1

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=3-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=40G 

# Request CPU
#SBATCH --cpus-per-task=20

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run baypass using PC1 and PC2 of the morphology shape data and using an omega file. 

#--------------------------------------------------------------------------------

# Load modules 
module load gcc/13.3.0-xp3epyt
baypass=/gpfs1/home/e/l/elongman/software/baypass_public/sources/g_baypass

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/morphometrics

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "baypass" ]
then echo "Working baypass folder exist"; echo "Let's move on."; date
else echo "Working baypass folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/morphometrics/baypass; date
fi

# Change directory 
cd $WORKING_FOLDER/data/processed/morphometrics/baypass

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "CV_1" ]
then echo "Working CV_1 folder exist"; echo "Let's move on."; date
else echo "Working CV_1 folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/morphometrics/baypass/CV_1; date
fi

#--------------------------------------------------------------------------------

# Change directory 
cd $WORKING_FOLDER/data/processed/morphometrics/baypass/CV_1

# Run baypass with morphometric data
$baypass -npop 19 \
-gfile $WORKING_FOLDER/data/processed/outlier_analyses/baypass/genobaypass \
-poolsizefile $WORKING_FOLDER/data/processed/outlier_analyses/baypass/poolsize \
-omegafile $WORKING_FOLDER/data/processed/outlier_analyses/baypass/omega/NC_baypass_mat_omega.out \
-efile $WORKING_FOLDER/data/processed/morphometrics/guide_files/Baypass_morphology_CV_1.txt \
-d0yij 4 \
-outprefix NC_baypass_morphology_CV_1 \
-nthreads 20

