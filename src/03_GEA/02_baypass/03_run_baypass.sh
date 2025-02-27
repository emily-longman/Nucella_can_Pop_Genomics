#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=baypass_run

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Request one task
#SBATCH --ntasks=8

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=30:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=60G 

# Request CPU
#SBATCH --cpus-per-task=1

# Submit job array
#SBATCH --array=1

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will perform the first step in running baypass. 

# Load modules 
module load gcc/13.3.0-xp3epyt
module load openmpi/5.0.5
baypass=/gpfs1/home/e/l/elongman/software/baypass_public/sources/g_baypass

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed

#--------------------------------------------------------------------------------

# Change directory 
cd $WORKING_FOLDER/GEA/baypass/baypass_chunks




# Split genotype file into 1000 chunks (you would need a script or command for this)
# Assume `chunk_$SLURM_ARRAY_TASK_ID.txt` are the split chunks

# Run BayPass for each chunk
/users/d/s/dsadler1/programmes/baypass_public-master/sources/./g_baypass \
-gfile genotypefilechunk_$SLURM_ARRAY_TASK_ID.txt \
-outprefix noomega_${outputfile}_$SLURM_ARRAY_TASK_ID -nthreads 8





# Run baypass - this will generate the omega file which will be used in the subsequent scripts
$baypass -npop 19 \
-gfile $WORKING_FOLDER/GEA/baypass/genobaypass \
-poolsizefile $WORKING_FOLDER/GEA/baypass/poolsize \
-d0yij 4 \
-outprefix NC_baypass \
-nthreads 8



# Split genotype file into 1000 chunks (you would need a script or command for this)
# Assume `chunk_$SLURM_ARRAY_TASK_ID.txt` are the split chunks

# Run BayPass for each chunk
/users/d/s/dsadler1/programmes/baypass_public-master/sources/./g_baypass -gfile genotypefilechunk_$SLURM_ARRAY_TASK_ID.txt -efile ${inputfile} -outprefix noomega_${outputfile}_$SLURM_ARRAY_TASK_ID -nthreads 8