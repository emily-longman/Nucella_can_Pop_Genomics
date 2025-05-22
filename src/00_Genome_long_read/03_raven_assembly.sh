#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=raven 

# Specify partition
#SBATCH --partition=bigmemwk

# Request nodes
#SBATCH --cpus-per-task=1 
#SBATCH --nodes=1 # on one node

# Reserve walltime --time=<dd-hh:mm:ss>
#SBATCH --time=07-00:00:00 #<= this may depend on your resources

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=900G #<= this may depend on your resources

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------
# Call package (installed with conda)
module load python3.11-anaconda/2023.09-0
source ${ANACONDA_ROOT}/etc/profile.d/conda.sh
conda create --name raven #create and name the environment
conda activate raven #activate the environment
conda install -c bioconda raven-assembler # install the program

#--------------------------------------------------------------------------------

# Move to the directory where the output files will be saved
cd /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/

#input
infa=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/ONT_fltlong/Nuc.3500.fltlong.fastq

raven $infa

conda deactivate