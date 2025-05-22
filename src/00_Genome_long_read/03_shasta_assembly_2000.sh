#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=shasta 

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --cpus-per-task=1 
#SBATCH --nodes=1 # on one node

# Reserve walltime -- hh:mm:ss
#SBATCH --time=24:00:00 #<= this may depend on your resources

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=900G #<= this may depend on your resources


# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/shasta.%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

#This script was because I realized I skipped one of them in my array.

# Move to the directory where the output files will be saved
cd /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Shasta_assemblies/

#executable
shasta=/netfiles/nunezlab/Shared_Resources/Software/shasta/shasta-Linux-0.10.0

#input
infa=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/ONT_fltlong
conf=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/00_genome/shasta_config

#run shasta
$shasta --input $infa/Nuc.2000.fltlong.fastq --config $conf/Nanopore-May2022_2000.sh --assemblyDirectory ShastaRun2000