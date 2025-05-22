#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=braker

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20

# Reserve walltime -- hh:mm:ss --7 day limit 
#SBATCH --time=7-00:00:00 

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

#--------------------------------------------------------------------------------

#Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly

#This is the location of the reference genome.
REFERENCE=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly/Base_Genome_Softmask/N.canaliculata_assembly.fasta.softmasked

#Working folder is core folder where this pipeline is being run.
SCRIPTS_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/00_Genome_short_read

# Export path to braker sif
export BRAKER_SIF=$SCRIPTS_FOLDER/23_braker_singularity/braker3.sif

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER_SCRATCH

if [ -d "braker" ]
then echo "Working braker folder exist"; echo "Let's move on."; date
else echo "Working braker folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/braker; date
fi

# Move to working directory
cd $WORKING_FOLDER_SCRATCH/braker

if [ -d "braker_cDNA" ]
then echo "Working braker_cDNA folder exist"; echo "Let's move on."; date
else echo "Working braker_cDNA folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/braker/braker_cDNA; date
fi

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER_SCRATCH/braker/braker_cDNA

# Execute breaker
braker.pl \
--species=Nucella_canaliculata_cDNA_bam \
--genome=$REFERENCE \
--threads 20 \
--bam=$WORKING_FOLDER_SCRATCH/cDNA_bam/Nucella.cDNA.srt.rmdp.bam