#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=map_CDNA

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 # on one node
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss 
#SBATCH --time=7-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=60G

# Request CPU
#SBATCH --cpus-per-task=6

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# Map cDNA to reference genome using minimap2.

#--------------------------------------------------------------------------------

# Call package (https://github.com/lh3/minimap2)
minimap2=/gpfs1/home/e/l/elongman/software/minimap2-2.28_x64-linux/minimap2

#--------------------------------------------------------------------------------

#Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly

# This is the location where the reference genome and all its indexes are stored.
REFERENCE=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly/Base_Genome_Softmask/N.canaliculata_assembly.fasta.softmasked

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER_SCRATCH

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "cDNA_sam" ]
then echo "Working cDNA_sam folder exist"; echo "Let's move on."; date
else echo "Working cDNA_sam folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/cDNA_sam; date
fi

#--------------------------------------------------------------------------------

# Index and map cDNA to genome w/ minimap2
$minimap2 -ax map-ont \
$REFERENCE \
$WORKING_FOLDER_SCRATCH/cDNA_trim/Nucella.ONT.cDNA.barcode12_clean.fastq.gz \
> $WORKING_FOLDER_SCRATCH/cDNA_sam/Nucella.cDNA.sam