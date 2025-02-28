#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=grenedalf

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=30:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=5G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will calculate population genetics statistics using grenedalf (grenedalf=/gpfs1/home/e/l/elongman/software/grenedalf_v0.6.2_linux_x86_64). 

# Load modules 
grenedalf=/gpfs1/home/e/l/elongman/software/grenedalf_v0.6.2_linux_x86_64

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed

# This is the location where the reference genome is stored.
REFERENCE_FOLDER=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Oct2024/Crassostrea_softmask

# This is the path to the reference genome.
REFERENCE=$REFERENCE_FOLDER/N.canaliculata_assembly.fasta.softmasked.fa

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/pop_gen

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "grenedalf" ]
then echo "Working grenedalf folder exist"; echo "Let's move on."; date
else echo "Working grenedalf folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/pop_gen/grenedalf; date
fi

#--------------------------------------------------------------------------------

# Calculate diversity metric for each population using the final bam files

$grenedalf diversity \
--sam-path $WORKING_FOLDER/fastq_to_vcf/RGSM_final_bams/*.bam \
--pool-sizes $WORKING_FOLDER/pop_gen/guide_files/Grenedalf_pool_size.csv \
--window-type genome \
--window-average-policy valid-loci \
--filter-sample-min-count 2 \
--filter-sample-min-read-depth 20 \
--reference-genome-fasta $REFERENCE \
--out-dir $WORKING_FOLDER/pop_gen/grenedalf \
--no-extra-columns




