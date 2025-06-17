#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=sam_to_bam_cDNA

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 # on one node
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss 
#SBATCH --time=2-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=5G

# Request CPU
#SBATCH --cpus-per-task=6

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# Convert sam to bam for cDNA.

#--------------------------------------------------------------------------------

# Load modules  
spack load gcc@9.3.0
spack load samtools@1.10

#--------------------------------------------------------------------------------

#Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly

# This is the location where the reference genome and all its indexes are stored.
REFERENCE=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly/Base_Genome_Softmask/N.canaliculata_assembly.fasta.softmasked

#--------------------------------------------------------------------------------

# Define parameters
CPU=6
echo "using #CPUs ==" $CPU
JAVAMEM=18G # Java memory

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER_SCRATCH

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "mapping_stats" ]
then echo "Working mapping_stats folder exist"; echo "Let's move on."; date
else echo "Working mapping_stats folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/mapping_stats; date
fi

if [ -d "cDNA_bam" ]
then echo "Working cDNA_bam folder exist"; echo "Let's move on."; date
else echo "Working cDNA_bam folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/cDNA_bam; date
fi

#--------------------------------------------------------------------------------

# Extract summary stats on sam
samtools flagstat --threads $CPU \
$WORKING_FOLDER_SCRATCH/cDNA_sam/Nucella.cDNA.sam \
> $WORKING_FOLDER_SCRATCH/mapping_stats/Nucella.cDNA.flagstats.sam.txt
# Remember to take a look at the flagstat outputs to check for inconsistencies.

# Build bam files
samtools view -b --threads $CPU  \
$WORKING_FOLDER_SCRATCH/cDNA_sam/Nucella.cDNA.sam \
> $WORKING_FOLDER_SCRATCH/cDNA_bam/Nucella.cDNA.bam