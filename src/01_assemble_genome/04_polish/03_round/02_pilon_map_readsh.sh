#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Map_reads

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss 
#SBATCH --time=24:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=50G 

# Request CPU
#SBATCH --cpus-per-task=24

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will map the reads with bwa mem.

#--------------------------------------------------------------------------------

# Load modules  
module load gcc/13.3.0-xp3epyt
module load samtools/1.19.2-pfmpoam
bwa=/netfiles/nunezlab/Shared_Resources/Software/bwa-mem2-2.2.1_x64-linux/bwa-mem2.avx2

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

# This is the location where the reference genome from the first round of pilon and all its indexes are stored.
REFERENCE=$WORKING_FOLDER/data/processed/genome_assembly/pilon/polished_genome_round_2/polished_assembly.fasta

#--------------------------------------------------------------------------------

# Define parameters
QUAL=40 # Quality threshold for samtools

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to polish directory
cd $WORKING_FOLDER/data/processed/genome_assembly/pilon

if [ -d "polished_genome_round_3" ]
then echo "Working polished_genome_round_3 folder exist"; echo "Let's move on."; date
else echo "Working polished_genome_round_3 folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/genome_assembly/pilon/polished_genome_round_3; date
fi

# Move to polish directory
cd $WORKING_FOLDER/data/processed/genome_assembly/pilon/polished_genome_round_3

if [ -d "sams" ]
then echo "Working sams folder exist"; echo "Let's move on."; date
else echo "Working sams folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/genome_assembly/pilon/polished_genome_round_3/sams; date
fi

if [ -d "mapping_stats" ]
then echo "Working mapping_stats folder exist"; echo "Let's move on."; date
else echo "Working mapping_stats folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/genome_assembly/pilon/polished_genome_round_3/mapping_stats; date
fi

if [ -d "bams" ]
then echo "Working bams folder exist"; echo "Let's move on."; date
else echo "Working bams folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/genome_assembly/pilon/polished_genome_round_3/bams; date
fi

#--------------------------------------------------------------------------------

# Map reads to a reference

# This part will map reads to the reference genome. After reads have been mapped, they will be compressed into bam files, 
# sorted, and duplicates will be removed. I will also conduct an intermediary QC step with Qualimap. 
# Remember to take a look at the qualimap and the flagstat outputs to check for inconsistencies.

# Start pipeline

# Move to working directory
cd $WORKING_FOLDER/data/processed/genome_assembly/pilon/polished_genome_round_3

# Starting mapping
echo "Begin mapping" 
  
# I will conduct the mapping with BWA-MEM 2	
$bwa mem -M -t 24 $REFERENCE \
$WORKING_FOLDER/data/processed/genome_assembly/trim_AVITI/NC3_R1_clean.fq.gz \
$WORKING_FOLDER/data/processed/genome_assembly/trim_AVITI/NC3_R2_clean.fq.gz \
> $WORKING_FOLDER/data/processed/genome_assembly/pilon/polished_genome_round_3/sams/Ncan.sam

#--------------------------------------------------------------------------------

# I will now extract some summary stats
samtools flagstat --threads 24 \
$WORKING_FOLDER/data/processed/genome_assembly/pilon/polished_genome_round_3/sams/Ncan.sam \
> $WORKING_FOLDER/data/processed/genome_assembly/pilon/polished_genome_round_3/mapping_stats/Ncan.flagstats_raw.sam.txt

# Build bam files
samtools view -b -q $QUAL --threads 24  \
$WORKING_FOLDER/data/processed/genome_assembly/pilon/polished_genome_round_3/sams/Ncan.sam \
> $WORKING_FOLDER/data/processed/genome_assembly/pilon/polished_genome_round_3/bams/Ncan.bam

#--------------------------------------------------------------------------------

# Inform that sample is done

echo "pipeline completed" $(date)