#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Map_reads

# Specify partition
#SBATCH --partition=bigmem

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

# Load modules  
spack load gcc@9.3.0
spack load samtools@1.10

bwa=/netfiles/nunezlab/Shared_Resources/Software/bwa-mem2-2.2.1_x64-linux/bwa-mem2.avx2

#--------------------------------------------------------------------------------

#Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly

# This is the location where the reference genome from the first round of pilon and all its indexes are stored.
REFERENCE=$WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_3/polished_assembly.fasta

#--------------------------------------------------------------------------------

# Define parameters
QUAL=40 # Quality threshold for samtools

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to polish directory
cd $WORKING_FOLDER_SCRATCH/pilon

if [ -d "polished_genome_round_4" ]
then echo "Working polished_genome_round_4 folder exist"; echo "Let's move on."; date
else echo "Working polished_genome_round_4 folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_4; date
fi

# Move to polish directory
cd $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_4

if [ -d "sams" ]
then echo "Working sams folder exist"; echo "Let's move on."; date
else echo "Working sams folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_4/sams; date
fi

if [ -d "mapping_stats" ]
then echo "Working mapping_stats folder exist"; echo "Let's move on."; date
else echo "Working mapping_stats folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_4/mapping_stats; date
fi

if [ -d "bams" ]
then echo "Working bams folder exist"; echo "Let's move on."; date
else echo "Working bams folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_4/bams; date
fi

#--------------------------------------------------------------------------------

# Map reads to a reference

# This part will map reads to the reference genome. After reads have been mapped, they will be compressed into bam files, 
# sorted, and duplicates will be removed. I will also conduct an intermediary QC step with Qualimap. 
# Remember to take a look at the qualimap and the flagstat outputs to check for inconsistencies.

# Start pipeline

# Move to working directory
cd $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_4

# Starting mapping
echo "Begin mapping" 
  
# I will conduct the mapping with BWA-MEM 2	
$bwa mem -M -t 24 $REFERENCE \
$WORKING_FOLDER_SCRATCH/fastp/NC3_R1_clean.fastq.gz \
$WORKING_FOLDER_SCRATCH/fastp/NC3_R2_clean.fastq.gz \
> $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_4/sams/Ncan.sam

#--------------------------------------------------------------------------------

# I will now extract some summary stats
samtools flagstat --threads 24 \
$WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_4/sams/Ncan.sam \
> $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_4/mapping_stats/Ncan.flagstats_raw.sam.txt

# Build bam files
samtools view -b -q $QUAL --threads 24  \
$WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_4/sams/Ncan.sam \
> $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_4/bams/Ncan.bam

#--------------------------------------------------------------------------------

# Inform that sample is done

echo "pipeline completed" $(date)