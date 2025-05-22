#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=pilon_polish

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 # on one node
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss 
#SBATCH --time=2:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=40G

# Request CPU
#SBATCH --cpus-per-task=2

# Submit job array
#SBATCH --array=1-634%30

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# This script will polish the genome. However, given the large size of the genome and a large amount of memory is needed to run pilon. 
# To make the process more manageable, break up both the genome and bam file of the short reads into individual scaffolds then polish each piece.
# To do this utilize both an array and a while loop. 
# The previous script produced a guide file that groups scaffolds into 30 scaffold chunks, for a total of 634 partitions.
# For each partition, this script will loop over each scaffold name, break the genome and bam file into that scaffold then clean that scaffold.

# Call packages
spack load samtools@1.10
PILONJAR=/gpfs1/home/e/l/elongman/software/pilon-1.24.jar
seqtk=/gpfs1/home/e/l/elongman/software/seqtk/seqtk

#--------------------------------------------------------------------------------

#Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=$WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_2/polished_assembly.fasta

# This is the location of the cleaned and indexed bams
BAM=$WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_3/bams_clean/Ncan.srt.rmdp.bam

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_3

if [ -d "scaffolds" ]
then echo "Working scaffolds folder exist"; echo "Let's move on."; date
else echo "Working scaffolds folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_3/scaffolds; date
fi

#--------------------------------------------------------------------------------

## Import master partition file 
guide_file=$WORKING_FOLDER_SCRATCH/pilon/guide_file_array_3.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers. (dimensions: 2, 19014; 634 partitions)
# Scaffold name       # Partition
# ntLink-6475              1
# ntLink-6476              1
# ntLink-6477              1
# ntLink-6478              1
# ntLink-6479              1
# ....

#--------------------------------------------------------------------------------

# Determine partition to process 

# Echo slurm array task ID
echo ${SLURM_ARRAY_TASK_ID}

# Using the guide file, extract the scaffold names associated based on the Slurm array task ID for a given partition
awk '$2=='${SLURM_ARRAY_TASK_ID}'' $guide_file | awk '{print $1}' > ntLink.names.${SLURM_ARRAY_TASK_ID}.txt

#--------------------------------------------------------------------------------

# For the scaffolds in a given partition, generate a bam file and genome segement for each scaffold then polish each piece. 

# Cat file of scaffold names and start while loop
cat ntLink.names.${SLURM_ARRAY_TASK_ID}.txt | \
while read scaffold 
do echo ${scaffold}

# Break up the bam file into each scaffold
samtools view -b ${BAM} ${scaffold} > ${scaffold}.bam
# Index the bam file
samtools index ${scaffold}.bam

# Break up the genome into each scaffold
grep -EA 1 "^>${scaffold}$" ${REFERENCE} > ${scaffold}.fasta

# Use pilon to polish the genome 
java -Xmx49G -jar $PILONJAR \
--genome ${scaffold}.fasta \
--frags ${scaffold}.bam \
--diploid \
--output ${scaffold}.polished \
--outdir $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_3
# --frags for paired-end sequencing of DNA fragments, such as Illumina paired-end reads of fragment size <1000bp.

# Note: the output of pilon is interleaved (i.e., the DNA sequence is in chunks of 80bps a line), 
# you need to switch them to put all of the sequences on one line 

$seqtk seq ${scaffold}.polished.fasta > $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_3/scaffolds/${scaffold}.polished.fasta

# Housekeeping - remove intermediate files
rm ${scaffold}.bam
rm ${scaffold}.bam.bai
rm ${scaffold}.fasta
rm ${scaffold}.polished.fasta

done

#--------------------------------------------------------------------------------

# Housekeeping
rm ntLink.names.${SLURM_ARRAY_TASK_ID}.txt