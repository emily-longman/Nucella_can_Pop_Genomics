#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=braker_array

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20

# Reserve walltime -- hh:mm:ss --7 day limit 
#SBATCH --time=1-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=30G

# Submit job array
#SBATCH --array=1 #1-631%30

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will annotate the soft masked genome. However, given the large size of the genome, the length of time exceeds run limits on the VACC. 
# Thus, to make the process more manageable, the genome is broken up into individual scaffolds then braker is run on each scaffold.
# To do this utilize both an array and a while loop. 
# The previous script produced a guide file that groups scaffolds into 1900 scaffold chunks, for a total of 10 partitions.
# For each partition, this script will loop over each scaffold name, break the genome and bam file into that scaffold then clean that scaffold.

# Load modules
module load singularity/3.7.1
spack load samtools@1.10

#--------------------------------------------------------------------------------

#Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly

#This is the location of the soft-masked reference genome.
REFERENCE=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly/Base_Genome_Softmask/N.canaliculata_assembly.fasta.softmasked

#Working folder is core folder where this pipeline is being run.
SCRIPTS_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/00_Genome_short_read

# This is the location of the cleaned and indexed bams
BAM=$WORKING_FOLDER_SCRATCH/cDNA_bam/Nucella.cDNA.srt.rmdp.bam

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

if [ -d "braker_array" ]
then echo "Working braker_array folder exist"; echo "Let's move on."; date
else echo "Working braker_array folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/braker/braker_array; date
fi

#--------------------------------------------------------------------------------

## Import master partition file 
guide_file=$WORKING_FOLDER_SCRATCH/scaffold_names_array.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers. (dimensions: 2, 18919; 10 partitions)
# Scaffold name       # Partition
# Backbone_10001              1
# Backbone_10003              1
# Backbone_10004              1
# Backbone_10005              1
# Backbone_10006              1
# ....

#--------------------------------------------------------------------------------

# Determine partition to process 

# Echo slurm array task ID
echo ${SLURM_ARRAY_TASK_ID}

# Using the guide file, extract the scaffold names associated based on the Slurm array task ID for a given partition
awk '$2=='${SLURM_ARRAY_TASK_ID}'' $guide_file | awk '{print $1}' > partition.names.${SLURM_ARRAY_TASK_ID}.txt

#--------------------------------------------------------------------------------

# Change directory 
cd $WORKING_FOLDER_SCRATCH/braker/braker_array

# For the contigs in a given partition, generate a bam file and genome segement for each chunk then run braker on that chunk. 

$content='cat partition.names.${SLURM_ARRAY_TASK_ID}.txt'

# Break bam file based on the paritions
samtools view -bh ${BAM} $(cat partition.names.${SLURM_ARRAY_TASK_ID}.txt) | tr "\n" " " > output.bam

# Cat file of scaffold names and start while loop
cat partition.names.${SLURM_ARRAY_TASK_ID}.txt | \
while read scaffold 
do echo ${scaffold}


# Break up the bam file into each scaffold
samtools view -b ${BAM} ${scaffold} > ${scaffold}.bam
# Index the bam file
samtools index ${scaffold}.bam

# Break up the genome into each scaffold
grep -EA 1 "^>${scaffold}$" ${REFERENCE} > ${scaffold}.fasta

# Execute breaker 
braker.pl \
--species=Nucella_canaliculata_cDNA_array \
--genome=$REFERENCE \
--threads 20 \
--bam=$WORKING_FOLDER_SCRATCH/cDNA_bam/Nucella.cDNA.srt.rmdp.bam

# Use pilon to polish the genome 
java -Xmx49G -jar $PILONJAR \
--genome ${scaffold}.fasta \
--frags ${scaffold}.bam \
--diploid \
--output ${scaffold}.polished \
--outdir $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_5
# --frags for paired-end sequencing of DNA fragments, such as Illumina paired-end reads of fragment size <1000bp.

# Note: the output of pilon is interleaved (i.e., the DNA sequence is in chunks of 80bps a line), 
# you need to switch them to put all of the sequences on one line 

# Housekeeping - remove intermediate files
rm ${scaffold}.bam
rm ${scaffold}.bam.bai
rm ${scaffold}.fasta
rm ${scaffold}.polished.fasta

done