#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Clean_bams

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=24:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=60G 

# Request CPU
#SBATCH --cpus-per-task=4

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will clean the bam files, specifically it will filter, sort, remove duplicates and index. 
# I will also conduct an intermediary QC step with Qualimap. 

# Load modules  
spack load samtools@1.10

PICARD=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar
qualimap=/netfiles/nunezlab/Shared_Resources/Software/qualimap_v2.2.1/qualimap

#--------------------------------------------------------------------------------

#Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=$WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_2/polished_assembly.fasta

#--------------------------------------------------------------------------------

# Define parameters
QUAL=40 # Quality threshold for samtools
JAVAMEM=18G # Java memory

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_3

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "bams_clean" ]
then echo "Working bams_clean folder exist"; echo "Let's move on."; date
else echo "Working bams_clean folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_3/bams_clean; date
fi

if [ -d "bams_qualimap" ]
then echo "Working bams_qualimap folder exist"; echo "Let's move on."; date
else echo "Working bams_qualimap folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_3/bams_qualimap; date
fi

#--------------------------------------------------------------------------------

# Clean the bam files

# These steps will sort the bams, and duplicates will be removed. I will also conduct an intermediary QC step with Qualimap. 
# Remember to take a look at the qualimap and the flagstat outputs to check for inconsistencies.

# Start pipeline

# Move to working directory
cd $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_3

# Filter merged bam files with samtools view and add flags
samtools view \
-b \
-q $QUAL \
-f 0x0002 -F 0x0004 -F 0x0008 \
--threads 4  \
$WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_3/bams/Ncan.bam \
> $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_3/bams_clean/Ncan.bam
# -q = Skip alignments with MAPQ smaller than $QUAL (40)
# 0x0002 = read mapped in proper pair (0x2)*
# 0x0004 = read unmapped (0x4)
# 0x0008 = mate unmapped (0x8)*

# Sort with picard
# Notice that once a file has been sorted it is added the "srt" suffix
java -Xmx$JAVAMEM -jar $PICARD SortSam \
I=$WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_3/bams_clean/Ncan.bam \
O=$WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_3/bams_clean/Ncan.srt.bam \
SO=coordinate \
VALIDATION_STRINGENCY=SILENT

# Remove duplicates with picard
# Notice that once a file has duplicates removed it is added the "rmdp" suffix
java -Xmx$JAVAMEM -jar $PICARD MarkDuplicates \
I=$WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_3/bams_clean/Ncan.srt.bam \
O=$WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_3/bams_clean/Ncan.srt.rmdp.bam \
M=$WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_3/bams_clean/Ncan.dupstat.txt \
VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

# Index with samtools
samtools index $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_3/bams_clean/Ncan.srt.rmdp.bam

# Lets do QC on the bam file
$qualimap bamqc \
-bam $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_3/bams_clean/Ncan.srt.rmdp.bam \
-outdir $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_3/bams_qualimap/Qualimap_Ncan \
--java-mem-size=$JAVAMEM

# Clean intermediate files
rm $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_3/bams_clean/Ncan.bam
rm $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_3/bams_clean/Ncan.srt.bam

# Housekeeping
mv $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_3/bams_clean/Ncan.dupstat.txt \
$WORKING_FOLDER_SCRATCH/pilon/mapping_stats

#--------------------------------------------------------------------------------

# Inform that sample is done

echo "pipeline completed" $(date)