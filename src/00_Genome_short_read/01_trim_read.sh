#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=trim_genome_short_reads 

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --cpus-per-task=40 
#SBATCH --nodes=1 # on one node

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=24:00:00 #<= this may depend on your resources

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=64G #<= this may depend on your resources

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# This script cleans the raw data using fastp according to parameters set below

# cd to the location (path) to the fastq data:
cd /netfiles/pespenilab_share/Nucella/raw/Shortreads_genome_2024/Project_ESEL_NC3

# Define the sample code to anlayze
READ1=NC3_R1.fastq.gz
READ2=NC3_R2.fastq.gz

NAME1=$(echo $READ1 | sed "s/_R1/_R1_clean/g")
NAME2=$(echo $READ2 | sed "s/_R2/_R2_clean/g")

# print the input and output to screen 

echo $READ1 $READ2
echo $NAME1 $NAME2

# call fastp
/gpfs1/home/e/l/elongman/software/fastp -i ${READ1} -I ${READ2} \
-o /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly/fastp/${NAME1} \
-O /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly/fastp/${NAME2} \
--detect_adapter_for_pe \
--trim_front1 8 \
--trim_poly_g \
--thread 16 \
--cut_right \
--cut_window_size 6 \
--qualified_quality_phred 20 \
--length_required 35 \
--html /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly/fastp/NC3_fastp.html \
--json /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly/fastp/NC3_fastp.json
