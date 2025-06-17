#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=fastqc 

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=5 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=2:00:00

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=5G

# Submit job array
#SBATCH --array=1-76%20

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/FastQC.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu

#--------------------------------------------------------------------------------

# This script will initiate a pipeline which will do some quality QC on the reads. 

# Load modules 
module load gcc/13.3.0-xp3epyt
module load fastqc/0.12.1-qxseug5

#--------------------------------------------------------------------------------

# Define important file locations

# RAW_READS indicates the folder where the raw reads are stored.
RAW_READS=/netfiles/pespenilab_share/Nucella/raw/Population_genomics/All_shortreads

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed

# Name of pipeline
PIPELINE=fastqc

#--------------------------------------------------------------------------------

# Read guide files
# This is a file with the name all the samples to be processed. One sample name per line with all the info.
GUIDE_FILE=$WORKING_FOLDER/fastq_to_vcf/guide_files/QC_reads.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers.
##               Read            Population  Sample#  Lane#  Forward/Reverse
## ARA_S168_L006_R1_001.fastq.gz	ARA	       S168	   L006	        R1
## ARA_S168_L006_R2_001.fastq.gz	ARA	       S168	   L006	        R2
## BMR_S156_L006_R1_001.fastq.gz	BMR	       S156	   L006	        R1
## BMR_S156_L006_R2_001.fastq.gz	BMR        S156	   L006	        R2
## ...
## VD_S6_L008_R1_001.fastq.gz	    VD	        S6	   L008	        R1
## VD_S6_L008_R2_001.fastq.gz	    VD	        S6	   L008	        R2

#--------------------------------------------------------------------------------

# Determine sample to process, "i" and read files
i=`awk -F "\t" '{print $1}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo ${i}

#--------------------------------------------------------------------------------

# Generate Folders and Files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

# Move to working directory
cd $WORKING_FOLDER/fastq_to_vcf

# Generating new folder
if [ -d "qc_reads" ]
then echo "Working qc_reads folder exist"; echo "Let's move on"; date
else echo "Working qc_reads folder doesnt exist. Let's fix that"; mkdir $WORKING_FOLDER/fastq_to_vcf/qc_reads; date
fi

# Change directory
cd $WORKING_FOLDER/fastq_to_vcf/qc_reads

# Generating new folders 
if [ -d "fastqc" ]
then echo "Working fastqc folder exist"; echo "Let's move on"; date
else echo "Working fastqc folder doesnt exist. Let's fix that"; mkdir $WORKING_FOLDER/fastq_to_vcf/qc_reads/fastqc; date
fi

#--------------------------------------------------------------------------------

# Do QC on raw reads with fastqc

# Move to working directory
cd $WORKING_FOLDER

echo -e $i "is now processing"; date

# Lets do some QC on the reads
fastqc $RAW_READS/${i} \
--outdir $WORKING_FOLDER/fastq_to_vcf/qc_reads/fastqc

#--------------------------------------------------------------------------------

# Inform that sample is done

# This part of the pipeline will produce a notification stating the completion of the script. 

echo ${i} " completed"

echo "pipeline" ${PIPELINE} $(date)