#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=filter_vcf

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=2:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=5G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will filter the vcf file created in the previous script. 

# Load modules 
module load gcc/13.3.0-xp3epyt
module load vcftools/0.1.16

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed

# VCF is the path to the vcf file created in step 11 (note, must unzip file)
VCF=$WORKING_FOLDER/fastq_to_vcf/vcf_freebayes/N_can_pops.vcf

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/fastq_to_vcf

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "vcf_clean" ]
then echo "Working vcf_clean folder exist"; echo "Let's move on."; date
else echo "Working vcf_clean folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/fastq_to_vcf/vcf_clean; date
fi

#--------------------------------------------------------------------------------

# Filter vcf
vcftools --vcf $VCF --minQ 20 --recode --recode-INFO-all --out $WORKING_FOLDER/fastq_to_vcf/vcf_clean/N_can_pops_output_snps-only.vcf

# --minQ: Includes only sites with Quality value above this threshold.
# --recode --recode-INFO-all: These options can be used with the above recode options to define an INFO key name to keep in the output file. This option can be used multiple times to keep more of the INFO fields. The second option is used to keep all INFO values in the original file

