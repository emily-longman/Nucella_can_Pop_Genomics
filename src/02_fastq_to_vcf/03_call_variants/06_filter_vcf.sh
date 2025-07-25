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
#SBATCH --time=5:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=5G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will filter the vcf file created in the previous script. 

#--------------------------------------------------------------------------------

# Load modules 
module load gcc/13.3.0-xp3epyt
module load vcftools/0.1.16

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/fastq_to_vcf

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "vcf_clean" ]
then echo "Working vcf_clean folder exist"; echo "Let's move on."; date
else echo "Working vcf_clean folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/fastq_to_vcf/vcf_clean; date
fi

#--------------------------------------------------------------------------------

# Filter vcf

# Filter based on missing data (keep only variants in 75% of pops)
vcftools --gzvcf $WORKING_FOLDER/data/processed/fastq_to_vcf/vcf_freebayes/N.canaliculata_pops.vcf.gz \
--minQ 60 --maf 0.01 --max-missing 1.0 --remove-indels --recode --recode-INFO-all \
--out $WORKING_FOLDER/data/processed/fastq_to_vcf/vcf_clean/N.canaliculata_pops_filter_minQ60_maxmissing1.0

# --minQ: Includes only sites with Quality value above this threshold.
# --maf: Include only sites with a minor allele frequency great than or equal to this value.
# --max-missing: Exclude sites on the basis of the proportion of missing data (defined btween 0 and 1). 0 allows sites that are completely missing and 1 indicates no missing data.
# --recode --recode-INFO-all: Write a new VCF file that keeps all the INFO flags from the old vcf file.

#--------------------------------------------------------------------------------

# Check the quality of output vcf

# Generate a summary of the number of SNPs for each filter category
vcftools --vcf $WORKING_FOLDER/data/processed/fastq_to_vcf/vcf_clean/N.canaliculata_pops_filter_minQ60_maxmissing1.0.recode.vcf \
--FILTER-summary --out $WORKING_FOLDER/data/processed/fastq_to_vcf/vcf_clean/N.canaliculata_pops_filter_minQ60_maxmissing1.0.recode.vcf
# Generate a file containing the depth per site summed across all individuals
vcftools --vcf $WORKING_FOLDER/data/processed/fastq_to_vcf/vcf_clean/N.canaliculata_pops_filter_minQ60_maxmissing1.0.recode.vcf \
--depth --out $WORKING_FOLDER/data/processed/fastq_to_vcf/vcf_clean/N.canaliculata_pops_filter_minQ60_maxmissing1.0.recode.vcf