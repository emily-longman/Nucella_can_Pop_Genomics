#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=genotype_table

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

# This script will generate a genotype table of the variable loci. 

# Load modules 
module load gcc/13.3.0-xp3epyt
module load bcftools/1.19-iq5mwek

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/fastq_to_vcf

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "genotype_table" ]
then echo "Working genotype_table folder exist"; echo "Let's move on."; date
else echo "Working genotype_table folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/fastq_to_vcf/genotype_table; date
fi

#--------------------------------------------------------------------------------





#output genotype table (variable loci)
bcftools query -f '%CHROM %POS[\t%GT]\n' \
snpcall_all/snpcall.final.rmdup.chr.${chr}.flt3.vcf.gz > snpcall_all/varpos.final.rmdup.chr.${chr}.flt3.genotypes.txt