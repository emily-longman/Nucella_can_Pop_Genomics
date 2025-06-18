#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=cat.reads

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --ntasks-per-node=1 

# Reserve walltime -- hh:mm:ss
#SBATCH --time=5:00:00

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=80G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out 

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# Concatenate ONT reads; all pass filter

#--------------------------------------------------------------------------------

# Define important file locations

# RAW_DATA is the core folder where all of the raw data is located.
RAW_DATA=/netfiles/pespenilab_share/Nucella/raw

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/raw

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "ONT" ]
then echo "Working ONT folder exist"; echo "Let's move on."; date
else echo "Working ONT folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/raw/ONT; date
fi

#--------------------------------------------------------------------------------

# Concatenate reads from each ONT flow cell
cat $RAW_DATA/ONT/PROM0150_Sanford_NC3-PCR_FC2_11152023/20231115_1701_2E_PAS52172_fe64f551/fastq_pass/*.fastq.gz > $WORKING_FOLDER/data/raw/ONT/FC1_ONT_NC3.fastq.gz

cat $RAW_DATA/ONT/PROM0150_Sanford_NC3-PCR_05312023/20230531_1740_3G_PAQ05812_93506d2a/fastq_pass/*.fastq.gz > $WORKING_FOLDER/data/raw/ONT/FC2_ONT_NC3.fastq.gz

cat $RAW_DATA/ONT/PROM0172_Sanford_NC3-PCR_FC1_02282024/20240228_1655_2E_PAW14761_69937e95/fastq_pass/*.fastq.gz > $WORKING_FOLDER/data/raw/ONT/FC3_ONT_NC3.fastq.gz

cat $RAW_DATA/ONT/PROM0172_Sanford_NC3-PCR_FC2_02282024/20240228_1655_3G_PAU83574_150a145e/fastq_pass/*.fastq.gz > $WORKING_FOLDER/data/raw/ONT/FC4_ONT_NC3.fastq.gz

cat $RAW_DATA/ONT/PROM0172_Sanford_NC3-PCR_FC3_03062024/20240306_1657_1G_PAU75920_614ca090/fastq_pass/*.fastq.gz > $WORKING_FOLDER/data/raw/ONT/FC5_ONT_NC3.fastq.gz

#--------------------------------------------------------------------------------

# Concatenate reads into one fastq.gz file
cat $WORKING_FOLDER/data/raw/ONT/*_ONT_NC3.fastq.gz > $WORKING_FOLDER/data/raw/ONT/FC_all.ONT.nuc.fastq.gz
