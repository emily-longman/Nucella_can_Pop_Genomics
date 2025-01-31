#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=prep_bam

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=20:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=60G 

# Request CPU
#SBATCH --cpus-per-task=5

# Submit job array
#SBATCH --array=1-19

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/Bam_sync.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will prep the bam files. More specifically, it will add read group information and index the final bam files. 

# Load modules 
module load gcc/13.3.0-xp3epyt
module load samtools/1.19.2-pfmpoam

### program dependencies
PICARD=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed/fastq_to_bam

#--------------------------------------------------------------------------------

# Define parameters
JAVAMEM=18G
echo ${SLURM_ARRAY_TASK_ID}

#--------------------------------------------------------------------------------

# Read guide files
# This is a file with the name all the samples to be processed. One sample name per line with all the info.
GUIDE_FILE=$WORKING_FOLDER/guide_files/Merge_bams.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers.
## Population      Merged_name 1      Merged_name 2 
##   ARA	       ARA_S168_L006	  ARA_S13_L008
##   BMR	       BMR_S156_L006	  BMR_S1_L008
##   CBL	       CBL_S169_L006	  CBL_S14_L008
##   ...
##   VD	           VD_S161_L006	      VD_S6_L008

#--------------------------------------------------------------------------------

# Determine sample to process, "i" and read files
i=`awk -F "\t" '{print $1}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo ${i}

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "RGSM_final_bams" ]
then echo "Working RGSM_final_bams folder exist"; echo "Let's move on."; date
else echo "Working RGSM_final_bams folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/RGSM_final_bams; date
fi

#--------------------------------------------------------------------------------

# Read Information
Group_library="Longman_2024"
Library_platform="NovaSeq"
Group_platform="EKL2024"

# Force a uniform read group to the joint bam file
java -jar $PICARD AddOrReplaceReadGroups \
I=$WORKING_FOLDER/bams_merged/${i}.lanes_merged.bam \
O=$WORKING_FOLDER/RGSM_final_bams/${i}.bam \
RGID=${Group_library}.${i} \
RGLB=$Group_library \
RGPL=$Library_platform \
RGPU=$Group_platform \
RGSM=${i}

#--------------------------------------------------------------------------------

# Index bams with samtools
samtools index $WORKING_FOLDER/RGSM_final_bams/${i}.bam

#--------------------------------------------------------------------------------

date
echo "done"