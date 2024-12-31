#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Trim_reads

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=6

# Reserve walltime -- hh:mm:ss
#SBATCH --time=1:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=10G 

# Submit job array
#SBATCH --array=1-38%20

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/Trim_reads.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will trim the reads.

# Load modules  
fastp=/gpfs1/home/e/l/elongman/software/fastp

#--------------------------------------------------------------------------------

# Define important file locations

# RAW_READS indicates the folder where the raw reads are stored.
RAW_READS=/netfiles/pespenilab_share/Nucella/raw/Population_genomics/All_shortreads

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed/fastq_to_bam

# Name of pipeline
PIPELINE=Trim_reads

#--------------------------------------------------------------------------------

# Read guide files
# This is a file with the name all the samples to be processed. One sample name per line with all the info.
GUIDE_FILE=$WORKING_FOLDER/guide_files/Trim_map.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers.
##             Read 1                            Read 2             Population   Sample#   Lane#    Paired_name    
## ARA_S168_L006_R1_001.fastq.gz	ARA_S168_L006_R2_001.fastq.gz	    ARA 	   S168	   L006	   ARA_S168_L006
## BMR_S156_L006_R1_001.fastq.gz	BMR_S156_L006_R2_001.fastq.gz	    BMR	       S156    L006	   BMR_S156_L006
## CBL_S169_L006_R1_001.fastq.gz	CBL_S169_L006_R2_001.fastq.gz	    CBL	       S169	   L006	   CBL_S169_L006
## ...
## VD_S6_L008_R1_001.fastq.gz	    VD_S6_L008_R2_001.fastq.gz	        VD	        S6	   L008	    VD_S6_L008

#--------------------------------------------------------------------------------

# Determine sample to process, "i" and read files
i=`awk -F "\t" '{print $6}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
read1=`awk -F "\t" '{print $1}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
read2=`awk -F "\t" '{print $2}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo "Sample i:" ${i} "Read 1:" ${read1} "Read 2:" ${read2}

#--------------------------------------------------------------------------------

# This part of the pipeline will generate log files to record warnings and completion status

# Move to working directory
cd $WORKING_FOLDER

# This part of the script will check and generate, if necessary, all of the output folders used in the script

# Create logs directory
if [ -d "logs" ]
then echo "Working logs folder exist"; echo "Let's move on."; date
else echo "Working logs folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/logs; date
fi

# Move to logs directory
cd $WORKING_FOLDER/logs

echo $PIPELINE

if [[ -e "${PIPELINE}.completion.log" ]]
then echo "Completion log exist"; echo "Let's move on."; date
else echo "Completion log doesnt exist. Let's fix that."; touch $WORKING_FOLDER/logs/${PIPELINE}.completion.log; date
fi

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "trimmed_reads" ]
then echo "Working trimmed_reads folder exist"; echo "Let's move on."; date
else echo "Working trimmed_reads folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/trimmed_reads; date
fi

if [ -d "trimmed_reads_reports" ]
then echo "Working trimmed_reads_reports folder exist"; echo "Let's move on."; date
else echo "Working trimmed_reads_reports folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/trimmed_reads_reports; date
fi

#--------------------------------------------------------------------------------

# Start pipeline

# Move to working directory
cd $WORKING_FOLDER

# This script uses an array and matches left and right reads and cleans the raw data using fastp according to parameters set below
# Decide on the trimming parameters based on fastQC step done before this script.

echo "Trimming reads for sample:" ${i}

# Call fastp and do some light trimming
$fastp \
-i $RAW_READS/${read1} \
-I $RAW_READS/${read2} \
-o $WORKING_FOLDER/trimmed_reads/${i}_R1_clean.fq.gz \
-O $WORKING_FOLDER/trimmed_reads/${i}_R2_clean.fq.gz \
--detect_adapter_for_pe \
--trim_front1 8 \
--trim_poly_g \
--thread 6 \
--cut_right \
--cut_right_window_size 6 \
--qualified_quality_phred 20 \
--html $WORKING_FOLDER/trimmed_reads_reports/${i}_clean.html \
--json $WORKING_FOLDER/trimmed_reads_reports/${i}_clean.json

# i = read 1
# I = read 2
# o & O = outputs (.json and html for qc)

# For PE data you can specify adapter sequence auto-detection by specifying --detect_adapter_for_pe
# For PE data the front trimming settings are --trim_front1 (You can also do trimming on each read separately by specifying trim_front2; if not specified then trim_front2=trim_front1)
# Detect and trim polyG in read tails using --trim_poly_g
# Per read cutting by quality options:
### --cut_right: Move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, then stop
### --cut_right_window_size: The window size option for cut_right
### --qualified_quality_phred: The quality vlue that a base is qualified 

#--------------------------------------------------------------------------------
# Inform that sample is done

# This part of the pipeline will notify the completion of run i. 

echo ${i} " completed" >> $WORKING_FOLDER/logs/${PIPELINE}.completion.log

echo "pipeline completed" $(date)