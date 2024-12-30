#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Merge_bams

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=1:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=10G 

# Submit job array
#SBATCH --array=1-19%10

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/Merge_bams.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will gather all data for each population across the 2 lanes of sequencing.

# Load modules  
spack load samtools@1.10
qualimap=/netfiles/nunezlab/Shared_Resources/Software/qualimap_v2.2.1/qualimap

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed/fastq_to_bam

# Name of pipeline
PIPELINE=Merge_bams

#--------------------------------------------------------------------------------

# Define parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
QUAL=40 # Quality threshold for samtools
JAVAMEM=18G # Java memory

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
echo $i

#--------------------------------------------------------------------------------

# This part of the pipeline will generate log files to record warnings and completion status

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

if [ -d "bams_merged" ]
then echo "Working bams_merged folder exist"; echo "Let's move on."; date
else echo "Working bams_merged folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/bams_merged; date
fi

if [ -d "bams_merged_qualimap" ]
then echo "Working bams_merged_qualimap folder exist"; echo "Let's move on."; date
else echo "Working bams_merged_qualimap folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/bams_merged_qualimap; date
fi

#--------------------------------------------------------------------------------

# Here I will merge the bam outputs for the 2 lanes of sequencing. These will be named 'Lanes merged'

echo "I will merge these files" $WORKING_FOLDER/bams_clean/${i}_*.srt.rmdp.bam

# Make temporary linefile with list of input BAM files
ls $WORKING_FOLDER/bams_clean/${i}_*.srt.rmdp.bam > ${i}.guide.txt

# Merge the 2 sequencing lanes
samtools merge \
-b ${i}.guide.txt \
$WORKING_FOLDER/bams_merged/${i}.lanes_merged.bam

# Housekeeping: Remove the temporary guide file
rm ${i}.guide.txt

#--------------------------------------------------------------------------------

# Assess quality of the final file
$qualimap bamqc \
-bam $WORKING_FOLDER/bams_merged/${i}.lanes_merged.bam \
-outdir $WORKING_FOLDER/bams_merged_qualimap/Qualimap_LaneMerged_${i} \
--java-mem-size=$JAVAMEM

#--------------------------------------------------------------------------------

# Index bams with samtools
samtools index $WORKING_FOLDER/bams_merged/${i}.lanes_merged.bam

#--------------------------------------------------------------------------------

# Inform that sample is done

# This part of the pipeline will notify the completion of run i. 

echo ${i} " completed" >> $WORKING_FOLDER/logs/${PIPELINE}.completion.log

echo "pipeline completed" $(date)
