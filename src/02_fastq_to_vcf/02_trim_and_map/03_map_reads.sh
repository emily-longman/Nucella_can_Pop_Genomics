#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Map_reads

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss
#SBATCH --time=20:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=20G 

# Request CPU
#SBATCH --cpus-per-task=10

# Submit job array
#SBATCH --array=1-38%19

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will map reads to the masked reference genome using bwa mem. 
# After reads have been mapped, they will be compressed into bam files.

#--------------------------------------------------------------------------------

# Load modules  
module load gcc/13.3.0-xp3epyt
module load samtools/1.19.2-pfmpoam
bwa=/netfiles/nunezlab/Shared_Resources/Software/bwa-mem2-2.2.1_x64-linux/bwa-mem2.avx2

#--------------------------------------------------------------------------------

# Define important file locations

# RAW_READS indicates the folder where the raw reads are stored.
RAW_READS=/netfiles/pespenilab_share/Nucella/raw/Population_genomics/All_shortreads

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

# This is the location where the reference genome is stored.
REFERENCE_FOLDER=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Oct2024/Crassostrea_softmask

# This is the path to the reference genome.
REFERENCE=$REFERENCE_FOLDER/N.canaliculata_assembly.fasta.softmasked.fa

# Name of pipeline
PIPELINE=Map_reads

#--------------------------------------------------------------------------------

# Define parameters
CPU=10
echo "using #CPUs ==" $CPU
QUAL=40 # Quality threshold for samtools
JAVAMEM=18G # Java memory

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
echo ${i}

#--------------------------------------------------------------------------------

# This part of the pipeline will generate log files to record warnings and completion status

# Move to logs directory
cd $WORKING_FOLDER/data/processed/fastq_to_vcf/logs

echo $PIPELINE

if [[ -e "${PIPELINE}.completion.log" ]]
then echo "Completion log exist"; echo "Let's move on."; date
else echo "Completion log doesnt exist. Let's fix that."; touch $WORKING_FOLDER/data/processed/fastq_to_vcf/logs/${PIPELINE}.completion.log; date
fi

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/fastq_to_vcf

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "sams" ]
then echo "Working sams folder exist"; echo "Let's move on."; date
else echo "Working sams folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/fastq_to_vcf/sams; date
fi

if [ -d "mapping_stats" ]
then echo "Working mapping_stats folder exist"; echo "Let's move on."; date
else echo "Working mapping_stats folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/fastq_to_vcf/mapping_stats; date
fi

if [ -d "bams" ]
then echo "Working bams folder exist"; echo "Let's move on."; date
else echo "Working bams folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/fastq_to_vcf/bams; date
fi

#--------------------------------------------------------------------------------

# Map reads to a reference

# Move to working directory
cd $WORKING_FOLDER/data/processed/fastq_to_vcf

# Starting mapping
echo "Begin mapping" ${i}
  
# I will conduct the mapping with BWA-MEM 2	
$bwa mem -M -t $CPU $REFERENCE \
$WORKING_FOLDER/data/processed/fastq_to_vcf/trimmed_reads/${i}_R1_clean.fq.gz \
$WORKING_FOLDER/data/processed/fastq_to_vcf/trimmed_reads/${i}_R2_clean.fq.gz \
> $WORKING_FOLDER/data/processed/fastq_to_vcf/sams/${i}.sam

#--------------------------------------------------------------------------------

# Extract sam summary stats
samtools flagstat --threads $CPU \
$WORKING_FOLDER/data/processed/fastq_to_vcf/sams/${i}.sam \
> $WORKING_FOLDER/data/processed/fastq_to_vcf/mapping_stats/${i}.flagstats_raw.sam.txt
# Take a look at the flagstat outputs to check for inconsistencies.

#--------------------------------------------------------------------------------

# Build bam files
samtools view -b -q $QUAL --threads $CPU  \
$WORKING_FOLDER/data/processed/fastq_to_vcf/sams/${i}.sam \
> $WORKING_FOLDER/data/processed/fastq_to_vcf/bams/${i}.bam

#--------------------------------------------------------------------------------

# Inform that sample is done

# This part of the pipeline will notify the completion of run i. 

echo ${i} "completed" >> $WORKING_FOLDER/data/processed/fastq_to_vcf/logs/${PIPELINE}.completion.log

echo "pipeline completed" $(date)