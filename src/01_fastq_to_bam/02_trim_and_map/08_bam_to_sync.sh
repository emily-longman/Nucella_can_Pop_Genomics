#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=make_sync

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

# This script will convert the bam files to gsync files (genome-wide extension of the SYNC file format) and will mask regions based on minimum and max depth thresholds. 

# Load modules 
module load python3.10-anaconda/2023.03-1
module load openjdk/1.8.0
module load gcc/13.3.0-xp3epyt
module load samtools/1.19.2-pfmpoam@1.10

### program dependencies
PICARD=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar
gatk3=/netfiles/nunezlab/Shared_Resources/Software/gatk3/gatk/GenomeAnalysisTK.jar
tabix=/netfiles/nunezlab/Shared_Resources/Software/htslib/tabix
bgzip=/netfiles/nunezlab/Shared_Resources/Software/htslib/bgzip

# Load Mpileup2sync and 
Mpileup2Sync=/netfiles/nunezlab/Shared_Resources/Software/DESTv2/mappingPipeline/scripts/Mpileup2Sync.py
MaskSYNC_snape_complete=/netfiles/nunezlab/Shared_Resources/Software/DESTv2/mappingPipeline/scripts/MaskSYNC_snape_complete.py

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed/fastq_to_bam

# This is the location where the reference genome is stored.
REFERENCE_FOLDER=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Oct2024/Crassostrea_softmask

# This is the path to the reference genome.
REFERENCE=$REFERENCE_FOLDER/N.canaliculata_assembly.fasta.softmasked.fa

# This is the path to the gff file.
GFF=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/N.can.gff

# This is the path to the pickled genome file.
PICKLED=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/N.canaliculata_assembly.fasta.softmasked_pickled_ref

#--------------------------------------------------------------------------------

# Define parameters
JAVAMEM=59G
THREADS=5
echo ${SLURM_ARRAY_TASK_ID}

# Read Information
Group_library="Longman_2024"
Library_platform="illumina"
Group_platform="EKL2024"

# GATK parameters
base_quality_threshold=25
illumina_quality_coding=1.8
minIndel=5
maxsnape=0.9
max_cov=0.95 
min_cov=10 

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

if [ -d "syncfiles" ]
then echo "Working syncfiles folder exist"; echo "Let's move on."; date
else echo "Working syncfiles folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/syncfiles; date
fi

# Change directory
cd syncfiles

if [ -d "${i}" ]
then echo "Working ${i} folder exist"; echo "Let's move on."; date
else echo "Working ${i} folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/syncfiles/${i}; date
fi

#--------------------------------------------------------------------------------

# Force a uniform read group to the joint bam file
java -jar $PICARD AddOrReplaceReadGroups \
I=$WORKING_FOLDER/bams_merged/${i}.lanes_merged.bam \
O=$WORKING_FOLDER/RGSM_final_bams/${i}.bam \
RGLB=$Group_library \
RGPL=$Library_platform \
RGPU=$Group_platform \
RGSM=${i}

#--------------------------------------------------------------------------------

# Index bams with samtools
samtools index $WORKING_FOLDER/RGSM_final_bams/${i}.bam

#--------------------------------------------------------------------------------

# Use RealignerTargetCreator to identify and create a target intervals list
java -jar $gatk3 -T RealignerTargetCreator \
-nt $THREADS \
-R $REFERENCE \
-I $WORKING_FOLDER/RGSM_final_bams/${i}.bam \
-o $WORKING_FOLDER/syncfiles/${i}/${i}.hologenome.intervals

# Perform local realignment for the target intervals using IndelRealigner 
java -jar $gatk3 -T IndelRealigner \
-R $REFERENCE \
-I $WORKING_FOLDER/RGSM_final_bams/${i}.bam \
-targetIntervals $WORKING_FOLDER/syncfiles/${i}/${i}.hologenome.intervals \
-o $WORKING_FOLDER/syncfiles/${i}/${i}.contaminated_realigned.bam
###rm $output/$sample/${sample}.dedup.bam*

# Run mpileup step 
$samtools mpileup \
$WORKING_FOLDER/syncfiles/${i}/${i}.contaminated_realigned.bam \
-B \
-Q ${base_quality_threshold} \
-f $REFERENCE > $WORKING_FOLDER/syncfiles/${i}/${i}_mpileup.txt

# Transform Mpileup to Sync
python3 $Mpileup2Sync \
--mpileup $WORKING_FOLDER/syncfiles/${i}/${i}_mpileup.txt \
--ref $genome_pickled \
--output $WORKING_FOLDER/syncfiles/${i} \
--base-quality-threshold $base_quality_threshold \
--coding $illumina_quality_coding \
--minIndel $minIndel

# Prepare PoolSNP output
python3 $MaskSYNC_snape_complete \
--sync $WORKING_FOLDER/syncfiles/${i}/${i}.sync.gz \
--output $WORKING_FOLDER/syncfiles/${i}/${i} \
--indel $WORKING_FOLDER/syncfiles/${i}/${i}.indel \
--coverage $WORKING_FOLDER/syncfiles/${i}/${i}.cov \
--mincov $min_cov \
--maxcov $max_cov \
--te $GFF \
--maxsnape $maxsnape

#--------------------------------------------------------------------------------

# Housekeeping

# Rename files
mv $WORKING_FOLDER/syncfiles/${i}/${i}_masked.sync.gz $WORKING_FOLDER/syncfiles/${i}/${i}.masked.sync.gz
# Unzip sync files
gunzip $WORKING_FOLDER/syncfiles/${i}/${i}.masked.sync.gz

$bgzip $WORKING_FOLDER/syncfiles/${i}/${i}.masked.sync
$tabix -s 1 -b 2 -e 2 $WORKING_FOLDER/syncfiles/${i}/${i}.masked.sync

#--------------------------------------------------------------------------------

date
echo "done"