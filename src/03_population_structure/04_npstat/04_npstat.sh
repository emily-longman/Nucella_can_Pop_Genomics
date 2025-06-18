#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=npstat_array

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=20:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=10G 

# Submit job array
#SBATCH --array=1-379%30

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will calculate diversity statistics on each scaffold

# Load modules 
module load gcc/13.3.0-xp3epyt
module load samtools/1.19.2-pfmpoam
module load apptainer/1.3.4
npstat=/gpfs1/home/e/l/elongman/software/npstat/npstat_latest.sif

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed

# BAM_FOLDER is the folder where the cleaned bam files are.
BAM_FOLDER=$WORKING_FOLDER/fastq_to_vcf/RGSM_final_bams

# Guide file with list of populations
POPS=$WORKING_FOLDER/pop_gen/guide_files/Populations.txt

#--------------------------------------------------------------------------------

# Change directory
cd $WORKING_FOLDER/pop_gen

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "npstat" ]
then echo "Working npstat folder exist"; echo "Let's move on."; date
else echo "Working npstat folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/pop_gen/npstat; date
fi

# Change directory
cd $WORKING_FOLDER/pop_gen/npstat

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "Partition_${SLURM_ARRAY_TASK_ID}" ]
then echo "Working Partition_${SLURM_ARRAY_TASK_ID} folder exist"; echo "Let's move on."; date
else echo "Working Partition_${SLURM_ARRAY_TASK_ID} folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/pop_gen/npstat/Partition_${SLURM_ARRAY_TASK_ID}; date
fi

#--------------------------------------------------------------------------------

## Import master partition file 
guide_file=$WORKING_FOLDER/pop_gen/guide_files/npstat_scaffold_guide_file_array.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers. (dimensions: 18919, 2; 379 partitions)
# Scaffold name         # Partition
# Backbone_10001              1
# Backbone_10003              1
# Backbone_10004              1
# Backbone_10005              1
# ....

#--------------------------------------------------------------------------------

# Determine partition to process 

# Change directory
cd $WORKING_FOLDER/pop_gen/npstat/Partition_${SLURM_ARRAY_TASK_ID}

# Echo slurm array task ID
echo ${SLURM_ARRAY_TASK_ID}

# Using the guide file, extract the scaffold names associated based on the Slurm array task ID for a given partition
awk '$2=='${SLURM_ARRAY_TASK_ID}'' $guide_file | awk '{print $1}' > scaffold.names.${SLURM_ARRAY_TASK_ID}.txt

#--------------------------------------------------------------------------------

# For the scaffolds in a given partition, generate a bam file and genome segement for each scaffold then polish each piece. 

# Cat file of pop names and start while loop - loop over populations
cat $POPS | \
while read pop
do echo ${pop}
BAM=$BAM_FOLDER/${pop}.bam

# Cat file of scaffold names and start while loop - loop over scaffolds for a given pool/population
cat scaffold.names.${SLURM_ARRAY_TASK_ID}.txt | \
while read scaffold 
do echo ${scaffold}

# Break up the bam file into each scaffold
samtools view -b ${BAM} ${scaffold} > ${scaffold}.${pop}.bam

# Index the bam file
samtools index ${scaffold}.${pop}.bam

# Generate pileup format
samtools mpileup -r ${scaffold} ${scaffold}.${pop}.bam > ${scaffold}.${pop}.pileup

# Run npstat on a given population and scaffold.
apptainer run $npstat -n 40 -l 25000 -mincov 20 -maxcov 120 -minqual 30 -nolowfreq 5 ${scaffold}.${pop}.pileup

# n sample size: haploid sample size
# mincov minimum_coverage: filter on minimum coverage
# maxcov maximum_coverage: filter on maximum coverage
# minqual minimum_base_quality: filter on base quality
# nolowfeq m: filter on minimum allele count > m

# Add file name to stats text file
for f in ${scaffold}.${pop}.pileup.stats ; do sed -i "s/$/\t$f/" $f; done

# Housekeeping - remove intermediate files
rm ${scaffold}.${pop}.bam
rm ${scaffold}.${pop}.bam.bai
rm ${scaffold}.${pop}.pileup

# End scaffold loop
done

# Merge stats for each scaffold together
sort -u *${pop}.pileup.stats > npstat.${pop}.${SLURM_ARRAY_TASK_ID}.txt

# Remove header 
grep -v "window" npstat.${pop}.${SLURM_ARRAY_TASK_ID}.txt > npstat.${pop}.${SLURM_ARRAY_TASK_ID}.clean.txt

# Housekeeping - remove intermediate files
rm *.pileup.stats
rm npstat.${pop}.${SLURM_ARRAY_TASK_ID}.txt

# End population loop
done 

#--------------------------------------------------------------------------------

# Housekeeping - remove intermediate files
rm scaffold.names.${SLURM_ARRAY_TASK_ID}.txt
