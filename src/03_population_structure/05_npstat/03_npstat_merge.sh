#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=npstat_merge

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=5:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=10G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will merge the npstat diversity statistics together.

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

# Guide file with list of populations
POPS=$WORKING_FOLDER/guide_files/Populations.txt

# Number of partitions
Partitions=379

#--------------------------------------------------------------------------------

# For each population, merge the npstat files from each partition. 


# Change directory
cd $WORKING_FOLDER/data/processed/pop_structure/npstat

# Cat file of pop names and start while loop - loop over populations
cat $POPS | \
while read pop
do 
echo ${pop}

for i in $(seq 1 $Partitions); 
do 
echo ${i}

# Merge stats for each scaffold together
cat $WORKING_FOLDER/data/processed/pop_structure/npstat/Partition_${i}/npstat.${pop}.${i}.clean.txt >> $WORKING_FOLDER/data/processed/pop_structure/npstat/npstat.${pop}.txt

# End loop over partitions
done 

# End population loop
done 

