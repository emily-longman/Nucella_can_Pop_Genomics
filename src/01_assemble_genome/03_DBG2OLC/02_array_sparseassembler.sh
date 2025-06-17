#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=SparseAssembler_array

# Specify partition
#SBATCH --partition=bigmemwk

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss --7 day limit
#SBATCH --time=07-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=900G

# Submit job array
#SBATCH --array=1-16

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/SparseAssembler_array.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# SparseAssembler executable
SparseAssembler=/gpfs1/home/e/l/elongman/software/SparseAssembler

#--------------------------------------------------------------------------------

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/short_read_assembly
#--------------------------------------------------------------------------------

## Read guide files
# This is a guide file with all of the parameter combinations
# k = 31, 51, 71, 91
# NodeCovTh = 1, 2
# EdgeCovTh = 0, 1

GUIDE_FILE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/short_read_assembly/SparseAssembler/SparseAssembler_GuideFile.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers.
##   k       NodeCovTh     EdgeCovTh   
##   31          1            0
##   31          1            1
##   31          2            0
##   31          2            1
##   51          1            0
##   ....

#--------------------------------------------------------------------------------

# Determine parameter combination to process
k=$( cat $GUIDE_FILE  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $1 }' )
NCT=$( cat $GUIDE_FILE  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $2 }' )
ECT=$( cat $GUIDE_FILE  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $3 }' )

echo $k  $NCT  $ECT

#--------------------------------------------------------------------------------

# If you haven't done it yet, gunzip the files 
#gunzip $WORKING_FOLDER/fastp/*fastq.gz

#--------------------------------------------------------------------------------

# Make directory for each parameter combination
mkdir $WORKING_FOLDER/SparseAssembler/SparseAssembler_${k}_${NCT}_${ECT}

#--------------------------------------------------------------------------------

# Move to the directory where the output files will be saved
cd $WORKING_FOLDER/SparseAssembler/SparseAssembler_${k}_${NCT}_${ECT}

# Use SparseAssembler to construct short but accurate contigs  
$SparseAssembler \
LD 0 k ${k} g 15 \
NodeCovTh ${NCT} \
EdgeCovTh ${ECT} \
GS 2500000000 \
i1 $WORKING_FOLDER/fastp/NC3_R1_clean.fastq \
i2 $WORKING_FOLDER/fastp/NC3_R2_clean.fastq

#--------------------------------------------------------------------------------

# Inform that SparseAssembler is done

echo "done"
