#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=SparseAssembler

# Specify partition
#SBATCH --partition=bigmemwk

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss --7 day limit
#SBATCH --time=07-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=900G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

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

# If you haven't done it yet, gunzip the files 
#gunzip $WORKING_FOLDER/fastp/*fastq.gz

#--------------------------------------------------------------------------------

# Move to the directory where the output files will be saved
cd $WORKING_FOLDER/SparseAssembler

# Use SparseAssembler to construct short but accurate contigs  
$SparseAssembler \
LD 0 k 51 g 15 \
NodeCovTh 1 \
EdgeCovTh 0 \
GS 2500000000 \
i1 $WORKING_FOLDER/fastp/NC3_R1_clean.fastq \
i2 $WORKING_FOLDER/fastp/NC3_R2_clean.fastq

# In future try k = 31 and k = 71

#--------------------------------------------------------------------------------

# Inform that SparseAssembler is done

echo "done"