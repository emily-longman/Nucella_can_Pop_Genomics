#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=DBG2OLC.optimum.2_100_0.01_rerun

# Specify partition
#SBATCH --partition=bigmemwk

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss --7 day limit
#SBATCH --time=03-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=200G

# Request CPUs
#SBATCH --cpus-per-task=32

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# DBG2OLC executable
DBG2OLC=/gpfs1/home/e/l/elongman/software/DBG2OLC

#--------------------------------------------------------------------------------

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/short_read_assembly

# Long reads
#ONT_FOLDER=/netfiles/pespenilab_share/Nucella/raw/ONT
#If you haven't done it yet, unzip the files 
#gunzip $ONT_FOLDER/FC_all.ONT.nuc.fastq.gz
#If you haven't done it yet, unzip the files 
ONT=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/ONT_fltlong/Nuc.2000.fltlong.fastq

# Contigs from Sparse Assembler (using the SparseAssembler with the largest N50 (SparseAssembler_101_2_1))
Contigs=$WORKING_FOLDER/SparseAssembler/SparseAssembler_101_2_1/Contigs.txt

# Extract base files from previous DBG2OLC run
compressed_ONT=$WORKING_FOLDER/DBG2OLC/DBG2OLC_2_20_0.015_rerun/ReadsInfoFrom_Nuc.2000.fltlong.fastq

#--------------------------------------------------------------------------------

# Parameter combinations
# kmerCovTh = 2
# MinOverlap = 100
# AdaptiveTh = 0.01

A=0.01
M=100
K=2

label=KmC_${K}_MinOv_${M}_Adth_${A}

echo ${K} ${M} ${A}

echo ${label}

#--------------------------------------------------------------------------------

cd $WORKING_FOLDER/DBG2OLC

# Make directory for each DBG2OLC parameter combination
if [ -d "DBG2OLC_rerun_${label}" ]
then echo "Working DBG2OLC_rerun_${label} folder exist"; echo "Let's move on."; date
else echo "Working DBG2OLC_rerun_${label} folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/DBG2OLC/DBG2OLC_rerun_${label}; date
fi

#--------------------------------------------------------------------------------

# Use DBG2OLC to construct short but accurate contigs  

# Move to the directory where the output files will be saved
cd $WORKING_FOLDER/DBG2OLC/DBG2OLC_rerun_${label}

# Copy read information from previous run
cp $compressed_ONT ./

# change permissions
chmod 777 *

# Use DBG2OLC to construct short but accurate contigs  
$DBG2OLC \
LD 1 \
k 17 \
AdaptiveTh ${A} \
KmerCovTh ${K} \
MinOverlap ${M} \
RemoveChimera 1 \
Contigs $Contigs \
f $ONT