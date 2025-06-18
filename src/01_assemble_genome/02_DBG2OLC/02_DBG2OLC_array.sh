#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=DBG2OLC_array

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 #
#SBATCH --ntasks-per-node=1  

# Request CPUs per task
#SBATCH --cpus-per-task=1

# Reserve walltime -- hh:mm:ss --30 hour limit
#SBATCH --time=07-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=900G

# Submit job array
#SBATCH --array=1-12

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out 

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will use the best assembly from SparseAssembler as well as the 2000 filt ONT data as the inputs to DBG20LC (https://github.com/yechengxi/DBG2OLC). 
# It will then assess the assemblies with quast.
# NOTE: This script will run in an array structure and will produce multiple assemblies based on the parameters specified in the guide file.

# Load modules  
DBG2OLC=/gpfs1/home/e/l/elongman/software/DBG2OLC
quast=/netfiles/nunezlab/Shared_Resources/Software/quast-5.2.0/quast.py

#--------------------------------------------------------------------------------

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

#--------------------------------------------------------------------------------

# Input files 

# Filtered long reads (note: must unzip the filt length prior to using)
ONT=$WORKING_FOLDER/data/processed/genome_assembly/ONT_fltlong/Nuc.2000.fltlong.fastq

# Contigs from Sparse Assembler (using the SparseAssembler with the largest N50 (SparseAssembler_101_2_1))
Contigs=$WORKING_FOLDER/data/processed/genome_assembly/SparseAssembler/SparseAssembler_101_2_1/Contigs.txt

#--------------------------------------------------------------------------------

# Read guide files (split up into two guide files becuase so many assembly outputs) 
GUIDE_FILE=$WORKING_FOLDER/guide_files/DBG2OLC_GuideFile.txt

# This is a guide file with all of the parameter combinations
# kmerCovTh = 2, 4, 6, 8, 10
# MinOverlap = 30, 50, 100, 150
# AdaptiveTh = 0.001, 0.01

#Example: -- the headers are just for descriptive purposes. The actual file has no headers.
##   kmerCovTh   MinOverlap       AdaptiveTh   
##   2             30               0.001
##   2             50               0.001
##   2             100              0.001
##   2             150              0.001
##   2             30                0.01
##   2             50                0.01
##   ...

#--------------------------------------------------------------------------------

# Determine parameter combination to process
K=$( cat $GUIDE_FILE  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $1 }' )
M=$( cat $GUIDE_FILE  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $2 }' )
A=$( cat $GUIDE_FILE  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $3 }' )

label=DBG2OLC_KmC_${K}_MinOv_${M}_Adth_${A}

echo "label:" ${label} 
echo "K:" ${K} "M:" ${M} "A:" ${A}

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/genome_assembly

# This part of the script will check and generate, if necessary, all of the output folders used in the script

# Make SparseAssembler directory
if [ -d "DBG2OLC" ]
then echo "Working DBG2OLC folder exist"; echo "Let's move on."; date
else echo "Working DBG2OLC folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/genome_assembly/DBG2OLC; date
fi

# Change directory
cd $WORKING_FOLDER/data/processed/genome_assembly/DBG2OLC

# Make directory for each DBG2OLC parameter combination

if [ -d "${label}" ]
then echo "Working ${label} folder exist"; echo "Let's move on."; date
else echo "Working ${label} folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/genome_assembly/DBG2OLC/${label}; date
fi

#--------------------------------------------------------------------------------

# Generate folders for Quast output

# Change directory
cd $WORKING_FOLDER/data/processed/genome_assembly/Quast

# Make sparse_assembler directory 
if [ -d "DBG2OLC" ]
then echo "Working DBG2OLC folder exist"; echo "Let's move on."; date
else echo "Working DBG2OLC folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/genome_assembly/Quast/DBG2OLC; date
fi

# Change directory
cd $WORKING_FOLDER/data/processed/genome_assembly/Quast/DBG2OLC

# Make Quast directory for each parameter combination
if [ -d "${label}" ]
then echo "Working ${label} folder exist"; echo "Let's move on."; date
else echo "Working ${label} folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/genome_assembly/Quast/DBG2OLC/${label}; date
fi

#--------------------------------------------------------------------------------

# Use DBG2OLC to construct assembly.  

# Move to the directory where the output files will be saved
cd $WORKING_FOLDER/data/processed/genome_assembly/DBG2OLC/${label}

# Run DBG2OLC
$DBG2OLC \
k 17 \
KmerCovTh ${K} \
AdaptiveTh ${A} \
MinOverlap ${M} \
RemoveChimera 1 \
Contigs $Contigs \
f $ONT

#--------------------------------------------------------------------------------

# Rename output backbone file with label name

echo "renaming outputs"

mv backbone_raw.fasta ${label}.backbone_raw.fasta

#--------------------------------------------------------------------------------

# Run quast
$quast $WORKING_FOLDER/data/processed/genome_assembly/DBG2OLC/${label}/${label}.backbone_raw.fasta \
-o $WORKING_FOLDER/data/processed/genome_assembly/Quast/DBG2OLC/${label}

#--------------------------------------------------------------------------------

# Inform that DBG2OLC is done

echo "done"