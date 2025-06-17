#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=SparseAssembler_array

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss 
#SBATCH --time=07-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=900G

# Submit job array
#SBATCH --array=1-28

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out 

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will use the AVITI short reads to assemble contigs with SparseAssembler (https://github.com/yechengxi/SparseAssembler).
# Then it will assess the assemblies with quast.
# NOTE: This script will run in an array structure and will produce multiple assemblies based on the parameters specified in the guide file.

# Load modules  
SparseAssembler=/gpfs1/home/e/l/elongman/software/SparseAssembler
quast=/netfiles/nunezlab/Shared_Resources/Software/quast-5.2.0/quast.py

#--------------------------------------------------------------------------------

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

#--------------------------------------------------------------------------------

# Read guide file
GUIDE_FILE=$WORKING_FOLDER/guide_files/SparseAssembler_GuideFile.txt

# Parameters:
# k = 31, 51, 71, 91, 101, 111, 127
# NodeCovTh = 1, 2
# EdgeCovTh = 0, 1

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

echo "k:" ${k} "NCT:" ${NCT} "ECT:" ${ECT}

#--------------------------------------------------------------------------------

# If you haven't done it yet, gunzip the files 
gunzip $WORKING_FOLDER/data/processed/genome_assembly/trim_AVITI/*fastq.gz

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/genome_assembly

# This part of the script will check and generate, if necessary, all of the output folders used in the script

# Make SparseAssembler directory
if [ -d "SparseAssembler" ]
then echo "Working SparseAssembler folder exist"; echo "Let's move on."; date
else echo "Working SparseAssembler folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/genome_assembly/SparseAssembler; date
fi

# Change directory
cd $WORKING_FOLDER/data/processed/genome_assembly/SparseAssembler

# Make SparseAssembler directory for parameter combination
if [ -d "SparseAssembler_${k}_${NCT}_${ECT}" ]
then echo "Working SparseAssembler_${k}_${NCT}_${ECT} folder exist"; echo "Let's move on."; date
else echo "Working SparseAssembler_${k}_${NCT}_${ECT} folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/genome_assembly/SparseAssembler/SparseAssembler_${k}_${NCT}_${ECT}; date
fi

#--------------------------------------------------------------------------------

# Generate folders for Quast output

# Move to working directory
cd $WORKING_FOLDER/data/processed/genome_assembly

# Make Quast directory 
if [ -d "Quast" ]
then echo "Working Quast folder exist"; echo "Let's move on."; date
else echo "Working Quast folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/genome_assembly/Quast; date
fi

# Change directory
cd $WORKING_FOLDER/data/processed/genome_assembly/Quast

# Make sparse_assembler directory 
if [ -d "SparseAssembler" ]
then echo "Working SparseAssembler folder exist"; echo "Let's move on."; date
else echo "Working SparseAssembler folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/genome_assembly/Quast/SparseAssembler; date
fi

# Change directory
cd $WORKING_FOLDER/data/processed/genome_assembly/Quast/SparseAssembler

# Make sparse_assembler directory 
if [ -d "SparseAssembler_${k}_${NCT}_${ECT}" ]
then echo "Working SparseAssembler_${k}_${NCT}_${ECT} folder exist"; echo "Let's move on."; date
else echo "Working SparseAssembler_${k}_${NCT}_${ECT} folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/genome_assembly/Quast/SparseAssembler/SparseAssembler_${k}_${NCT}_${ECT}; date
fi

#--------------------------------------------------------------------------------

# Move to the directory where the output files will be saved
cd $WORKING_FOLDER/data/processed/genome_assembly/SparseAssembler/SparseAssembler_${k}_${NCT}_${ECT}

# Use SparseAssembler to construct short but accurate contigs  
$SparseAssembler \
LD 0 k ${k} g 15 \
NodeCovTh ${NCT} \
EdgeCovTh ${ECT} \
GS 2500000000 \
i1 $WORKING_FOLDER/data/processed/genome_assembly/trim_AVITI/NC3_R1_clean.fastq \
i2 $WORKING_FOLDER/data/processed/genome_assembly/trim_AVITI/NC3_R2_clean.fastq

#--------------------------------------------------------------------------------

# Run quast
$quast \
$WORKING_FOLDER/data/processed/genome_assembly/SparseAssembler/SparseAssembler_${k}_${NCT}_${ECT}/Contigs.txt \
-o $WORKING_FOLDER/data/processed/genome_assembly/Quast/SparseAssembler/SparseAssembler_${k}_${NCT}_${ECT}

#--------------------------------------------------------------------------------

# Inform that SparseAssembler is done

echo "done"
