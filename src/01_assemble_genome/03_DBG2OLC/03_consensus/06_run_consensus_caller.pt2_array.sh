#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=consensus_pt2_array

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss 
#SBATCH --time=5:00:00 

# Request memory for the entire job 
#SBATCH --mem=40G

# Request CPUs
#SBATCH --cpus-per-task=32

# Submit job array
#SBATCH --array=1-931%15

# Name output of this job using %x=job-name and %j=job-id
#SBATCH -o ./slurmOutput/%x.%A_%a.out 

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run the second step in the consensus script. Specifically it will run the consensus script "split_and_run_sparc.pt2_array.sh".

#--------------------------------------------------------------------------------

# Load modules
module load blasr
export PATH=/gpfs1/sw/rh9/pkgs/stacks/gcc/13.3.0/python/2.7.18/bin:$PATH

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

# Script for split_and_run_sparc.pt2_array.sh
# NOTE: need to change permissions on supplementary consensus scripts prior to running.
sprun_pt2=$WORKING_FOLDER/src/01_assemble_genome/03_DBG2OLC/03_consensus/consensus_scripts/split_and_run_sparc.pt2_array.sh

#--------------------------------------------------------------------------------

# Input files for consensus: 
#(1) backbone_raw.fasta by DBG2OLC
backbone=$WORKING_FOLDER/data/processed/genome_assembly/DBG2OLC/DBG2OLC_KmC_2_MinOv_100_Adth_0.01/DBG2OLC_KmC_2_MinOv_100_Adth_0.01.backbone_raw.fasta
#(2) DBG2OLC_Consensus_info.txt by DBG2OLC
cons_info=$WORKING_FOLDER/data/processed/genome_assembly/DBG2OLC/DBG2OLC_KmC_2_MinOv_100_Adth_0.01/DBG2OLC_Consensus_info.txt
#(3) DBG contigs (in fasta format)
#Contigs=$WORKING_FOLDER/data/processed/genome_assembly/SparseAssembler/SparseAssembler_101_2_1/Contigs.txt
#(4) ONT reads (in fasta format)
#ONT_FA=$WORKING_FOLDER/data/processed/genome_assembly/consensus/Nuc.2000.fltlong.FQtoFA.fasta

#--------------------------------------------------------------------------------

# Import master partition file (created in step 05_guide_file_pt2.R)
GUIDE_FILE=$WORKING_FOLDER/data/processed/genome_assembly/consensus/guide_file_array.txt

echo ${SLURM_ARRAY_TASK_ID}

# Move to working directory
cd $WORKING_FOLDER/data/processed/genome_assembly/consensus

# Using the guide file, extract the rows associated based on the Slurm array task ID
awk '$2=='${SLURM_ARRAY_TASK_ID}'' $GUIDE_FILE | awk '{print $1}' > backbone.names.${SLURM_ARRAY_TASK_ID}.txt

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER/data/processed/genome_assembly/consensus

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "final_assembly_array" ]
then echo "Working final_assembly_array folder exist"; echo "Let's move on."; date
else echo "Working final_assembly_array folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/genome_assembly/consensus/final_assembly_array; date
fi

if [ -d "Logs" ]
then echo "Working Logs folder exist"; echo "Let's move on."; date
else echo "Working Logs folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/genome_assembly/consensus/Logs; date
fi

#--------------------------------------------------------------------------------

# We need to open a lot of files to distribute the above file into lots of smaller files

# change permissions for consensus
chmod 777 *

# change permissions for accompanying scripts
cd $WORKING_FOLDER/src/01_assemble_genome/03_DBG2OLC/03_consensus/consensus_scripts
chmod 777 *

#--------------------------------------------------------------------------------

# Run consensus

# Move to working directory
cd $WORKING_FOLDER/data/processed/genome_assembly/consensus

# Run consensus (i.e. run split_and_run_sparc.pt2_array.sh)
$sprun_pt2 \
$backbone \
$cons_info \
ctg_ont.fasta \
$WORKING_FOLDER/data/processed/genome_assembly/consensus/consensus_chunked \
2 \
32 > Logs/cns_log_pt2_array.${SLURM_ARRAY_TASK_ID}.txt 2>&1
# 2>&1 redirects stderr to stdout

#--------------------------------------------------------------------------------

# Remove backbone names file
rm backbone.names.${SLURM_ARRAY_TASK_ID}.txt

echo "done"