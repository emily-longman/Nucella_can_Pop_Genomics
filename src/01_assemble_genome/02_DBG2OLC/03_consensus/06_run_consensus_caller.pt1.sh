#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=consensus_pt1

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss 
#SBATCH --time=28:00:00 

# Request memory for the entire job 
#SBATCH --mem=50G

# Submit job array
#SBATCH --array=1-559%15

# Name output of this job using %x=job-name and %j=job-id
#SBATCH -o ./slurmOutput/%x.%A_%a.out 

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run the first step in the consensus script. Specifically, it will run the consensus script "split_and_run_sparc.pt1".

#--------------------------------------------------------------------------------

# Load modules 
module load blasr
export PATH=/gpfs1/sw/rh9/pkgs/stacks/gcc/13.3.0/python/2.7.18/bin:$PATH

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

# Script for split_and_run_sparc.pt1
# NOTE: need to change permissions on supplementary consensus scripts prior to running.
sprun_pt1=$WORKING_FOLDER/src/01_assemble_genome/02_DBG2OLC/03_consensus/consensus_scripts/split_and_run_sparc.pt1.sh

#--------------------------------------------------------------------------------

# Input files for consensus: 
#(1) backbone_raw.fasta by DBG2OLC
#backbone=$WORKING_FOLDER/data/processed/genome_assembly/DBG2OLC/DBG2OLC_KmC_2_MinOv_100_Adth_0.01/DBG2OLC_KmC_2_MinOv_100_Adth_0.01.backbone_raw.fasta
#(2) DBG2OLC_Consensus_info.txt by DBG2OLC
#cons_info=$WORKING_FOLDER/data/processed/genome_assembly/DBG2OLC/DBG2OLC_KmC_2_MinOv_100_Adth_0.01/DBG2OLC_Consensus_info.txt
#(3) DBG contigs (in fasta format)
#Contigs=$WORKING_FOLDER/data/processed/genome_assembly/SparseAssembler/SparseAssembler_101_2_1/Contigs.txt
#(4) ONT reads (in fasta format)
#ONT_FA=$WORKING_FOLDER/data/processed/genome_assembly/consensus/Nuc.2000.fltlong.FQtoFA.fasta

#--------------------------------------------------------------------------------

# Partition guide file (first row is column headers)
GUIDE_FILE=$WORKING_FOLDER/data/processed/genome_assembly/consensus/dat.win.partitions.txt

echo ${SLURM_ARRAY_TASK_ID}

# Specify the first and last contig for a given partition 
init_bck=$(cat ${GUIDE_FILE} | sed '1d' | awk '{print $1}' | sed "${SLURM_ARRAY_TASK_ID}q;d")
final_bck=$(cat ${GUIDE_FILE} | sed '1d' | awk '{print $2}' | sed "${SLURM_ARRAY_TASK_ID}q;d")
echo "init_bck:" ${init_bck} "final_bck:" ${final_bck}
shiftn=`expr $init_bck - 1`
echo ${shiftn}

#--------------------------------------------------------------------------------

# Generate Folders and files

# Change directory
cd $WORKING_FOLDER/data/processed/genome_assembly/consensus

# This part of the script will check and generate, if necessary, all of the output folders used in the script

# Make consensus chunked 
if [ -d "consensus_chunked" ]
then echo "Working consensus_chunked folder exist"; echo "Let's move on."; date
else echo "Working consensus_chunked folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/genome_assembly/consensus/consensus_chunked; date
fi

#--------------------------------------------------------------------------------

# Run consensus

# Change directory
cd $WORKING_FOLDER/data/processed/genome_assembly/consensus

# change permissions
chmod 777 *

# Run consensus (i.e. run split_and_run_sparc.pt1.sh)
$sprun_pt1 \
gen_chunks/gen_chunks.${init_bck}.${final_bck}.fasta \
chunks/chunk.${init_bck}.${final_bck}.txt \
ctg_ont.fasta \
$WORKING_FOLDER/data/processed/genome_assembly/consensus/consensus_chunked \
2 \
32 \
$shiftn > cns_log_pt1.txt 2>&1
# 2>&1 redirects stderr to stdout

echo "done"