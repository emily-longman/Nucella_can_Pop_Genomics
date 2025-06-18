#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=cat_contigs_ONT

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1  

# Request CPUs per task
#SBATCH -c 1

# Reserve walltime -- hh:mm:ss --30 hour limit
#SBATCH --time=10:00:00 

# Request memory for the entire job 
#SBATCH --mem=60G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out 

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will concatenate the SparseAssembler contigs and filtered ONT reads.

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

#--------------------------------------------------------------------------------

# Input files for consensus: 
#(1) backbone_raw.fasta by DBG2OLC
#backbone=$WORKING_FOLDER/data/processed/genome_assembly/DBG2OLC/DBG2OLC_KmC_2_MinOv_100_Adth_0.01/DBG2OLC_KmC_2_MinOv_100_Adth_0.01.backbone_raw.fasta
#(2) DBG2OLC_Consensus_info.txt by DBG2OLC
#cons_info=$WORKING_FOLDER/data/processed/genome_assembly/DBG2OLC/DBG2OLC_KmC_2_MinOv_100_Adth_0.01/DBG2OLC_Consensus_info.txt
#(3) DBG contigs (in fasta format)
Contigs=$WORKING_FOLDER/data/processed/genome_assembly/SparseAssembler/SparseAssembler_101_2_1/Contigs.txt
#(4) ONT reads (in fasta format) - converted in prior step
ONT_FA=$WORKING_FOLDER/data/processed/genome_assembly/consensus/Nuc.2000.fltlong.FQtoFA.fasta

#--------------------------------------------------------------------------------

# Change to consensus directory
cd $WORKING_FOLDER/data/processed/genome_assembly/consensus

# Cat contigs and the raw reads for consensus 
cat $Contigs $ONT_FA > ctg_ont.fasta