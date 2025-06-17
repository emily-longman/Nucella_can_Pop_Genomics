#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=partition_genome

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss 
#SBATCH --time=10:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=50G

# Submit job array
#SBATCH --array=1-559%20

# Name output of this job using %x=job-name and %j=job-id
#SBATCH -o ./slurmOutput/%x.%A_%a.out 

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will partition the contigs file and backbone into 50 contig chunks. 
# To do so, it will use the guide file produced from the previous R script.
# NOTE: There are 27,922 contigs in the assembly. If the assembly is broken into 50 contig chunks then there are 559 chunks. 

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

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

# Import master partition file
GUIDE_FILE=$WORKING_FOLDER/data/processed/genome_assembly/consensus/dat.win.partitions.txt

echo ${SLURM_ARRAY_TASK_ID}

# Specify the first and last contig for a given partition 
init_bck=$(cat ${GUIDE_FILE} | sed '1d' | awk '{print $1}' | sed "${SLURM_ARRAY_TASK_ID}q;d")
final_bck=$(cat ${GUIDE_FILE} | sed '1d' | awk '{print $2}' | sed "${SLURM_ARRAY_TASK_ID}q;d")
echo "init_bck:" ${init_bck} "final_bck:" ${final_bck}

# Specify backbone name for a given partition
var1=$(echo "^>Backbone_${init_bck}$" ) 
var2=$(echo "^>Backbone_${final_bck}$" ) 
echo "var1:" ${var1} "var2:" ${var2}

#--------------------------------------------------------------------------------

# Generate Folders and files

# Change directory
cd $WORKING_FOLDER/data/processed/genome_assembly/consensus

# This part of the script will check and generate, if necessary, all of the output folders used in the script

# Make genome chunks
if [ -d "gen_chunks" ]
then echo "Working gen_chunks folder exist"; echo "Let's move on."; date
else echo "Working gen_chunks folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/genome_assembly/consensus/gen_chunks; date
fi

# Make info chunks
if [ -d "chunks" ]
then echo "Working chunks folder exist"; echo "Let's move on."; date
else echo "Working chunks folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/genome_assembly/consensus/chunks; date
fi

#--------------------------------------------------------------------------------

# Break the backbone up into partitions
cat $backbone | \
sed -n "/$var1/, /$var2/p" | \
sed '$d'  > $WORKING_FOLDER/data/processed/genome_assembly/consensus/gen_chunks/gen_chunks.$init_bck.$final_bck.fasta

# Break the contig info into partitions 
cat $cons_info | \
sed -n "/$var1/, /$var2/p" | \
sed '$d'  > $WORKING_FOLDER/data/processed/genome_assembly/consensus/chunks/chunk.$init_bck.$final_bck.txt

#--------------------------------------------------------------------------------

# Notify that the process is complete

echo "done"