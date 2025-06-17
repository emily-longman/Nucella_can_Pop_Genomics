#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=partition_genome

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss 
#SBATCH --time=10:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=50G

# Submit job array
#SBATCH --array=1-560%20

# Name output of this job using %x=job-name and %j=job-id
#SBATCH -o ./slurmOutput/partition_genome.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# This script will partition the genome into chunks. 

# It will use the txt file produced from the  previous R script which calculated the number of chunks necessary to split the assembly into 50 contig chunks

# There are 27922 contigs in the assembly. If the assembly is broken into 50 contig chunks then there are 559 chunks. 

#--------------------------------------------------------------------------------

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly
#WORKING_FOLDER_NETFILES=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/short_read_assembly

#--------------------------------------------------------------------------------

# Input files for consensus: 
#(1) backbone_raw.fasta by DBG2OLC
backbone=$WORKING_FOLDER_SCRATCH/DBG2OLC/DBG2OLC_KmC_2_MinOv_100_Adth_0.01/DBG2OLC_KmC_2_MinOv_100_Adth_0.01.backbone_raw.fasta
#(2) DBG2OLC_Consensus_info.txt by DBG2OLC
cons_info=$WORKING_FOLDER_SCRATCH/DBG2OLC/DBG2OLC_KmC_2_MinOv_100_Adth_0.01/DBG2OLC_Consensus_info.txt
#(3) DBG contigs (in fasta format)
#Contigs=$WORKING_FOLDER_NETFILES/SparseAssembler/SparseAssembler_101_2_1/Contigs.txt
#(4) ONT reads (in fasta format)
#ONT_FA=$WORKING_FOLDER_SCRATCH/consensus/Nuc.2000.fltlong.FQtoFA.fasta

#--------------------------------------------------------------------------------

# Import master partition file
guide=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly/consensus/dat.win.partitions.txt

echo ${SLURM_ARRAY_TASK_ID}

init_bck=$(cat ${guide} | sed '1d' | awk '{print $1}' | sed "${SLURM_ARRAY_TASK_ID}q;d")
final_bck=$(cat ${guide} | sed '1d' | awk '{print $2}' | sed "${SLURM_ARRAY_TASK_ID}q;d")
echo $init_bck $final_bck

var1=$(echo "^>Backbone_${init_bck}$" ) 
var2=$(echo "^>Backbone_${final_bck}$" ) 
echo $var1 $var2

#--------------------------------------------------------------------------------

# Move to consensus directory
cd $WORKING_FOLDER_SCRATCH/consensus/

# Make genome chunks
if [ -d "gen_chunks" ]
then echo "Working gen_chunks folder exist"; echo "Let's move on."; date
else echo "Working gen_chunks folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/consensus/gen_chunks; date
fi
cat $backbone | \
sed -n "/$var1/, /$var2/p" | \
sed '$d'  > $WORKING_FOLDER_SCRATCH/consensus/gen_chunks/gen_chunks.$init_bck.$final_bck.fasta

# Make info chunks
if [ -d "chunks" ]
then echo "Working chunks folder exist"; echo "Let's move on."; date
else echo "Working chunks folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/consensus/chunks; date
fi
cat $cons_info | \
sed -n "/$var1/, /$var2/p" | \
sed '$d'  > $WORKING_FOLDER_SCRATCH/consensus/chunks/chunk.$init_bck.$final_bck.txt

#--------------------------------------------------------------------------------

echo "done"