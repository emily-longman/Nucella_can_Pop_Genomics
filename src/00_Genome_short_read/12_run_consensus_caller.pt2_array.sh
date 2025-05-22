#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=consensus_pt2_array

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss 
#SBATCH --time=5:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=40G

# Request CPUs
#SBATCH --cpus-per-task=32

# Submit job array
#SBATCH --array=1-931%15

# Name output of this job using %x=job-name and %j=job-id
#SBATCH -o ./slurmOutput/consensus_pt2.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# This script will run the second half of the consensus script. 

#--------------------------------------------------------------------------------

# Load modules
module load blasr
spack load python@2.7.18

#--------------------------------------------------------------------------------

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly
WORKING_FOLDER_NETFILES=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/short_read_assembly

#--------------------------------------------------------------------------------
# Change to consensus directory
#cd $WORKING_FOLDER_SCRATCH/consensus

# Input files for consensus: 
#(1) backbone_raw.fasta by DBG2OLC
backbone=$WORKING_FOLDER_SCRATCH/DBG2OLC/DBG2OLC_KmC_2_MinOv_100_Adth_0.01/DBG2OLC_KmC_2_MinOv_100_Adth_0.01.backbone_raw.fasta
#(2) DBG2OLC_Consensus_info.txt by DBG2OLC
cons_info=$WORKING_FOLDER_SCRATCH/DBG2OLC/DBG2OLC_KmC_2_MinOv_100_Adth_0.01/DBG2OLC_Consensus_info.txt
#(3) DBG contigs (in fasta format)
#Contigs=$WORKING_FOLDER_NETFILES/SparseAssembler/SparseAssembler_101_2_1/Contigs.txt
#(4) ONT reads (in fasta format) - converted in step 9 part 1
#ONT_FA=$WORKING_FOLDER_SCRATCH/consensus/Nuc.2000.fltlong.FQtoFA.fasta

# Cat contigs and the raw reads for consensus - cat in step 9 part 2
#cat $Contigs $ONT_FA > ctg_ont.fasta

#--------------------------------------------------------------------------------

## Import master partition file 
guide_file=$WORKING_FOLDER_SCRATCH/consensus/guide_file_array.txt

echo ${SLURM_ARRAY_TASK_ID}

# Move to working directory
cd $WORKING_FOLDER_SCRATCH/consensus

# Using the guide file, extract the rows associated based on the Slurm array task ID
awk '$2=='${SLURM_ARRAY_TASK_ID}'' $guide_file | awk '{print $1}' > backbone.names.${SLURM_ARRAY_TASK_ID}.txt

#--------------------------------------------------------------------------------
# Move to working directory
cd $WORKING_FOLDER_SCRATCH/consensus

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "final_assembly_array" ]
then echo "Working final_assembly_array folder exist"; echo "Let's move on."; date
else echo "Working final_assembly_array folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/consensus/final_assembly_array; date
fi

if [ -d "Logs" ]
then echo "Working Logs folder exist"; echo "Let's move on."; date
else echo "Working Logs folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/consensus/Logs; date
fi

#--------------------------------------------------------------------------------

# We need to open a lot of files to distribute the above file into lots of smaller files

# change permissions for consensus
chmod 777 *

# change permissions for accompanying scripts
cd /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/00_Genome_short_read/12_consensus_scripts_extra
chmod 777 *

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER_SCRATCH/consensus

# Run consensus

### Run consensus (i.e. run split_and_run_sparc.pt2.array.sh)
sprun_pt2=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/00_Genome_short_read/12_consensus_scripts_extra/split_and_run_sparc.pt2_array.sh

$sprun_pt2 \
$backbone \
$cons_info \
ctg_ont.fasta \
$WORKING_FOLDER_SCRATCH/consensus/consensus_dir_chunked \
2 \
32 > Logs/cns_log_pt2_array.${SLURM_ARRAY_TASK_ID}.txt 2>&1
# 2>&1 redirects stderr to stdout

#--------------------------------------------------------------------------------

# Remove backbone names file
rm backbone.names.${SLURM_ARRAY_TASK_ID}.txt

echo "done"