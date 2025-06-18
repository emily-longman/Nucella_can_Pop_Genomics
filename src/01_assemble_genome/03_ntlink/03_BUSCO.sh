#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=busco_ntlink

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1  

# Request CPUs per task
#SBATCH -c 1

# Reserve walltime -- hh:mm:ss --30 hour limit
#SBATCH --time=30:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=60G

# Submit job array
#SBATCH --array=1-2

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/busco.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# This script runs on BUSCO on the best ntlink assembly. You must cd to that directory before running.

#--------------------------------------------------------------------------------

# Call package (installed with conda)
module load python3.11-anaconda/2024.02-1
source ${ANACONDA_ROOT}/etc/profile.d/conda.sh
#conda create --name busco_env #create and name the environment
conda activate busco_env #activate the environment
#conda install -c conda-forge -c bioconda busco=5.7.1 # install the program

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

#--------------------------------------------------------------------------------

# Run BUSCO in directory where given assembly is.

INPUT=$WORKING_FOLDER/data/processed/genome_assembly/ntlink/ntlink_k.24_w.150/final_assembly.fasta.k24.w150.z1000.ntLink.ntLink.ntLink.ntLink.ntLink.gap_fill.fa.k24.w150.z1000.ntLink.scaffolds.gap_fill.fa

echo $INPUT

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER/data/processed/genome_assembly/ntlink

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

# Make Quast directory 
if [ -d "BUSCO" ]
then echo "Working BUSCO folder exist"; echo "Let's move on."; date
else echo "Working BUSCO folder doesnt exist. Let's fix that."; mkdir cd $WORKING_FOLDER/data/processed/genome_assembly/ntlink/BUSCO; date
fi

#--------------------------------------------------------------------------------

# Change directory
cd $WORKING_FOLDER/data/processed/genome_assembly/ntlink/BUSCO

# Run BUSCO for both the eukaryota and mollusca lineages

if [ ${SLURM_ARRAY_TASK_ID} == 1 ]
then
echo "eval 1"
OUTPUT=./Nucella_consensus_eukaryota
LINEAGE=/netfiles/nunezlab/BUSCO_Lineages/busco_downloads/lineages/eukaryota_odb10
elif [ ${SLURM_ARRAY_TASK_ID} == 2 ]
then
echo "eval 2"
OUTPUT=./Nucella_consensus_mollusca
LINEAGE=/netfiles/nunezlab/BUSCO_Lineages/busco_downloads/lineages/mollusca_odb10
fi

busco -m genome -i $INPUT -o $OUTPUT -l $LINEAGE

#--------------------------------------------------------------------------------

conda deactivate