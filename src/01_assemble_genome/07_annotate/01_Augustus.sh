#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Augustus_softmask

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss
#SBATCH --time=4-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=30G

# Request CPU
#SBATCH --cpus-per-task=5

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out  # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Run Augustus (https://github.com/Gaius-Augustus/Augustus) using apptainer to annotate the genome.

#--------------------------------------------------------------------------------

# Load modules
module load apptainer/1.3.4
AUGUSTUS=/gpfs1/cont/augustus/augustus-3.5.0.sif 

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

# This is the location of the reference genome.
REFERENCE_FULL_ADDRESS=$WORKING_FOLDER/data/processed/genome_assembly/Crassostrea_softmask/N.canaliculata_assembly.fasta.softmasked.fa
REFERENCE=N.canaliculata_assembly.fasta.softmasked.fa

#--------------------------------------------------------------------------------

# Generate Folders and files

# Change directory
cd $WORKING_FOLDER/data/processed/genome_assembly

if [ -d "Augustus" ]
then echo "Working Augustus folder exist"; echo "Let's move on."; date
else echo "Working Augustus folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/genome_assembly/Augustus; date
fi

#--------------------------------------------------------------------------------

# Define Parameters
# Note: AUGUSTUS has currently been trained on species specific training sets to predict genes in the following species.
SPECIES=fly
PROJECT=N.canaliculata

#--------------------------------------------------------------------------------

# Change directory
cd $WORKING_FOLDER/data/processed/genome_assembly/Augustus

# Move copy of Reference to Augustus directory
cp $REFERENCE_FULL_ADDRESS $WORKING_FOLDER/data/processed/genome_assembly/Augustus

# Run Augustus on genome using apptainer
apptainer run \
--home $WORKING_FOLDER/data/processed/genome_assembly/Augustus \
$AUGUSTUS augustus \
--strand=both \
--gff3=on \
--species=${SPECIES} \
${REFERENCE} > \
${PROJECT}.genepred.softmask.gff3
