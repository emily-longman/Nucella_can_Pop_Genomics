#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Augustus_hardmask

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

# Run Augustus on apptainer.

#Load modules
module load apptainer/1.3.4
AUGUSTUS=/gpfs1/cont/augustus/augustus-3.5.0.sif 

#--------------------------------------------------------------------------------

# Define important file locations

# HOME is core folder where this pipeline is being run.
HOME=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/annotation

# This is the location of the reference genome.
REFERENCE_FULL_ADDRESS=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Oct2024/Crassostrea_mask/N.canaliculata_assembly.fasta.masked
REFERENCE=N.canaliculata_assembly.fasta.masked

#--------------------------------------------------------------------------------

# Generate Folders and files

# Change directory
cd $HOME

if [ -d "Augustus" ]
then echo "Working Augustus folder exist"; echo "Let's move on."; date
else echo "Working Augustus folder doesnt exist. Let's fix that."; mkdir $HOME/Augustus; date
fi

#--------------------------------------------------------------------------------

# Define Parameters
# Note: AUGUSTUS has currently been trained on species specific training sets to predict genes in the following species.
SPECIES=fly
PROJECT=N.canaliculata

#--------------------------------------------------------------------------------

# Change directory
cd $HOME/Augustus

# Move copy of Reference to Augustus directory
cp $REFERENCE_FULL_ADDRESS $HOME/Augustus

# Run Augustus on genome using apptainer
apptainer run \
--home $HOME/Augustus \
$AUGUSTUS augustus \
--strand=both \
--gff3=on \
--species=${SPECIES} \
${REFERENCE} > \
${PROJECT}.genepred.hardmask.gff3
