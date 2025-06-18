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
# NOTE: prior to running this step the genome was move to netfiles (i.e. long term storage).

#--------------------------------------------------------------------------------

# Load modules
module load apptainer/1.3.4
AUGUSTUS=/gpfs1/cont/augustus/augustus-3.5.0.sif 

#--------------------------------------------------------------------------------

# Define important file locations

# NETFILES_FOLDER is the folder where the genome assembly is located.
NETFILES_FOLDER=/netfiles/pespenilab_share/Nucella/processed/Base_Genome

# This is the location of the reference genome.
REFERENCE_FULL_ADDRESS=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Oct2024/Crassostrea_softmask/N.canaliculata_assembly.fasta.softmasked.fa
REFERENCE=N.canaliculata_assembly.fasta.softmasked.fa

#--------------------------------------------------------------------------------

# Generate Folders and files

# Change directory
cd $NETFILES_FOLDER

if [ -d "annotation" ]
then echo "Working annotation folder exist"; echo "Let's move on."; date
else echo "Working annotation folder doesnt exist. Let's fix that."; mkdir $NETFILES_FOLDER/annotation; date
fi

cd $NETFILES_FOLDER/annotation

if [ -d "Augustus" ]
then echo "Working Augustus folder exist"; echo "Let's move on."; date
else echo "Working Augustus folder doesnt exist. Let's fix that."; mkdir $NETFILES_FOLDER/annotation/Augustus; date
fi

#--------------------------------------------------------------------------------

# Define Parameters
# Note: AUGUSTUS has currently been trained on species specific training sets to predict genes in the following species.
SPECIES=fly
PROJECT=N.canaliculata

#--------------------------------------------------------------------------------

# Change directory
cd $NETFILES_FOLDER/annotation/Augustus

# Move copy of Reference to Augustus directory
cp $REFERENCE_FULL_ADDRESS $NETFILES_FOLDER/annotation/Augustus

# Run Augustus on genome using apptainer
apptainer run \
--home $NETFILES_FOLDER/annotation/Augustus \
$AUGUSTUS augustus \
--strand=both \
--gff3=on \
--species=${SPECIES} \
${REFERENCE} > \
${PROJECT}.genepred.softmask.gff3
