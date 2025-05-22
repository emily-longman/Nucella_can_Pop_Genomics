#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=SNPeff

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --7 day limit 
#SBATCH --time=30:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=20G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Generate SnpEff database

# Load software  
module load snpeff/5.2c

#--------------------------------------------------------------------------------

# Define important file locations

NETFILES_DIR=/netfiles/pespenilab_share/Nucella/processed
REFERENCE=$NETFILES_DIR/Base_Genome/Base_Genome_Oct2024/Crassostrea_softmask/N.canaliculata_assembly.fasta.softmasked

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $NETFILES_DIR

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "N.canaliculata_snpeff_March_2025" ]
then echo "Working N.canaliculata_snpeff_March_2025 folder exist"; echo "Let's move on."; date
else echo "Working N.canaliculata_snpeff_March_2025 folder doesnt exist. Let's fix that."; mkdir $NETFILES_DIR/N.canaliculata_snpeff_March_2025; date
fi

#--------------------------------------------------------------------------------
# Move a copy of the genome, gff, and config file to directory

# Change directory 
#cd $NETFILES_DIR/N.canaliculata_snpeff_March_2025

# Move a copy of the genome and rename as sequences.fa, then move another copy of the genome
#scp $REFERENCE $NETFILES_DIR/N.canaliculata_snpeff_March_2025
#mv N.canaliculata_assembly.fasta.softmasked sequences.fa
#scp $REFERENCE $NETFILES_DIR/N.canaliculata_snpeff_March_2025

# Move a copy of the gff and rename as genes.gff
#scp $NETFILES_DIR/Base_Genome/annotation/Augustus/N.canaliculata.genepred.softmask.gff3 $NETFILES_DIR/N.canaliculata_snpeff_March_2025
#mv N.canaliculata.genepred.softmask.gff3 genes.gff

# Move a copy of the config file 
#scp $NETFILES_DIR/N.can_genome_Dec2024/snpEff.config $NETFILES_DIR/N.canaliculata_snpeff_March_2025

# Note: edit config file to add N.canaliculata_snpeff_March_2025 directory at end: 
### N.canaliculata_snpeff_March_2025
## N.canaliculata_snpeff_March_2025.genome : N.canaliculata_snpeff_March_2025

#--------------------------------------------------------------------------------

# Create protein file using AGAT (used RHEL7)

# Load AGAT using singularity
#module load singularity
#AGAT=/netfiles/nunezlab/Shared_Resources/Software/AGAT/agat_1.0.0--pl5321hdfd78af_0.sif

# Run AGAT using singularity to make a protein file
#singularity run $AGAT

#agat_sp_extract_sequences.pl --gff /netfiles/pespenilab_share/Nucella/processed/N.canaliculata_snpeff_March_2025/genes.gff \
#-f /netfiles/pespenilab_share/Nucella/processed/N.canaliculata_snpeff_March_2025/sequences.fa \
#-p -o /netfiles/pespenilab_share/Nucella/processed/N.canaliculata_snpeff_March_2025/protein.fa

#--------------------------------------------------------------------------------

# Change directory
cd $NETFILES_DIR/N.canaliculata_snpeff_March_2025

# Define important file locations
DATA_DIR=$NETFILES_DIR
PARAM=$NETFILES_DIR/N.canaliculata_snpeff_March_2025/snpEff.config

snpeff build -dataDir $DATA_DIR -c $PARAM -gff3 -noCheckCds -v N.canaliculata_snpeff_March_2025
