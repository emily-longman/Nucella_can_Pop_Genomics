#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=SNPeff

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --7 day limit 
#SBATCH --time=20:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=10G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# SNPeff code

# Create directory 
#cd /netfiles/pespenilab_share/Nucella/processed
#mkdir N.can_genome_Dec2024

# Re-save genome as sequences.fa then gzip
#awk '/^>/ {$0=$1} 1' N.canaliculata_assembly.fasta.softmasked > sequences.fa
#gzip sequences.fa

# Move gtf file to directory and rename
#scp braker.gtf /netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024/
#mv braker.gtf genes.gtf

# Convert gtf to gff. Get script from: https://github.com/nextgenusfs/augustus/blob/master/scripts/gtf2gff.pl
#perl /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/00_Genome_short_read/24_gtf_to_gff.pl \
#< /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly/braker/braker_cDNA/braker/braker.gtf \
#-o myfile.gff

## Re-do conversion in Jan 2025
##perl /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/00_Genome_short_read/24_gtf_to_gff.pl \
##< /netfiles/pespenilab_share/Nucella/processed/Base_Genome/annotation/braker/braker_cDNA/braker/braker.gtf \
##-o N.can.gff

# Move gff and rename
#scp myfile.gff /netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024/
#mv myfile.gff N.can.gff

# Create protein file using AGAT
#module load singularity
#cd /netfiles/nunezlab/Shared_Resources/Software/AGAT
#singularity run agat_1.0.0--pl5321hdfd78af_0.sif
#agat_sp_extract_sequences.pl --gff /netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024/N.can.gff \
#-f /netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024/N.canaliculata_assembly.fasta.softmasked -p -o /netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024/protein.fa

# Move SNPeff config file to directory and add Nucella directory (N.can_genome_Dec2024) to end of script

#--------------------------------------------------------------------------------

# Create directory 
#cd /netfiles/pespenilab_share/Nucella/processed
#mkdir N.can_genome_Dec2024

# Re-save genome as name of directory .fa and as sequences.fa
#scp /netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Oct2024/Crassostrea_softmask/N.canaliculata_assembly.fasta.softmasked /netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024
#cd /netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024
#mv N.canaliculata_assembly.fasta.softmasked N.can_genome_Dec2024.fa
#chmod 777 N.can_genome_Dec2024.fa
#scp /netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Oct2024/Crassostrea_softmask/N.canaliculata_assembly.fasta.softmasked /netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024
#mv N.canaliculata_assembly.fasta.softmasked sequences.fa
#chmod 777 sequences.fa

# Move gtf file to directory and rename to "genes.gtf"
#scp braker.gtf /netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024/
#mv braker.gtf genes.gtf
#chmod 777 genes.gtf

# Build directory

# Load SNPeff
module load snpeff

# Define important file locations
datdir=/netfiles/pespenilab_share/Nucella/processed
param=/netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024/snpEff.config

snpeff build -dataDir $datdir -c $param  -gtf22 -noCheckCds -v N.can_genome_Dec2024
