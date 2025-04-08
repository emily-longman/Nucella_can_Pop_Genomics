#### Create guide file for array via an interactive session.

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

# PROTEIN_FILE is the N. canaliculata protein file generated from the gff
PROTEIN_FILE=/netfiles/pespenilab_share/Nucella/processed/N.canaliculata_snpeff_March_2025/protein.fa

#--------------------------------------------------------------------------------

# Currently the protein file is interleaven (i.e., the DNA sequence is in chunks of 80bps a line).
# For easy of running a loop and array, make it not interleaven

# Load modules 
seqtk=/gpfs1/home/e/l/elongman/software/seqtk/seqtk

# Use seqtk to make protein file not interleaven
$seqtk seq $PROTEIN_FILE > $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/Nucella_canaliculata_protein.fa

#--------------------------------------------------------------------------------

# Get list of protein names
 
# Change directory
cd $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology

# Make guide files directory
mkdir guide_files

# Change directory
cd $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/guide_files

# Get all of the scaffold/contig names 
grep ">" $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/Nucella_canaliculata_protein.fa > protein_file_names.txt

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

# Modify in R
module load R/4.4.1  
R

# Load packages
install.packages(c('data.table', 'tidyverse', 'groupdata2'))
library(data.table)
library(tidyverse)
library(groupdata2)

#--------------------------------------------------------------------------------

# Read in the data 
data <- fread("protein_file_names.txt", header = F)
dim(data) # 40634  4

# Break up the protein names up into chunks
group(data, n=45, method = "greedy") -> protein_file_names

# Write guide file
write.table(protein_file_names, "protein_file_names_array.txt", col.names = F, row.names = F, quote = F)
# Note the guide file has dimensions:  40634, 5 (903 partitions with 45 protein names in each)

q()