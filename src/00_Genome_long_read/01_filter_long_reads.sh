#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=fltlong 

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --cpus-per-task=40 
#SBATCH --nodes=1 # on one node

# Reserve walltime -- hh:mm:ss
#SBATCH --time=24:00:00 #<= this may depend on your resources

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=80G #<= this may depend on your resources

# Submit job array
#  #SBATCH --array=1-4

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/fltlong.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to


filtlong=/gpfs1/home/e/l/elongman/software/Filtlong/bin/filtlong
input=/netfiles/pespenilab_share/Nucella/raw/ONT/FC_all.ONT.nuc.fastq.gz
#input=./FC_all.ONT.nuc.fastq.gz

#arr=(1000 2000 3500 5000)
#L="${arr[$SLURM_ARRAY_TASK_ID]}"
#echo $L

$filtlong --min_length 1000 $input | gzip > Nuc.1000.fltlong.fastq.gz
$filtlong --min_length 2000 $input | gzip > Nuc.2000.fltlong.fastq.gz
$filtlong --min_length 3500 $input | gzip > Nuc.3500.fltlong.fastq.gz
$filtlong --min_length 5000 $input | gzip > Nuc.5000.fltlong.fastq.gz
$filtlong --min_length 7500 $input | gzip > Nuc.5000.fltlong.fastq.gz
$filtlong --min_length 10000 $input | gzip > Nuc.5000.fltlong.fastq.gz

# Code from Filtlong github
# filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 input.fastq.gz | gzip > output.fastq.gz

# $filtlong \
# --min_length $L \
# $input | gzip > Nuc.$L.fltlong.fastq.gz