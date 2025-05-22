#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=fastqc_cDNA

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 # on one node
#SBATCH --ntasks-per-node=5  

# Reserve walltime -- hh:mm:ss 
#SBATCH --time=7-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=40G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# Call packages
spack load fastqc@0.11.7 # Call fastqc package 

#--------------------------------------------------------------------------------

#Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly

#RAW cDNA indicates the folder where the raw reads are stored.
RAW_READS=/netfiles/pespenilab_share/Nucella/raw/cDNA/Nucella.ONT.cDNA.barcode12.fastq.gz

#--------------------------------------------------------------------------------

cd $WORKING_FOLDER_SCRATCH

# Make fastqc directory 
if [ -d "cDNA_fastqc_barcode12" ]
then echo "Working cDNA_fastqc_barcode12 folder exist"; echo "Let's move on."; date
else echo "Working cDNA_fastqc_barcode12 folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/cDNA_fastqc_barcode12; date
fi

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER_SCRATCH

# Lets do some QC on the reads
fastqc $RAW_READS \
--outdir $WORKING_FOLDER_SCRATCH/cDNA_fastqc_barcode12