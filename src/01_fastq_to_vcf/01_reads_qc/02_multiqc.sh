#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=multiqc 

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=2 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=8:00:00

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=20G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu

#--------------------------------------------------------------------------------

# This script will initiate a pipeline which will do some quality QC on the reads.
# This script will use the ouput fastqc files produced in the prior step to generate a multiqc report. 

# Load modules 
# Call package (installed with conda)
module load python3.11-anaconda/2024.02-1
source ${ANACONDA_ROOT}/etc/profile.d/conda.sh
#conda create --name multiqc #If you haven't already done so, create and name the environment
source activate multiqc #activate the environment
#conda install -c bioconda multiqc # If you haven't already done so, install the program
conda activate multiqc 

#--------------------------------------------------------------------------------

# Define important file locations

# RAW_READS indicates the folder where the raw reads are stored.
RAW_READS=/netfiles/pespenilab_share/Nucella/raw/Population_genomics/All_shortreads

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed/fastq_to_bam

# Name of pipeline
PIPELINE=multiqc

#--------------------------------------------------------------------------------

# Generate Folders and Files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

# Change directory
cd $WORKING_FOLDER/qc_reads

# Generating new folder 
if [ -d "multiqc" ]
then echo "Working multiqc folder exist"; echo "Let's move on"; date
else echo "Working multiqc folder doesnt exist. Let's fix that"; mkdir $WORKING_FOLDER/qc_reads/multiqc; date
fi

#--------------------------------------------------------------------------------

# Run multiqc on the fastqc files.
# First run multiqc on all of the reads. Subsequently run multiqc on each lane separately.

# Move to working directory
cd $WORKING_FOLDER/qc_reads

# Run multiqc on all of the reads
multiqc $WORKING_FOLDER/qc_reads/fastqc \
-n multiqc_report_all.html \
-o multiqc

# Run multiqc on each lane individually

# Run multiqc on L006 
multiqc $WORKING_FOLDER/qc_reads/fastqc \
-n multiqc_report_L006.html \
--ignore "*L008*" \
-o multiqc

# Run multiqc on L008
multiqc $WORKING_FOLDER/qc_reads/fastqc \
-n multiqc_report_L008.html \
--ignore "*L006*" \
-o multiqc

#--------------------------------------------------------------------------------

# This part of the pipeline will produce a notification stating the completion of the script. 

echo "pipeline" ${PIPELINE} $(date)

# Deactivate conda
conda deactivate