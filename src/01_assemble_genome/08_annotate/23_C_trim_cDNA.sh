#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=trim_cDNA

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

# Trim and filter ONT cDNA using nanoq (https://github.com/esteinig/nanoq)

#--------------------------------------------------------------------------------

# Call package (installed with conda)
module load python3.11-anaconda/2023.09-0
source ${ANACONDA_ROOT}/etc/profile.d/conda.sh
#conda create --name nanoq #If you haven't already done so, create and name the environment
conda activate nanoq #activate the environment
#conda install -c conda-forge -c bioconda nanoq # If you haven't already done so, install the program
 
#--------------------------------------------------------------------------------

#Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly

#RAW cDNA indicates the folder where the raw reads are stored.
RAW_READS=/netfiles/pespenilab_share/Nucella/raw/cDNA/Nucella.ONT.cDNA.barcode12.fastq.gz

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER_SCRATCH

if [ -d "cDNA_trim" ]
then echo "Working cDNA_trim folder exist"; echo "Let's move on."; date
else echo "Working cDNA_trim folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/cDNA_trim; date
fi

#--------------------------------------------------------------------------------

# Trim and filter ONT cDNA using nanoq

nanoq \
-i $RAW_READS \
-o $WORKING_FOLDER_SCRATCH/cDNA_trim/Nucella.ONT.cDNA.barcode12_clean.fastq.gz \
-l 1000 \
-q 10 \
-S 10 \
-j

# Minimum read length (-l)
# Minimum average read quality (-q)
# Fixed number of bases can be trimmed from the start (-S) 
# Output summary report in JSON format (-j)

#--------------------------------------------------------------------------------

conda deactivate
