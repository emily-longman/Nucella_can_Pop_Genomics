#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=ntLink.shasta.2000 

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --cpus-per-task=40 
#SBATCH --nodes=1 # on one node

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=28:00:00 #<= this may depend on your resources

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=95G #<= this may depend on your resources

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/ntLinkNu.2000.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

# ---------------------------

# Load ntlink
module load python3.11-anaconda/2023.09-0
source ${ANACONDA_ROOT}/etc/profile.d/conda.sh
conda create --name ntlink #create and name the environment
source activate ntlink #activate the environment
#conda install -c bioconda -c conda-forge ntlink # Only need to install once
conda activate ntlink 

#Working folder
WORKING_FOLDER=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Shasta_assemblies/ntlink

# ---------------------------

# Generate Folders and files

cd ${WORKING_FOLDER}

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "ntlink_2000" ]
then echo "Working ntlink_2000 folder exist"; echo "Let's move on."; date
else echo "Working ntlink_2000 folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/ntlink_2000; date
fi

#Output folder
OUTPUT=$WORKING_FOLDER/ntlink_2000

# Move to the directory where the output files will be saved
cd ${OUTPUT}

# ---------------------------

link_pool=/netfiles/pespenilab_share/Nucella/raw/ONT/FC_all.ONT.nuc.fastq.gz
assembly=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Shasta_assemblies/ShastaRun2000/Assembly.fasta 

# Move the assembly to the Base_genome file before scafolding
cp $assembly ${OUTPUT}

# Run ntLink_rounds
ntLink_rounds run_rounds_gaps \
target=Assembly.fasta \
reads=$link_pool k=24 w=150 t=40 rounds=8

# Deactivate conda
conda deactivate

echo "done"
