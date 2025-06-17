#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=ntLink.consensus.6round

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --cpus-per-task=40 
#SBATCH --nodes=1 # on one node

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=28:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=60G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# Load ntlink
module load python3.11-anaconda/2023.09-0
source ${ANACONDA_ROOT}/etc/profile.d/conda.sh
conda create --name ntlink #create and name the environment
source activate ntlink #activate the environment
#conda install -c bioconda -c conda-forge ntlink # Only need to install once
conda activate ntlink 

#--------------------------------------------------------------------------------

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly

#--------------------------------------------------------------------------------

# Generate Folders and files

cd $WORKING_FOLDER_SCRATCH/consensus

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "ntlink_test_6round" ]
then echo "Working ntlink_test_6round folder exist"; echo "Let's move on."; date
else echo "Working ntlink_test_6round folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/ntlink_test_6round; date
fi

#Output folder
OUTPUT=$WORKING_FOLDER_SCRATCH/ntlink_test_6round

# Move to the directory where the output files will be saved
cd ${OUTPUT}

##--------------------------------------------------------------------------------

link_pool=/netfiles/pespenilab_share/Nucella/raw/ONT/FC_all.ONT.nuc.fastq.gz
assembly=$WORKING_FOLDER_SCRATCH/consensus/final_assembly.fasta

# Move the assembly to the Base_genome file before scafolding
cp $assembly ${OUTPUT}

# Run ntLink_rounds
ntLink_rounds run_rounds_gaps \
target=final_assembly.fasta \
reads=$link_pool k=24 w=150 t=40 rounds=6

# Deactivate conda
conda deactivate

echo "done"
