#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=ntLink_array

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --cpus-per-task=1 
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=25:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=60G 

# Submit job array
#SBATCH --array=1-12

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/ntLink_array.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# This script will use ntLink (https://github.com/bcgsc/ntLink) to scaffold the assembly.

#--------------------------------------------------------------------------------

# Load ntlink
module load python3.11-anaconda/2024.02-1
source ${ANACONDA_ROOT}/etc/profile.d/conda.sh
#conda create --name ntlink #create and name the environment
source activate ntlink #activate the environment
#conda install -c bioconda -c conda-forge ntlink # Only need to install once
conda activate ntlink 

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

#--------------------------------------------------------------------------------

# Input files

# Raw ONT material
ONT=$WORKING_FOLDER/data/raw/ONT/FC_all.ONT.nuc.fastq.gz

# Genome from consensus
assembly=$WORKING_FOLDER/data/processed/genome_assembly/consensus/final_assembly.fasta

#--------------------------------------------------------------------------------

# Read guide files
GUIDE_FILE=$WORKING_FOLDER_SCRATCH/ntlink/ntlink_guide_file_2.txt

# Parameters:
# k = 20, 24, 30
# w = 60, 75, 100, 150

#Example: -- the headers are just for descriptive purposes. The actual file has no headers.
##   k            w        
##   20          60            
##   20          75           
##   20          100        
##   20          150         
##   24          60          
##   ...

#--------------------------------------------------------------------------------

# Determine parameter combination to process
k=$( cat $GUIDE_FILE  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $1 }' )
w=$( cat $GUIDE_FILE  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $2 }' )

label=k.${k}_w.${w}

echo ${k} ${w} ${label}

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/genome_assembly

# This part of the script will check and generate, if necessary, all of the output folders used in the script

# Make ntlink directory
if [ -d "ntlink" ]
then echo "Working ntlink folder exist"; echo "Let's move on."; date
else echo "Working ntlink folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/genome_assembly/ntlink; date
fi

# Change directory
cd $WORKING_FOLDER/data/processed/genome_assembly/ntlink

if [ -d "ntlink_${label}" ]
then echo "Working ntlink_${label} folder exist"; echo "Let's move on."; date
else echo "Working ntlink_${label} folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/genome_assembly/ntlink/ntlink_${label}; date
fi

#--------------------------------------------------------------------------------

# Specify output folder 

# Output folder
OUTPUT=$WORKING_FOLDER/data/processed/genome_assembly/ntlink/ntlink_${label}

# Move to the directory where the output files will be saved
cd ${OUTPUT}

#--------------------------------------------------------------------------------

# Move the assembly to the Base_genome file before scafolding
cp $assembly ${OUTPUT}

# Run ntLink_rounds
ntLink_rounds run_rounds_gaps \
target=final_assembly.fasta \
reads=$ONT k=${k} w=${w} t=10 rounds=6

#--------------------------------------------------------------------------------

# Deactivate conda
conda deactivate
echo "done"