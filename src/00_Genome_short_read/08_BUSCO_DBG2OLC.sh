#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=busco

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 # on one node
#SBATCH --ntasks-per-node=1  

# Request CPUs per task
#SBATCH -c 1

# Reserve walltime -- hh:mm:ss --30 hour limit
#SBATCH --time=30:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=60G

# Submit job array
#SBATCH --array=1-2

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/busco.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# This script runs on BUSCO on an assembly in the current directory. Thus must cd to that directory before running.

#--------------------------------------------------------------------------------

# Call package (installed with conda)
module load python3.11-anaconda/2023.09-0
source ${ANACONDA_ROOT}/etc/profile.d/conda.sh
#conda create --name busco_env #create and name the environment
conda activate busco_env #activate the environment
#conda install -c conda-forge -c bioconda busco=5.7.1 # install the program

#--------------------------------------------------------------------------------

# Run BUSCO in directory where given assembly is.

INPUT=./DBG2OLC_KmC_*_MinOv_*_Adth_*.backbone_raw.fasta

#Input for test run of DBG2OLC
#INPUT=./backbone_raw.fasta

echo $INPUT

#--------------------------------------------------------------------------------

# Run BUSCO for both the eukaryota and mollusca lineages

if [ ${SLURM_ARRAY_TASK_ID} == 1 ]
then
echo "eval 1"
OUTPUT=./Nucella_DBG2OLC_eukaryota
LINEAGE=/netfiles/nunezlab/BUSCO_Lineages/busco_downloads/lineages/eukaryota_odb10
elif [ ${SLURM_ARRAY_TASK_ID} == 2 ]
then
echo "eval 2"
OUTPUT=./Nucella_DBG2OLC_mollusca
LINEAGE=/netfiles/nunezlab/BUSCO_Lineages/busco_downloads/lineages/mollusca_odb10
fi

busco -m genome -i $INPUT -o $OUTPUT -l $LINEAGE


#--------------------------------------------------------------------------------

conda deactivate