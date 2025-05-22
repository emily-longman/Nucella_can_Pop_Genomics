#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Quast_SparseAssembler_array

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss --7 day limit
#SBATCH --time=8:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=50G

# Submit job array
#SBATCH --array=1-16

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# This script will run Quast on each SparseAssembler contig.txt file

# Quast executable
quast=/netfiles/nunezlab/Shared_Resources/Software/quast-5.2.0/quast.py

#--------------------------------------------------------------------------------

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/short_read_assembly

#--------------------------------------------------------------------------------

## Read guide files
# This is a guide file with all of the parameter combinations
# k = 31, 51, 71, 91
# NodeCovTh = 1, 2
# EdgeCovTh = 0, 1

GUIDE_FILE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/short_read_assembly/SparseAssembler/SparseAssembler_GuideFile.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers.
##   k       NodeCovTh     EdgeCovTh   
##   31          1            0
##   31          1            1
##   31          2            0
##   31          2            1
##   51          1            0
##   ....

#--------------------------------------------------------------------------------

# Determine parameter combination to process
k=$( cat $GUIDE_FILE  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $1 }' )
NCT=$( cat $GUIDE_FILE  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $2 }' )
ECT=$( cat $GUIDE_FILE  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $3 }' )

echo $k  $NCT  $ECT

#--------------------------------------------------------------------------------
cd $WORKING_FOLDER/SparseAssembler

# Make Quast directory 
if [ -d "Quast" ]
then echo "Working Quast folder exist"; echo "Let's move on."; date
else echo "Working Quast folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/SparseAssembler/Quast; date
fi

cd $WORKING_FOLDER/SparseAssembler/Quast

# Make Quast directory for each parameter combination
if [ -d "SparseAssembler_${k}_${NCT}_${ECT}" ]
then echo "Working SparseAssembler_${k}_${NCT}_${ECT} folder exist"; echo "Let's move on."; date
else echo "Working SparseAssembler_${k}_${NCT}_${ECT} folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/SparseAssembler/Quast/SparseAssembler_${k}_${NCT}_${ECT}; date
fi

#--------------------------------------------------------------------------------

# Run quast
$quast $WORKING_FOLDER/SparseAssembler/SparseAssembler_${k}_${NCT}_${ECT}/Contigs.txt \
-o $WORKING_FOLDER/SparseAssembler/Quast/SparseAssembler_${k}_${NCT}_${ECT}