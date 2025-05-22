#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Quast_DBG2OLC

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss --7 day limit
#SBATCH --time=8:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=50G

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

# Determine parameter combination to process
#kCT=$( cat $GUIDE_FILE  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $1 }' )
#MO=$( cat $GUIDE_FILE  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $2 }' )
#AT=$( cat $GUIDE_FILE  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $3 }' )

AT=0.015
MO=20
kCT=2

echo ${kCT}  ${MO}  ${AT}

#--------------------------------------------------------------------------------
cd $WORKING_FOLDER/DBG2OLC

# Make Quast directory 
if [ -d "Quast" ]
then echo "Working Quast folder exist"; echo "Let's move on."; date
else echo "Working Quast folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/DBG2OLC/Quast; date
fi

cd $WORKING_FOLDER/DBG2OLC/Quast

# Make Quast directory for each parameter combination
if [ -d "DBG2OLC_${kCT}_${MO}_${AT}" ]
then echo "Working DBG2OLC_${kCT}_${MO}_${AT} folder exist"; echo "Let's move on."; date
else echo "Working DBG2OLC_${kCT}_${MO}_${AT} folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/DBG2OLC/Quast/DBG2OLC_${kCT}_${MO}_${AT}; date
fi


#--------------------------------------------------------------------------------

# Run quast
$quast $WORKING_FOLDER/DBG2OLC/DBG2OLC_${kCT}_${MO}_${AT}/backbone_raw.fasta \
-o $WORKING_FOLDER/DBG2OLC/Quast/DBG2OLC_${kCT}_${MO}_${AT}