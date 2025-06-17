#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=DBG2OLC.rerun

# Specify partition
#SBATCH --partition=bigmemwk

# Request nodes
#SBATCH --nodes=3 
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss --7 day limit
#SBATCH --time=07-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=900G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# DBG2OLC executable
DBG2OLC=/gpfs1/home/e/l/elongman/software/DBG2OLC

#--------------------------------------------------------------------------------

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/short_read_assembly

# Long reads
#ONT_FOLDER=/netfiles/pespenilab_share/Nucella/raw/ONT
#If you haven't done it yet, unzip the files 
#gunzip $ONT_FOLDER/FC_all.ONT.nuc.fastq.gz
#If you haven't done it yet, unzip the files 
ONT=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/ONT_fltlong/Nuc.2000.fltlong.fastq

#--------------------------------------------------------------------------------

### Read guide files
# This is a guide file with all of the parameter combinations
# kmerCovTh = 2, 5, 10
# MinOverlap = 30, 50, 100, 150
# AdaptiveTh = 0.001, 0.01

#GUIDE_FILE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/short_read_assembly/DBG2OLC/DBG2OLC_GuideFile.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers.
##   kmerCovTh   MinOverlap       AdaptiveTh   
##   2             30               0.001
##   2             50               0.001
##   2             100              0.001
##   2             150              0.001
##   2             30                0.01
##   2             50                0.01
##   ...
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

# Make directory for each DBG2OLC parameter combination
if [ -d "DBG2OLC_${kCT}_${MO}_${AT}_rerun" ]
then echo "Working DBG2OLC_${kCT}_${MO}_${AT}_rerun folder exist"; echo "Let's move on."; date
else echo "Working DBG2OLC_${kCT}_${MO}_${AT}_rerun folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/DBG2OLC/DBG2OLC_${kCT}_${MO}_${AT}_rerun; date
fi

#--------------------------------------------------------------------------------

# Move to the directory where the output files will be saved
cd $WORKING_FOLDER/DBG2OLC/DBG2OLC_${kCT}_${MO}_${AT}_rerun

# Using the SparseAssembler with the largest N50 (SparseAssembler_101_2_1)
# k= 101, NodeCovTh = 2, EdgeCovTh =1 (the last parameter didn't have an effect)

# Use DBG2OLC to construct short but accurate contigs  
$DBG2OLC \
k 17 \
AdaptiveTh ${AT} \
KmerCovTh ${kCT} \
MinOverlap ${MO} \
RemoveChimera 1 \
Contigs $WORKING_FOLDER/SparseAssembler/SparseAssembler_101_2_1/Contigs.txt \
f $ONT