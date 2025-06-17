#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=DBG2OLC_array

# Specify partition
#SBATCH --partition=bigmemwk

# Request nodes
#SBATCH --nodes=1 # on one node
#SBATCH --ntasks-per-node=1  

# Request CPUs per task
#SBATCH --cpus-per-task=1

# Reserve walltime -- hh:mm:ss --30 hour limit
#SBATCH --time=07-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=900G

# Submit job array
#SBATCH --array=1-12

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/DBG2OLC_array.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# DBG2OLC executable
DBG2OLC=/gpfs1/home/e/l/elongman/software/DBG2OLC
# Quast executable
quast=/netfiles/nunezlab/Shared_Resources/Software/quast-5.2.0/quast.py

#--------------------------------------------------------------------------------

#Working folder is core folder where this pipeline is being run.
OUTPUT_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly
WORKING_FOLDER=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/short_read_assembly

#--------------------------------------------------------------------------------

# Filtered long reads
ONT=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/ONT_fltlong/Nuc.2000.fltlong.fastq

# Contigs from Sparse Assembler (using the SparseAssembler with the largest N50 (SparseAssembler_101_2_1))
Contigs=$WORKING_FOLDER/SparseAssembler/SparseAssembler_101_2_1/Contigs.txt

# Extract base files from previous DBG2OLC run
compressed_ONT=$WORKING_FOLDER/DBG2OLC/DBG2OLC_2_20_0.015/ReadsInfoFrom_Nuc.2000.fltlong.fastq

#--------------------------------------------------------------------------------

### Read guide files (split up into two guide files becuase so many assembly outputs) 
#NOTE: as switch guide files will need to update the array

# This is a guide file with all of the parameter combinations
# kmerCovTh = 2, 4, 6, 8, 10
# MinOverlap = 30, 50, 100, 150
# AdaptiveTh = 0.001, 0.01, 0.02

GUIDE_FILE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/short_read_assembly/DBG2OLC/DBG2OLC_GuideFile_4.txt

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
K=$( cat $GUIDE_FILE  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $1 }' )
M=$( cat $GUIDE_FILE  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $2 }' )
A=$( cat $GUIDE_FILE  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $3 }' )

label=DBG2OLC_KmC_${K}_MinOv_${M}_Adth_${A}

echo $label $K $M $A

#--------------------------------------------------------------------------------
# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

# Make directory for each DBG2OLC parameter combination
cd $OUTPUT_FOLDER/DBG2OLC
if [ -d "${label}" ]
then echo "Working ${label} folder exist"; echo "Let's move on."; date
else echo "Working ${label} folder doesnt exist. Let's fix that."; mkdir $OUTPUT_FOLDER/DBG2OLC/${label}; date
fi

# Make Quast directory for each parameter combination
cd $OUTPUT_FOLDER/DBG2OLC/Quast
if [ -d "${label}" ]
then echo "Working ${label} folder exist"; echo "Let's move on."; date
else echo "Working ${label} folder doesnt exist. Let's fix that."; mkdir $OUTPUT_FOLDER/DBG2OLC/Quast/${label}; date
fi

#--------------------------------------------------------------------------------

# Use DBG2OLC to construct short but accurate contigs  

# Move to the directory where the output files will be saved
cd $OUTPUT_FOLDER/DBG2OLC/${label}

# Copy read information from previous run
cp $compressed_ONT ./

# Run DBG2OLC
$DBG2OLC \
LD 1 \
KmerCovTh $K \
AdaptiveTh $A \
MinOverlap $M \
RemoveChimera 1 \
Contigs $Contigs \
k 17 \
f $ONT

#--------------------------------------------------------------------------------

# Rename output backbone file with label name

echo "renaming outputs"

mv backbone_raw.fasta ${label}.backbone_raw.fasta

#--------------------------------------------------------------------------------

# Run quast
$quast $OUTPUT_FOLDER/DBG2OLC/${label}/${label}.backbone_raw.fasta \
-o $OUTPUT_FOLDER/DBG2OLC/Quast/${label}

#--------------------------------------------------------------------------------

echo "done"