#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Quast_ntlink_array

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss --7 day limit
#SBATCH --time=2:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=50G

# Submit job array
#SBATCH --array=1-4

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# This script will run Quast on the final consensus assembly 

# Quast executable
quast=/netfiles/nunezlab/Shared_Resources/Software/quast-5.2.0/quast.py

#--------------------------------------------------------------------------------

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly

#--------------------------------------------------------------------------------

## Read guide files
# This is a guide file with all of the parameter combinations
# k = 20, 24, 30
# w = 60, 75, 100, 150
# r = 6, 8, 10

GUIDE_FILE=$WORKING_FOLDER_SCRATCH/ntlink/ntlink_guide_file_2_subset.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers.
##   k            w        r
##   20          60        6    
##   20          75        6   
##   20          100       6   
##   20          150       6  
##   24          60        6  
##   ...

#--------------------------------------------------------------------------------

# Determine parameter combination to process
k=$( cat $GUIDE_FILE  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $1 }' )
w=$( cat $GUIDE_FILE  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $2 }' )
r=$( cat $GUIDE_FILE  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $3 }' )

label=k.${k}_w.${w}_r.${r}

echo ${k} ${w} ${r} ${label}

#--------------------------------------------------------------------------------

cd $WORKING_FOLDER_SCRATCH/ntlink

# Make Quast directory 
if [ -d "Quast" ]
then echo "Working Quast folder exist"; echo "Let's move on."; date
else echo "Working Quast folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/ntlink/Quast; date
fi

cd $WORKING_FOLDER_SCRATCH/ntlink/Quast
# Make Quast directory for each parameter combination
if [ -d "ntlink_${label}" ]
then echo "Working Quast_ntlink_${label} folder exist"; echo "Let's move on."; date
else echo "Working Quast_ntlink_${label} folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/ntlink/Quast/Quast_ntlink_${label}; date
fi

#--------------------------------------------------------------------------------

# Assembly name
ASSEMBLY=$WORKING_FOLDER_SCRATCH/ntlink/ntlink_${label}/final_assembly.fasta.k${k}.w${w}.z1000.ntLink.ntLink.ntLink.ntLink.ntLink.gap_fill.fa.k${k}.w${w}.z1000.ntLink.scaffolds.gap_fill.fa

# Run quast
$quast $ASSEMBLY \
-o $WORKING_FOLDER_SCRATCH/ntlink/Quast/Quast_ntlink_${label}