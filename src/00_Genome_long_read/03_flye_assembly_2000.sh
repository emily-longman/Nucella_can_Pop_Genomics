#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=flye_Nucella_2000

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --cpus-per-task=1 
#SBATCH --nodes=1 

#SBATCH -c 40 #number of cores

# Reserve walltime --time=<dd-hh:mm:ss>
#SBATCH --time=03-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=400G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

### Make flye executable
flye=/gpfs1/home/e/l/elongman/software/Flye/bin/flye

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed

#--------------------------------------------------------------------------------

# Generate Folders and Files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

# Change directory
cd $WORKING_FOLDER

# Generating new folder 
if [ -d "flye_2000" ]
then echo "Working flye_2000 folder exist"; echo "Let's move on"; date
else echo "Working flye_2000 folder doesnt exist. Let's fix that"; mkdir $WORKING_FOLDER/flye_2000; date
fi

#--------------------------------------------------------------------------------

# Move to the directory where the output files will be saved
cd $WORKING_FOLDER/flye_2000

#input
ont=$WORKING_FOLDER/Nuc.2000.fltlong.fastq
echo $ont

## out folder
out=$WORKING_FOLDER/flye_2000

###
#size=2.5g
cpu=40

#--------------------------------------------------------------------------------

###run flye
python $flye --nano-raw $ont \
--out-dir $out \
--threads $cpu \
--no-alt-contigs 


echo "done"