#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=flye_Nucella 

# Specify partition
#SBATCH --partition=bigmemwk

# Request nodes
#SBATCH --cpus-per-task=1 
#SBATCH --nodes=1 # on one node

#SBATCH -c 40 #number of cores

# Reserve walltime --time=<dd-hh:mm:ss>
#SBATCH --time=07-00:00:00 #<= this may depend on your resources

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=900G #<= this may depend on your resources

# Submit job array
#SBATCH --array=0-2

# Name output of this job using %x=job-name and %j=job-id
#SBATCH -o ./slurmOutput/flye.%A_%a.out # Standard output
#SBATCH -e ./slurmOutput/flye.%A_%a.err # Standard error

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

### Make flye executable
flye=/gpfs1/home/e/l/elongman/software/Flye/bin/flye

## reads
arr=(15000 25000 50000)
L="${arr[$SLURM_ARRAY_TASK_ID]}"
echo $L

# Move to the directory
cd /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/ONT_fltlong

# Unzip files
gunzip Nuc.*.fastq.gz

# Move to the directory where the output files will be saved
cd /netfiles/pespenilab_share/Nucella/processed/Base_Genome/Flye_assembly

#input
ont=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/ONT_fltlong/Nuc.$L.fltlong.fastq
echo $ont

## out folder
out=./nucella_${L}

###
#size=2.5g
cpu=40

###run flye
python $flye --nano-raw $ont \
--out-dir $out \
--threads $cpu \
--no-alt-contigs 


echo "done"