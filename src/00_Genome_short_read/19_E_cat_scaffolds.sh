#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=cat_scaffolds_pilon

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss 
#SBATCH --time=5:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=60G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly

scaffolds_dir=$WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_4/scaffolds

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_4

# Cat files together to produce a final consensus
for file in $(find ${scaffolds_dir} -name "*.polished.fasta"); do
cmd="cat ${file};"
eval $cmd
done > $WORKING_FOLDER_SCRATCH/pilon/polished_genome_round_4/polished_assembly.fasta