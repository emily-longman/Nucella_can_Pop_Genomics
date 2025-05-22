#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=repeat_masker

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 # on one node
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss 
#SBATCH --time=6-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=60G

# Request CPU
#SBATCH --cpus-per-task=6

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# This script will use Repeat Masker to screen the genome for repeats and low complexity DNA sequences. 
# It will output an annotation of the repeats, as well as a modified version of the genome in which all annotated repeats have been masked (replaced with Ns).

# Call packages
module load python3.11-anaconda/2023.09-0
RepeatMasker=/netfiles/nunezlab/Shared_Resources/Software/RepeatMasker/RepeatMasker

#--------------------------------------------------------------------------------

#Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly

#This is the location where the reference genome. (note: copied a final version from pilon to repeatmasker directory)
REFERENCE=$WORKING_FOLDER_SCRATCH/repeatmasker/polished_assembly.fasta

#--------------------------------------------------------------------------------

# Closest species in famdb library

#cd /netfiles/nunezlab/Shared_Resources/Software/RepeatMasker
#python famdb.py -i ./Libraries/famdb lineage -d "Spiralia"

# Closest species is Littorina saxatilis

#--------------------------------------------------------------------------------

# Change directory
cd $WORKING_FOLDER_SCRATCH/repeatmasker

# Use RepeatMasker to mask repeats

#RepeatMasker [-options] <seqfiles(s) in fasta format>
$RepeatMasker \
-e nhmmer \
-pa 6 \
-gff \
-species "Littorina saxatilis" \
$REFERENCE

# Note: maybe add -dir in future so doesn't overlap - also can input a masked genome and do a second round
# RepeatMasker -pa 8 -e ncbi -species tetrapoda -dir Tetrapoda_mask Boa_constrictor_SGA_7C.scaffolds.BovB.fa

# Sequence comparison are performed by NHMMEr - a profile Hidden Markov Model aligner
#The script creates a .gff file with the annotation in 'General Feature Finding' format. 