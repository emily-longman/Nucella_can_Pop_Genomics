#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=ntLinkNu.Flye

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --cpus-per-task=40 
#SBATCH --nodes=1 # on one node

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=28:00:00 #<= this may depend on your resources

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=65G #<= this may depend on your resources

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

# ---------------------------

### Just need to use ntLink_rounds (Bennet already helped me install it into my path)


module load python3.11-anaconda/2023.09-0
source ${ANACONDA_ROOT}/etc/profile.d/conda.sh
conda create --name ntlink #create and name the environment
source activate ntlink #activate the environment
conda install -c bioconda -c conda-forge ntlink
#conda install -c bioconda -c conda-forge btllib
conda activate ntlink 

#Try to use Joaquin's version 
#conda activate /netfiles/nunezlab/Shared_Resources/Software/miniconda3_envs/ntlink
# Got an error abox btllib and a dependency so tried to install separately
#conda install -c bioconda -c conda-forge btllib

#Output folder
OUTPUT=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Flye_assembly/ntlink_flye_25000

# Move to the directory where the output files will be saved
cd ${OUTPUT}

# ---------------------------

link_pool=/netfiles/pespenilab_share/Nucella/raw/ONT/FC_all.ONT.nuc.fastq.gz
asm=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Flye_assembly/nucella_25000/assembly.fasta 

# Move the assembly to the Base_genome file before scafolding
cp $asm ${OUTPUT}

# Run ntLink_rounds
ntLink_rounds run_rounds_gaps \
target=assembly.fasta \
reads=$link_pool k=32 w=100 t=40 rounds=10

echo "done"