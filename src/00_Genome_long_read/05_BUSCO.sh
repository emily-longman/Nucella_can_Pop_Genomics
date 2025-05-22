#!/usr/bin/env bash
#
#SBATCH -J busco # A single job name for the array
#SBATCH -c 1
#SBATCH -N 1 # on one node
#SBATCH -t 72:00:00 #<= this may depend on your resources
#SBATCH --mem 80G #<= this may depend on your resources
#SBATCH -o ./slurmOutput/busco.%A_%a.out # Standard output
#SBATCH -p bluemoon
#SBATCH --array=1-2


source /gpfs1/home/j/c/jcnunez/miniconda3/etc/profile.d/conda.sh

## Load most up to date version 
#module load python3.11-anaconda/2023.09-0
#module list
## create 
#source ${ANACONDA_ROOT}/etc/profile.d/conda.sh

conda activate busco_env

INPUT=./Assembly.fasta.k24.w150.z1000.ntLink.8rounds.fa


if [ ${SLURM_ARRAY_TASK_ID} == 1 ]
then
echo "eval 1"
OUTPUT=./NucV2_eu
LINEAGE=/netfiles/nunezlab/BUSCO_Lineages/busco_downloads/lineages/eukaryota_odb10
elif [ ${SLURM_ARRAY_TASK_ID} == 2 ]
then
echo "eval 2"
OUTPUT=./NucV2_mol
LINEAGE=/netfiles/nunezlab/BUSCO_Lineages/busco_downloads/lineages/mollusca_odb10
fi

busco -m genome -i $INPUT -o $OUTPUT -l $LINEAGE