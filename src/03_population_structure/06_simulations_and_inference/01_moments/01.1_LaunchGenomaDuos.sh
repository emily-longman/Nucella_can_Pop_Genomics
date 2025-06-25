#!/bin/sh
#
#SBATCH -J genoDuo
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=5:00:00
#SBATCH -o ./slurmOut/genom.%A_%a.out # Standard output
#SBATCH -p general
#SBATCH --array=1-19

module load Rgeospatial

guide=/gpfs2/scratch/jcnunez/Nucella_scra/Duos_guide.txt

Rscript --vanilla 1.Genomalicious.Duos.Nucella.R \
${SLURM_ARRAY_TASK_ID} \
$guide
