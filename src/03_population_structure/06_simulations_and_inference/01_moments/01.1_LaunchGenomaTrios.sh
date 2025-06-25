#!/bin/sh
#
#SBATCH -J genoTrios
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=5:00:00
#SBATCH -o ./slurmOut/genom.%A_%a.out # Standard output
#SBATCH -p general
#SBATCH --array=1-3

module load Rgeospatial

guide=/gpfs2/scratch/jcnunez/Nucella_scra/trios_guide.txt

Rscript --vanilla 1.Genomalicious.Trios.Nucella.R \
${SLURM_ARRAY_TASK_ID} \
$guide
