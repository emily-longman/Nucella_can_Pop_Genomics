#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=1:00:00
#SBATCH -o ./SlurmOut/slim.%A_%a.out # Standard output
#SBATCH -p general
#SBATCH --array=1-500

repid=${SLURM_ARRAY_TASK_ID}
root=/gpfs2/scratch/jcnunez/Nucella_scra/slim_glacial_refugiaSouthExpansion
mkdir slim_glacial_refugiaSouthExpansion

module load slim

#####
    
    slim \
	-d "repId=${repid}" \
	-d "root='${root}'" \
     2.Glaciation_refugiaSouthExpansion.slim
     
