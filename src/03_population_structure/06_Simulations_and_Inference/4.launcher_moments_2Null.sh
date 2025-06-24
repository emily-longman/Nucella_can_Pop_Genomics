#!/bin/sh
#
#SBATCH -J null2
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=5:00:00
#SBATCH -o ./slurmOut/admix.%A_%a.out # Standard output
#SBATCH -p general
#SBATCH --array=0-18

module load python3.12-anaconda/2024.06-1
source ${ANACONDA_ROOT}/etc/profile.d/conda.sh
conda activate moments_jcbn
moment_script="3.nucella.run_moments_2Null.py"

  python $moment_script \
	${SLURM_ARRAY_TASK_ID}

conda deactivate

### Print the time
  echo "ended at"  `date`
