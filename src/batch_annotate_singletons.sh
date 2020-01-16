#!/bin/sh
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --ntasks=1
#SBATCH --mem=1GB
#SBATCH --time=00:15:00
#SBATCH --job-name=annoSing
#SBATCH --partition=nomosix
#SBATCH --array=1-22
#SBATCH --requeue

srun Rscript /net/snowwhite/home/beckandy/research/singleton_case_control/src/annotate_singletons.r ${SLURM_ARRAY_TASK_ID}
