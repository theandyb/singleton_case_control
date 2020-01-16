#!/bin/sh
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --ntasks=1
#SBATCH --mem=3000MB
#SBATCH --time=01:10:00
#SBATCH --job-name=sampleSingletons
#SBATCH --partition=nomosix
#SBATCH --array=1-6
#SBATCH --requeue

srun Rscript /net/snowwhite/home/beckandy/research/singleton_case_control/src/global_7mer.R ${SLURM_ARRAY_TASK_ID}
