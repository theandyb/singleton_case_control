#!/bin/sh
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --ntasks=1
#SBATCH --mem=1GB
#SBATCH --time=00:10:00
#SBATCH --job-name=sampleSingletons
#SBATCH --partition=nomosix
#SBATCH --array=1-22
#SBATCH --requeue

srun Rscript /net/snowwhite/home/beckandy/research/singleton_case_control/src/aggregate_counts.r /net/snowwhite/home/beckandy/research/singleton_case_control/data/chr${SLURM_ARRAY_TASK_ID}.csv /net/snowwhite/home/beckandy/research/singleton_case_control/data/control_tables/chr${SLURM_ARRAY_TASK_ID}.RData
