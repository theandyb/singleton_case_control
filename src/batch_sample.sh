#!/bin/sh
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --ntasks=1
#SBATCH --mem=2GB
#SBATCH --time=01:00:00
#SBATCH --job-name=sampleSingletons
#SBATCH --partition=nomosix
#SBATCH --array=1-22
#SBATCH --requeue

srun python /net/snowwhite/home/beckandy/research/singleton_case_control/src/sampling.py -s /net/snowwhite/home/beckandy/research/smaug-redux/summaries/filtered/chr${SLURM_ARRAY_TASK_ID}.singletons -f /net/snowwhite/home/beckandy/research/singleton_case_control/reference_data/human_g1k_v37/chr${SLURM_ARRAY_TASK_ID}.fasta -o /net/snowwhite/home/beckandy/research/singleton_case_control/data/adaptive/chr${SLURM_ARRAY_TASK_ID}.csv ${SLURM_ARRAY_TASK_ID}