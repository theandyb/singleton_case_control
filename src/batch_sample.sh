#!/bin/sh
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --ntasks=1
#SBATCH --mem=800
#SBATCH --time=36:00:00
#SBATCH --job-name=sampleSingletons
#SBATCH --partition=nomosix
#SBATCH --requeue

srun python /net/snowwhite/home/beckandy/research/singleton_case_control/src/sampling.py -s /net/snowwhite/home/beckandy/research/smaug-redux/summaries/chr17_full.singletons -f /net/snowwhite/home/beckandy/research/smaug-redux/reference_data/human_g1k_v37/chr17.fasta.gz -o /net/snowwhite/home/beckandy/research/singleton_case_control/data/chr17.csv 17
