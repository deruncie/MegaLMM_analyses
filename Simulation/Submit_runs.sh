#!/bin/bash -l
#SBATCH -D /home/deruncie/projects/BSFG/Method_comparison
#SBATCH -o /home/deruncie/projects/BSFG/Method_comparison/logs/out-%j.txt
#SBATCH -e /home/deruncie/projects/BSFG/Method_comparison/logs/error-%j.txt
#SBATCH -J MC


module load R

Rscript Method_comparisons.R $SLURM_ARRAY_TASK_ID > logs/Method_comparisons_$SLURM_ARRAY_TASK_ID.Rout

