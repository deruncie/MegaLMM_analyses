#!/bin/bash -l
#SBATCH -D /home/deruncie/projects/BSFG/Method_comparison
#SBATCH -o /home/deruncie/projects/BSFG/Method_comparison/logs/out-%j.txt
#SBATCH -e /home/deruncie/projects/BSFG/Method_comparison/logs/error-%j.txt
#SBATCH -J MC


module load R

Rscript Method_comparisons_clean_changeK.R $SLURM_ARRAY_TASK_ID > logs/Method_comparisons_clean_changeK_$SLURM_ARRAY_TASK_ID.Rout

