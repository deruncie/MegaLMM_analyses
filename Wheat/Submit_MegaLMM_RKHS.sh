#!/bin/bash -l
#SBATCH -D /home/deruncie/projects/MegaLMM/Krause
#SBATCH -o /home/deruncie/projects/MegaLMM/Krause/logs/out-%j.txt
#SBATCH -e /home/deruncie/projects/MegaLMM/Krause/logs/error-%j.txt
#SBATCH -J MegaLMM_R


module load R

Rscript MegaLMM_Krause_RKHS_collectResults.R $SLURM_ARRAY_TASK_ID > logs/MegaLMM_Krause_RKHS_collectResults_$SLURM_ARRAY_TASK_ID.Rout

