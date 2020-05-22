#!/bin/bash -l
#SBATCH -D /home/deruncie/projects/BSFG/Krause
#SBATCH -o /home/deruncie/projects/BSFG/Krause/logs/out-%j.txt
#SBATCH -e /home/deruncie/projects/BSFG/Krause/logs/error-%j.txt
#SBATCH -J BSFG_Kr


module load R

Rscript BSFG_Krause_collectResults.R $SLURM_ARRAY_TASK_ID > logs/BSFG_Krause_collectResults_$SLURM_ARRAY_TASK_ID.Rout

