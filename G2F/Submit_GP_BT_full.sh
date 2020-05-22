#!/bin/bash -l
#SBATCH -D /home/deruncie/projects/BSFG/G2F
#SBATCH -o /home/deruncie/projects/BSFG/G2F/logs/out-%j.txt
#SBATCH -e /home/deruncie/projects/BSFG/G2F/logs/error-%j.txt
#SBATCH -J GP_bT_full


module load R

Rscript GP_byTrait_full.R $SLURM_ARRAY_TASK_ID > logs/GP_byTrait_full_$SLURM_ARRAY_TASK_ID.Rout

