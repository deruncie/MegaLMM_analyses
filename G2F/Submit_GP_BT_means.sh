#!/bin/bash -l
#SBATCH -D /home/deruncie/projects/BSFG/G2F
#SBATCH -o /home/deruncie/projects/BSFG/G2F/logs/out-%j.txt
#SBATCH -e /home/deruncie/projects/BSFG/G2F/logs/error-%j.txt
#SBATCH -J GP_bT_mean


module load R

Rscript GP_byTrait_means.R $SLURM_ARRAY_TASK_ID > logs/GP_byTrait_means_$SLURM_ARRAY_TASK_ID.Rout

