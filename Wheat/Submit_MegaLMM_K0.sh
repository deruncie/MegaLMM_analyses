#!/bin/bash -l
#SBATCH -D /group/runciegrp2/Projects/MegaLMM_revisions/Wheat
#SBATCH -o /group/runciegrp2/Projects/MegaLMM_revisions/Wheat/logs/out-%j.txt
#SBATCH -e /group/runciegrp2/Projects/MegaLMM_revisions/Wheat/logs/error-%j.txt
#SBATCH -J MegaLMM_Kr

set -e

module load R

Rscript MegaLMM_Krause_K_no2nd.R 10 ${SLURM_ARRAY_TASK_ID} > logs/MegaLMM_Krause_K_no2nd_10_${SLURM_ARRAY_TASK_ID}.Rout
