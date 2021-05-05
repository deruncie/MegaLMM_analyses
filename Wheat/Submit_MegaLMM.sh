#!/bin/bash -l
#SBATCH -D /group/runciegrp2/Projects/MegaLMM_revisions/Wheat
#SBATCH -o /group/runciegrp2/Projects/MegaLMM_revisions/Wheat/logs/out-%j.txt
#SBATCH -e /group/runciegrp2/Projects/MegaLMM_revisions/Wheat/logs/error-%j.txt
#SBATCH -J MegaLMM_Kr

set -e

module load R

Rscript prep_data_trial.R ${SLURM_ARRAY_TASK_ID}

foldid=1
for foldid in {1..5};do
	Rscript univariate_GP.R ${SLURM_ARRAY_TASK_ID} ${foldid} > logs/univariate_${SLURM_ARRAY_TASK_ID}_${foldid}.Rout
	Rscript H_mat_Krause.R ${SLURM_ARRAY_TASK_ID} ${foldid} > logs/H_mat_Krause_${SLURM_ARRAY_TASK_ID}_${foldid}.Rout
	Rscript MegaLMM_Krause_K.R ${SLURM_ARRAY_TASK_ID} ${foldid} > logs/MegaLMM_Krause_K_${SLURM_ARRAY_TASK_ID}_${foldid}.Rout
	Rscript MegaLMM_Krause_K_no2nd.R ${SLURM_ARRAY_TASK_ID} ${foldid} > logs/MegaLMM_Krause_K_no2nd_${SLURM_ARRAY_TASK_ID}_${foldid}.Rout
	Rscript MegaLMM_Krause_RKHS.R ${SLURM_ARRAY_TASK_ID} ${foldid} > logs/MegaLMM_Krause_RKHS_${SLURM_ARRAY_TASK_ID}_${foldid}.Rout
done

# for foldid in {1..5};do
# 	Rscript MegaLMM_Krause_geno.R ${SLURM_ARRAY_TASK_ID} ${foldid} > logs/MegaLMM_Krause_geno_${SLURM_ARRAY_TASK_ID}_${foldid}.Rout
# done