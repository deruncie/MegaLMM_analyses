# MegaLMM Benchmarking analyses

## Data download and processing
The following files were downloaded from: `https://data.cimmyt.org/dataset.xhtml?persistentId=hdl:11529/10548109`

- `Krause_et_al_2018_Yield_iid_BLUPs.csv`- `Krause_et_al_2018_Yield_BLUEs.csv`- `Krause_et_al_2018_Hyper_BLUEs_Individual_Time_Points.csv`- `Krause_et_al_2018_Genotypes.csv`

Calculation of the kinship matrices was done with the script: `data_prep_Krause.R` (note the commented lines must be run one time to generate the K matrix)

## Analyses

The following scripts run the analyses:

- `H_mat_Krause.R`- `MegaLMM_Krause_geno.R`- `MegaLMM_Krause_RKHS.R`- `MegaLMM_Krause.R`- `univariate_GP.R`

These can be parallelized on a cluster using the `*.sh` files.