# MegaLMM Benchmarking analyses

## Data download and processing
Expression data downloaded from: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80744`, file: `GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv`

The kinship matrix downloaded from: `http://1001genomes.org/data/GMI-MPI/releases/v3.1/SNP_matrix_imputed_hdf5`.

Data processing was done with `Prep_data.R`

## Analyses

The script `Method_comparisons_clean.R` runs all analyses. It can be parallelized with `Submit_runs.sh`

The script `Collect_results.R` makes the figures