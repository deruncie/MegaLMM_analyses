# MegaLMM analyses of the Genome2Field data

## Data download and proccessing

### Phenotype data:

| File                            | Location                                                                                                                                                                     |
|---------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `g2f_2014_hybrid_no_outliers.csv` |       `https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/_deprecated_Carolyn_Lawrence_Dill_G2F_Nov_2016_V2/a._2014_hybrid_phenotypic_data`        |
| `g2f_2015_hybrid_data_clean.csv`  |       `https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/GenomesToFields_2014_2017_v1/G2F_Planting_Season_2015_v2/a._2015_hybrid_phenotypic_data` |
| `g2f_2016_hybrid_data_clean.csv`  |       `https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/GenomesToFields_2014_2017_v1/G2F_Planting_Season_2016_v2/a._2016_hybrid_phenotypic_data` |
| `g2f_2017_hybrid_data_clean.csv`  |       `https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/GenomesToFields_2014_2017_v1/G2F_Planting_Season_2017_v1/a._2017_hybrid_phenotypic_data` |

### Field data:
| File                               | Location                                                                                                                                                                |
|------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `g2f_2014_field_characteristics.csv` |       `https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/GenomesToFields_2014_2017_v1/G2F_Planting_Season_2014_v4/z._2014_supplemental_info` |
| `g2f_2015_field_metadata.csv`        |       `https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/GenomesToFields_2014_2017_v1/G2F_Planting_Season_2015_v2/z._2015_supplemental_info` |
| `g2f_2016_field_metadata.csv`        |       `https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/GenomesToFields_2014_2017_v1/G2F_Planting_Season_2016_v2/z._2016_supplemental_info` |
| `g2f_2017_field_metadata.csv`        |       `https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/GenomesToFields_2014_2017_v1/G2F_Planting_Season_2017_v1/z._2017_supplemental_info` |

### Genotype data:
Data downloaded from: `http://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/Carolyn_Lawrence_Dill_G2F_Nov_2016_V.3/b._2014_gbs_data`

Load `g2f_2014_zeaGBSv27.imp.h5` and export `g2f_2014_zeaGBSv27_TaxaList.txt` and `g2f_2014_zeaGBSv27_TaxaSummary.txt` using Data-> Get Taxa List and Data-> Geno Summary

### Prep data:
Run `G2F_prep.R`

Run `GP_byTrait_setup.R`, but stop after line 29, go to TASSEL, load the GBS data and `Data/Hybrid_genotypes.txt`, run Data-> Create Hybrid Genotypes, then compute Kinship matrix, saving the result as: `Data/g2f_2014_zeaGBSv27_CenteredIBS_allYears.txt`. Then continue with `GP_byTrait_setup.R`

## Analyses

The analysis scripts are:

- `GP_byTrait_alternatives.R` (univariate analysis and phenix)- `GP_byTrait_full.R` (MegaLMM)- `GP_byTrait_means.R` (univariate on the trait means

These can be paralelized with the corresponding `.sh` scripts, though the directories need to be adjusted.

Finally, results can be collected and figures made by running:
`GP_byTrait_results.R`