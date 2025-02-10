# CLABSI_missing_data

The project shares the code used to compare different missing data imputation methods (median/mode, LOCF, regression imputation, MICE, mixed effects model and missForest) for 7-days CLABSI prediction.

The playground directory contains R scripts used to:

create baseline models (BASE_model_MF_xx_xx_days)

create dynamic models (DYN_model_MF_xx_xx_days)

sample size calculation (sample_size_calculation)

baseline model performance summaries (summary_table_BASE)

dynamic model performance summaries (summary_table_DYN)

The R directory contains functions used to build and evaluate the models

BASE_DYN_functions contains the functions used to build and evaluate the models

config_playgorund_model_building contains some general configuration loaded before model building
