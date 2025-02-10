# CLABSI_missing_data

The project shares the code used to compare different missing data imputation methods (median/mode, LOCF, regression imputation, MICE, mixed effects model and missForest) for 7-days CLABSI prediction.

The playground directory contains R scripts used to:

impute missing data without incorporating missing indicators (xx_clean)

impute missing data with missing indicators (xx_mi)


The R directory contains functions used to build and evaluate the models

BASE_DYN_functions contains the functions used to build and evaluate the models

impute_config_clean and make_cols_binary_drop contain some general configurations and functions loaded before model building
