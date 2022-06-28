

# This code sources in the right order all the code chunks needed to reproduce all the analyses that have been performed so far.

##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")


##### clean and reformat raw data
source("Code/010_data_reformating.R")

##### compute descriptive statistics
source("Code/020_descriptive_statistics.R")

##### runs linear models 
source("Code/030_linear_models_pop.R")
source("Code/031_linear_models_pop_lines_estimates.R")
source("Code/032_linear_models_lat.R")
source("Code/033_linear_models_lon.R")
source("Code/034_linear_models_alt.R")

##### runs meta analyses
source("Code/040_meta_analyses_pop.R")

##### runs cox me models for a subset of traits
#source("Code/050_survival_analyses.R")

##### plots survival curves for a subset of traits
source("Code/051_survival_curves.R")

##### combines p values from linear models and meta analyses 
source("Code/060_models_summary_table.R")

##### makes lab correlations scatterplot matrices 
source("Code/070_lab_correlations_pop.R")

