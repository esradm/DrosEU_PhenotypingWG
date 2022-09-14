
##########################################################################
######################  LMER and META SUMMARY TABLE ######################
##########################################################################




##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)

##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")

##### load data
p_pop <- read.csv("LinearModelsPop/all_models_pop_pvalues.csv")
p_lat <- read.csv("LinearModelsLat/all_models_lat_pvalues.csv")
p_lon <- read.csv("LinearModelsLon/all_models_lon_pvalues.csv")
p_alt <- read.csv("LinearModelsAlt/all_models_alt_pvalues.csv")
p_coxme <-  read.csv("SurvivalAnalyses/all_coxmes_pop_pvalues.csv")
p_metas <- read.csv("MetaAnalyses/all_models_pop_meta_pvalues.csv")


##### create output directory
sum_dir <- "InfoTables"
dir.create(sum_dir, showWarnings = F) 

##### combine p values
p <- bind_rows(p_pop, p_lat, p_lon, p_alt, p_coxme)
p$Model[p$Model == "lm"] <- "lmer" # quick fix
p <- pivot_wider(p, names_from = Predictor, values_from = P, names_prefix = "P_")
p <- left_join(p,  p_metas %>%
                 dplyr::select(-c(Q, P, Min_lab, Max_lab, P_bh)) %>%
                 dplyr::rename(P_meta_adj = P_bonf, Model = Models) %>%
                 mutate(Model = str_replace(Model, "lmers", "lmer")))
p <- arrange(p, desc(Model), Trait, Lab, Sex)


##### export p values
write.csv(p, file.path(sum_dir, "models_metas_pvalues_summary_table.csv"), row.names = F)


