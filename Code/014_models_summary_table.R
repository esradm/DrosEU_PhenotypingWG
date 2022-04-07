
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
p_pop <- read.csv("LinearModelsPop/all_lmers_pop_pvalues.csv")
p_lat <- read.csv("LinearModelsLat/all_lmers_lat_pvalues.csv")
p_lon <- read.csv("LinearModelsLon/all_lmers_lon_pvalues.csv")
p_alt <- read.csv("LinearModelsAlt/all_lmers_alt_pvalues.csv")
p_metas <- read.csv("MetaAnalyses/all_metas_pvalues.csv")

##### create output directory
sum_dir <- "SummaryLmerMeta"
dir.create(sum_dir, showWarnings = F) 



p_lmers <- inner_join(p_pop, p_lat) %>% 
  inner_join(p_lon) %>% 
  inner_join(p_alt)

p_lmers$Lab <- str_match(p_lmers$Trait, '([^_]+)(?:_[^_]+){0}$')[,1]
p_lmers$Trait <- str_match(p_lmers$Trait, '(.*[^_]+)(?:_[^_]+){1}$')[,2]

p_lmers <- relocate(p_lmers, Lab)

p_metas <- dplyr::select(p_metas, Meta, P_bonf) %>%
  dplyr::rename(Trait = Meta, P_meta = P_bonf) %>%
  mutate(Trait = str_replace(Trait, "_meta", ""))


p <- left_join(p_lmers, p_metas) %>%
  arrange(Trait, Lab)

write.csv(p, file.path(sum_dir, "lmers_metas_pvalues_summary_table.csv"), row.names = F)


