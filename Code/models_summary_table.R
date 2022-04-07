
#########################################################################
######################  RUN ALL LMMs FOR Altitude ######################
#########################################################################




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


p_lmers <- inner_join(p_pop, p_lat) %>% 
  inner_join(p_lon) %>% 
  inner_join(p_alt)


p_lmers$T <- str_match(p_lmers$Trait, '(.*[^_]+)(?:_[^_]+){1}$')[,2]


p_lmers$Lab <- 

separate(p_lmers, Trait, into = c("Lab", "Trait"), sep = "_", remove = T, extra = "merge")





separate(x,a,into = c("b","c"),sep = "_",remove = FALSE,extra = "merge")
