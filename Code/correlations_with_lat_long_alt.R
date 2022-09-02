
###############################################################
###################### SOME CORRELATIONS ######################
###############################################################


##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)

##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")


##### load data
pops <- readRDS("InfoTables/DrosEU_Populations.rds")
droseu <- readRDS("Data/droseu_master_list_2022-05-02.rds")
pop_comp <- read.csv("MetaAnalyses/all_metas_pop_compound_estimates.csv")
lines_comp <- readRDS("LinearModelsPop/all_models_line_compound_random_effects_list.rds")



##### pop level

pop_geo <- inner_join(select(pop_comp, Trait, Sex, Population, Mstar), pops$by_lat) %>%
  group_split(Trait, Sex)


via <- filter(pop_geo, Trait == "Via")

cor.test(via$Mstar, via$Latitude)
