


############################################################
###################### HERITABILITIES ######################
############################################################




##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(lme4)
#library(lsmeans)
#library(afex)
#library(multcomp)



##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")


##### source functions
#source("Functions/lab_correlations_functions.R")

##### load data
droseu <- readRDS("Data/droseu_master_list_2022-05-02.rds")

##### create output directory
#desc_dir <- "DescriptiveStatistics"
#dir.create(desc_dir, showWarnings = F) 


##### functions

  
H2_lmm <- function(data, phenot, genot) {
  if(length(unlist(unique(data[,genot]))) == nrow(data)){
    H2 <- "NA"
  } else {
  lmm <- lmer(as.formula(paste0(phenot, "~1 + 1|", genot)), data = data)	
  lmm.varcor <- as.data.frame(summary(lmm)$varcor)
  Vg <- lmm.varcor$vcov[1]
  Ve <- lmm.varcor$vcov[2]
  Vp <- Vg + Ve
  H2 <- Vg/Vp
  H2
  }
}


#####


via <- group_by(droseu$via, Supervisor.PI, Line, ReplicateVial) %>%
  summarise(via = mean(ProportionEggtoAdultSurvival))

h2_via <- H2_lmm(via, "via", "Line")

via_pi <- split(via, f = via$Supervisor.PI)

h2_wal_pi <- lapply(via_pi, H2_lmm, "via", "Line")




##### 

wa <- group_by(droseu$wa, Supervisor.PI, Line, ReplicateVial) %>%
  summarise(WAL = mean(CentroidSizeLeft_micrometers), WAR = mean(CentroidSizeRight_micrometers))

h2_wal <- H2_lmm(wa, "WAL", "Line")
h2_war <- H2_lmm(wa, "WAR", "Line")

wa_pi <- split(wa, f = wa$Supervisor.PI)

h2_wal_pi <- lapply(wa_pi, H2_lmm, "WAL", "Line")
h2_war_pi <- lapply(wa_pi, H2_lmm, "WAR", "Line")


