


############################################################
###################### HERITABILITIES ######################
############################################################



###### some things to consider
# for between labs H2, instead of using mean Line values, one could consider using Line random effects extracted from individual labs linear mixed models



##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(lme4)

##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")


##### load data
droseu <- readRDS("Data/droseu_master_list_2022-05-02.rds")
line_re <- read.csv("LinearModelsPop/all_models_line_random_effects.csv")

##### functions

# this function, based on a linear mixed model (lmer) in which genotypes (i.e. lines) are included as random effect variables, allows to extract genetic and environmental variances to compute broad sense heritability, H2.
# usage: H2_lmm(data = via_df, phenot = "via", genot = "Line")
# data: data.frame with columns for trait values and genotypes
# phenot: trait values
# genot: genotype ids

H2_lmm <- function(data, phenot, genot) {
  require("lme4")
  # make sure that H can be computed by checking that genotypes have been measured more than once. If not the case H cannot be computed and NA is returned.
  if(length(unlist(unique(data[,genot]))) == nrow(data)){ 
    H2 <- "NA"
  } else {
    # build and run the lmer model with genot as random effect variable
    lmm <- lmer(as.formula(paste0(phenot, "~1 + 1|", genot)), data = data)	
    # extract variance components
    lmm.varcor <- as.data.frame(summary(lmm)$varcor)
    # variance explained by genot
    Vg <- lmm.varcor[lmm.varcor$grp == genot, "vcov"]
    # variance explained by environment (residuals)
    Ve <- lmm.varcor[lmm.varcor$grp == "Residual", "vcov"]
    # total variance
    Vp <- Vg + Ve
    # calculate H2
    H2 <- Vg/Vp; H2
  }
}


##### example with viability

### calculate within lab H2 using Line values. This works only for labs that have several measurements per line (all labs except Schmidt)

# split to a list by PI
via_df_pi <- split(droseu$via, f = droseu$via$Supervisor.PI)

# apply H2_lmm to each element of the list
h2_via_pi <- unlist(lapply(via_df_pi, H2_lmm, "ProportionEggtoAdultSurvival", "Line"))
h2_via_pi <- data.frame(Lab = names(h2_via_pi), H2 = h2_via_pi)


### calculate between labs H2 with per lab mean Line values

# aggregate line vials to get mean line values for each lab
via_df_line <- group_by(droseu$via, Supervisor.PI, Line) %>%
  summarise(via = mean(ProportionEggtoAdultSurvival), n = n())

# apply H2_lmm
h2_via <- H2_lmm(via_df_line, "via", "Line")


### calculate between labs H2 with Line random effects
h2_via_re <- H2_lmm(filter(line_re, Trait == "Via"), "Estimate", "Line")






##### deprecated

#wa <- group_by(droseu$wa, Supervisor.PI, Line, ReplicateVial) %>%
#  summarise(WAL = mean(CentroidSizeLeft_micrometers), WAR = mean(CentroidSizeRight_micrometers), n = n())

#wa_line <- group_by(wa, Supervisor.PI, Line)%>%
#  summarise(WAL = mean(WAL), WAR = mean(WAR), n = n())

#h2_wal <- H2_lmm(wa_line, "WAL", "Line")
#h2_war <- H2_lmm(wa_line, "WAR", "Line")

#wa_pi <- split(wa, f = wa$Supervisor.PI)

#h2_wal_pi <- lapply(wa_pi, H2_lmm, "WAL", "Line")
#h2_war_pi <- lapply(wa_pi, H2_lmm, "WAR", "Line")


#lmer(DW_micrograms ~ Population + (1|Line), data = droseu$dw)

