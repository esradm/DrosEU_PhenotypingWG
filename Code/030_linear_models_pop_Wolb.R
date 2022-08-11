

###########################################################################
######################  RUN ALL LMMs FOR POPULATIONS ######################
###########################################################################

### Update by Martin to include also Wolbachia as a covariate and test for model fit using stepwise regression



##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(lme4)
library(lsmeans)
library(afex)
library(multcomp)
library(multcompView)
library(foreach)

## additional library stats
library(stats)

##### set working directory
#setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")
setwd("/media/inter/mkapun/projects/DrosEU_PhenotypingWG/")


##### source functions
source("Code/functions.R")

##### load data
droseu <- readRDS("Data/droseu_master_list_2022-05-02.rds")

## read Wolbachia dataset

Wolb <- read.table("Data/Wolbachia.txt",
  header=T,
  na.strings="NA")

Wolb$Line <- toupper(Wolb$Line)


Wolb.cons <- na.omit(Wolb)

##### create output directory
#lmer_dir <- "LinearModelsPop"

lmer_dir <- "LinearModelsPop_Wolb"
dir.create(lmer_dir, showWarnings = F)

############# VIABILITY #############

# create output directory
out_dir <- "Viability"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F)

# initialize output list
Via_lmers_pop <- list()

Viab <- droseu$via %>%
  inner_join(Wolb.cons,by=c("Line"))

# Gibert
Via_lmers_pop$Via_Gibert_lmer_pop <- lmer(ProportionEggtoAdultSurvival_asin ~ Population + Wolbachia +(1|Line:Population) + (1|Batch), data = filter(Viab, Supervisor.PI == "Gibert"))

# Grath, Batch is removed
Via_lmers_pop$Via_Grath_lmer_pop <- lmer(ProportionEggtoAdultSurvival_asin ~ Population + Wolbachia +(1|Line:Population), data = filter(Viab, Supervisor.PI == "Grath"))

# Hoedjes, Batch is removed because of singularity warnings
Via_lmers_pop$Via_Hoedjes_lmer_pop <- lmer(ProportionEggtoAdultSurvival_asin ~ Population + Wolbachia +(1|Line:Population), data = filter(Viab, Supervisor.PI == "Hoedjes"))

# Schmidt, LM because no Line replication
Via_lmers_pop$Via_Schmidt_lm_pop <- lm(ProportionEggtoAdultSurvival_asin ~ Population, data = filter(Viab, Supervisor.PI == "Schmidt"))

# StamenkovicRadak
Via_lmers_pop$Via_StamenkovicRadak_lmer_pop <- lmer(ProportionEggtoAdultSurvival_asin ~ Population + Wolbachia +(1|Line:Population) + (1|Batch), data = filter(Viab, Supervisor.PI == "StamenkovicRadak"))

# Zwaan, Batch is removed because of singularity warnings
Via_lmers_pop$Via_Zwaan_lmer_pop <- lmer(ProportionEggtoAdultSurvival_asin ~ Population + Wolbachia +(1|Line:Population), data = filter(Viab, Supervisor.PI == "Zwaan"))

# save output list
saveRDS(Via_lmers_pop, file = file.path(lmer_dir, out_dir, "Via_lmers_pop.rds"))





############# DEVELOPMENT TIME #############

# create output directory
out_dir <- "DevelopmentTime"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F)



DTP <- droseu$dtp %>%
  inner_join(Wolb.cons,by=c("Line"))
### Egg-to-Pupae

# initialize output list
DT_lmers_pop <- list()

# Schmidt
DT_lmers_pop$DT_P_Schmidt_lmer_pop <- lmer(DT_EggPupa ~ Population + Wolbachia +(1|Line:Population), data = filter(DTP, Supervisor.PI == "Schmidt"))


### Egg-to-Adult

DTA <- droseu$dta %>%
  inner_join(Wolb.cons,by=c("Line"))

## Females

# Gibert
DT_lmers_pop$DT_A_F_Gibert_lmer_pop <- lmer(DT_EggAdult ~ Population + Wolbachia +(1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(DTA, Supervisor.PI == "Gibert" & Sex == "F"))

# Grath, Batch and Rep removed because of singular fit
DT_lmers_pop$DT_A_F_Grath_lmer_pop <- lmer(DT_EggAdult ~ Population + Wolbachia +(1|Line:Population), data = filter(DTA, Supervisor.PI == "Grath" & Sex == "F"))

# Hoedjes
DT_lmers_pop$DT_A_F_Hoedjes_lmer_pop <- lmer(DT_EggAdult ~ Population + Wolbachia +(1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(DTA, Supervisor.PI == "Hoedjes" & Sex == "F"))

# Schmidt
DT_lmers_pop$DT_A_F_Schmidt_lmer_pop <- lmer(DT_EggAdult ~ Population + Wolbachia +(1|Line:Population), data = filter(DTA, Supervisor.PI == "Schmidt" & Sex == "F"))

# StamenkovicRadak
DT_lmers_pop$DT_A_F_StamenkovicRadak_lmer_pop <- lmer(DT_EggAdult ~ Population + Wolbachia +(1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(DTA, Supervisor.PI == "StamenkovicRadak" & Sex == "F"))

# Zwaan
DT_lmers_pop$DT_A_F_Zwaan_lmer_pop <- lmer(DT_EggAdult ~ Population + Wolbachia +(1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(DTA, Supervisor.PI == "Zwaan" & Sex == "F"))


## Males

# Gibert
DT_lmers_pop$DT_A_M_Gibert_lmer_pop <- lmer(DT_EggAdult ~ Population + Wolbachia +(1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(DTA, Supervisor.PI == "Gibert" & Sex == "M"))

# Grath, Batch and Rep removed because of singular fit
DT_lmers_pop$DT_A_M_Grath_lmer_pop <- lmer(DT_EggAdult ~ Population + Wolbachia +(1|Line:Population), data = filter(DTA, Supervisor.PI == "Grath" & Sex == "M"))

# Hoedjes
DT_lmers_pop$DT_A_M_Hoedjes_lmer_pop <- lmer(DT_EggAdult ~ Population + Wolbachia +(1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(DTA, Supervisor.PI == "Hoedjes" & Sex == "M"))

# Schmidt
DT_lmers_pop$DT_A_M_Schmidt_lmer_pop <- lmer(DT_EggAdult ~ Population + Wolbachia +(1|Line:Population), data = filter(DTA, Supervisor.PI == "Schmidt" & Sex == "M"))

# StamenkovicRadak
DT_lmers_pop$DT_A_M_StamenkovicRadak_lmer_pop <- lmer(DT_EggAdult ~ Population + Wolbachia +(1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(DTA, Supervisor.PI == "StamenkovicRadak" & Sex == "M"))

# Zwaan
DT_lmers_pop$DT_A_M_Zwaan_lmer_pop <- lmer(DT_EggAdult ~ Population + Wolbachia +(1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(DTA, Supervisor.PI == "Zwaan" & Sex == "M"))

# save output list
saveRDS(DT_lmers_pop, file = file.path(lmer_dir, out_dir, "DT_lmers_pop.rds"))







############# DRY WEIGHT #############

# create output directory
out_dir <- "DryWeight"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F)

# initialize output list
DW_lmers_pop <- list()

DW <- droseu$dw %>%
  inner_join(Wolb.cons,by=c("Line"))


## Females

# Colinet
DW_lmers_pop$DW_F_Colinet_lmer_pop <- lmer(DW_micrograms ~ Population + Wolbachia +(1|Line:Population) + (1|Batch), data = filter(DW, Supervisor.PI == "Colinet" & Sex == "F"))

# Hoedjes
DW_lmers_pop$DW_F_Hoedjes_lmer_pop <- lmer(DW_micrograms ~ Population + Wolbachia +(1|Line:Population) + (1|Batch), data = filter(DW, Supervisor.PI == "Hoedjes" & Sex == "F"))

# Onder
DW_lmers_pop$DW_F_Onder_lmer_pop <- lmer(DW_micrograms ~ Population + Wolbachia +(1|Line:Population) + (1|Batch), data = filter(DW, Supervisor.PI == "Onder" & Sex == "F"))

## Males

# Colinet, singular fit, removed Batch
DW_lmers_pop$DW_M_Colinet_lmer_pop <- lmer(DW_micrograms ~ Population + Wolbachia +(1|Line:Population), data = filter(DW, Supervisor.PI == "Colinet" & Sex == "M"))

# Hoedjes
DW_lmers_pop$DW_M_Hoedjes_lmer_pop <- lmer(DW_micrograms ~ Population + Wolbachia +(1|Line:Population) + (1|Batch), data = filter(DW, Supervisor.PI == "Hoedjes" & Sex == "M"))

# Onder
DW_lmers_pop$DW_M_Onder_lmer_pop <- lmer(DW_micrograms ~ Population + Wolbachia +(1|Line:Population) + (1|Batch), data = filter(DW, Supervisor.PI == "Onder" & Sex == "M"))

# save output list
saveRDS(DW_lmers_pop, file = file.path(lmer_dir, out_dir, "DW_lmers_pop.rds"))







############# THORAX LENGTH #############

# create output directory
out_dir <- "ThoraxLength"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F)

# initialize output list
TL_lmers_pop <- list()

TL <- droseu$tl %>%
  inner_join(Wolb.cons,by=c("Line"))


## Females

# Kozeretska
TL_lmers_pop$TL_F_Kozeretska_lmer_pop <- lmer(TL_micrometers ~ Population + Wolbachia +(1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(TL, Supervisor.PI == 'Kozeretska' & Sex == "F"))

# Posnien, Batch and Rep removed
TL_lmers_pop$TL_F_Posnien_lmer_pop <- lmer(TL_micrometers ~ Population + Wolbachia +(1|Line:Population), data = filter(TL, Supervisor.PI == 'Posnien' & Sex == "F"))

# Ritchie
TL_lmers_pop$TL_F_Ritchie_lmer_pop <- lmer(TL_micrometers ~ Population + Wolbachia +(1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(TL, Supervisor.PI == 'Ritchie' & Sex == "F"))

# Schmidt, Batch and Rep removed
TL_lmers_pop$TL_F_Schmidt_lmer_pop <- lmer(TL_micrometers ~ Population + Wolbachia +(1|Line:Population), data = filter(TL, Supervisor.PI == 'Schmidt' & Sex == "F"))

## Males

# Kozeretska
TL_lmers_pop$TL_M_Kozeretska_lmer_pop <- lmer(TL_micrometers ~ Population + Wolbachia +(1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(TL, Supervisor.PI == 'Kozeretska' & Sex == "M"))

# Posnien, Batch and Rep removed
TL_lmers_pop$TL_M_Posnien_lmer_pop <- lmer(TL_micrometers ~ Population + Wolbachia +(1|Line:Population), data = filter(TL, Supervisor.PI == 'Posnien' & Sex == "M"))

# Ritchie
TL_lmers_pop$TL_M_Ritchie_lmer_pop <- lmer(TL_micrometers ~ Population + Wolbachia +(1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(TL, Supervisor.PI == 'Ritchie' & Sex == "M"))

# save output list
saveRDS(TL_lmers_pop, file = file.path(lmer_dir, out_dir, "TL_lmers_pop.rds"))







############# WING AREA #############

# create output directory
out_dir <- "WingArea"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F)

# initialize output list
WA_lmers_pop <- list()

WA <- droseu$wa %>%
  inner_join(Wolb.cons,by=c("Line"))



## Females left

# Onder
WA_lmers_pop$WA_L_F_Onder_lmer_pop <- lmer(CentroidSizeLeft_micrometers ~ Population + Wolbachia +(1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(WA, Supervisor.PI == "Onder" & Sex == "F"))

# Posnien, removed Batch and Rep
WA_lmers_pop$WA_L_F_Posnien_lmer_pop <- lmer(CentroidSizeLeft_micrometers ~ Population + Wolbachia +(1|Line:Population), data = filter(WA, Supervisor.PI == "Posnien" & Sex == "F"))

# Ritchie
WA_lmers_pop$WA_L_F_Ritchie_lmer_pop <- lmer(CentroidSizeLeft_micrometers ~ Population + Wolbachia +(1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line:Population), data = filter(WA, Supervisor.PI == "Ritchie" & Sex == "F"))

# StamenkovicRadak, warning nearly unidentifiable model, removed Rep because it explains the least, no more warnings, same output as when "nearly unidentifiable"
WA_lmers_pop$WA_L_F_StamenkovicRadak_lmer_pop <- lmer(CentroidSizeLeft_micrometers ~ Population + Wolbachia +(1|Line:Population), data = filter(WA, Supervisor.PI == "StamenkovicRadak" & Sex == "F"))


## Males left

# Onder, failed to converge, removing Batch or Rep fixes it. Removed batch as it is the variable that explains the least
WA_lmers_pop$WA_L_M_Onder_lmer_pop <- lmer(CentroidSizeLeft_micrometers ~ Population + Wolbachia +(1|Line:Population) + (1|ReplicateVial:Line:Population), data = filter(WA, Supervisor.PI == "Onder" & Sex == "M"))

# Posnien, removed Batch and Rep
WA_lmers_pop$WA_L_M_Posnien_lmer_pop <- lmer(CentroidSizeLeft_micrometers ~ Population + Wolbachia +(1|Line:Population), data = filter(WA, Supervisor.PI == "Posnien" & Sex == "M"))

# Ritchie
WA_lmers_pop$WA_L_M_Ritchie_lmer_pop <- lmer(CentroidSizeLeft_micrometers ~ Population + Wolbachia +(1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line:Population), data = filter(WA, Supervisor.PI == "Ritchie" & Sex == "M"))

# StamenkovicRadak, warning nearly unidentifiable model, removed Rep because it explains the least, no more warnings, same output as when "nearly unidentifiable"
WA_lmers_pop$WA_L_M_StamenkovicRadak_lmer_pop <- lmer(CentroidSizeLeft_micrometers ~ Population + Wolbachia +(1|Line:Population) + (1|Batch), data = filter(WA, Supervisor.PI == "StamenkovicRadak" & Sex == "M"))


## Females right

# Onder
WA_lmers_pop$WA_R_F_Onder_lmer_pop <- lmer(CentroidSizeRight_micrometers ~ Population + Wolbachia +(1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line:Population), data = filter(WA, Supervisor.PI == "Onder" & Sex == "F"))

# Posnien, removed Batch and Rep
WA_lmers_pop$WA_R_F_Posnien_lmer_pop <- lmer(CentroidSizeRight_micrometers ~ Population + Wolbachia +(1|Line:Population), data = filter(WA, Supervisor.PI == "Posnien" & Sex == "F"))

# Ritchie
WA_lmers_pop$WA_R_F_Ritchie_lmer_pop <- lmer(CentroidSizeRight_micrometers ~ Population + Wolbachia +(1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(WA, Supervisor.PI == "Ritchie" & Sex == "F"))

# StamenkovicRadak, warning nearly unidentifiable model, removed Rep because it explains the least, no more warnings, same output as when "nearly unidentifiable"
WA_lmers_pop$WA_R_F_StamenkovicRadak_lmer_pop <- lmer(CentroidSizeRight_micrometers ~ Population + Wolbachia +(1|Line:Population), data = filter(WA, Supervisor.PI == "StamenkovicRadak" & Sex == "F"))


## Males right

# Onder, failed to converge, removing Batch or Rep fixes it. Removed Batch as it is the variable that explains the least
WA_lmers_pop$WA_R_M_Onder_lmer_pop <- lmer(CentroidSizeRight_micrometers ~ Population + Wolbachia +(1|Line:Population) + (1|ReplicateVial:Line), data = filter(WA, Supervisor.PI == "Onder" & Sex == "M"))

# Posnien, removed Batch and Rep
WA_lmers_pop$WA_R_M_Posnien_lmer_pop <- lmer(CentroidSizeRight_micrometers ~ Population + Wolbachia +(1|Line:Population), data = filter(WA, Supervisor.PI == "Posnien" & Sex == "M"))

# Ritchie
WA_lmers_pop$WA_R_M_Ritchie_lmer_pop <- lmer(CentroidSizeRight_micrometers ~ Population + Wolbachia +(1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line:Population), data = filter(WA, Supervisor.PI == "Ritchie" & Sex == "M"))

# StamenkovicRadak, warning nearly unidentifiable model, removed Rep because it explains the least, no more warnings, same output as when "nearly unidentifiable"
WA_lmers_pop$WA_R_M_StamenkovicRadak_lmer_pop <- lmer(CentroidSizeRight_micrometers ~ Population + Wolbachia +(1|Line:Population), data = filter(WA, Supervisor.PI == "StamenkovicRadak" & Sex == "M"))


# save output list
saveRDS(WA_lmers_pop, file = file.path(lmer_dir, out_dir, "WA_lmers_pop.rds"))





############# FECUNDITY #############

# create output directory
out_dir <- "Fecundity"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F)

# initialize output list
Fec_lmers_pop <- list()

FEC <- droseu$fec %>%
  inner_join(Wolb.cons,by=c("Line"))


# Billeter, no Batch
Fec_lmers_pop$Fec_Billeter_lmer_pop <- lmer(NumberOfAdultsEclosed ~ Population + Wolbachia +(1|Line:Population), data = filter(FEC, Supervisor.PI == "Billeter"))

# Fricke
Fec_lmers_pop$Fec_Fricke_lmer_pop <- lmer(NumberOfAdultsEclosed ~ Population + Wolbachia +(1|Line:Population) + (1|Batch), data = filter(FEC, Supervisor.PI == "Fricke"))


# save output list
saveRDS(Fec_lmers_pop, file = file.path(lmer_dir, out_dir, "Fec_lmers_pop.rds"))






############# LIFESPAN #############

# create output directory
out_dir <- "Lifespan"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F)

# initialize output list
LS_lmers_pop <- list()

LSL <- droseu$lsl %>%
  inner_join(Wolb.cons,by=c("Line"))


## Females

# Parsch
LS_lmers_pop$LS_F_Parsch_lmer_pop <- lmer(LSL_AgeAtDeath_days ~ Population + Wolbachia +(1|Batch) + (1|Line:Population) + (1|Line:ReplicateVial), data = filter(LSL, Censor == "0" & Supervisor.PI == "Parsch" & Sex == "F"))

# Pasyukova, failed to converge, removed Batch as it explains the least, same output as when does not converge
LS_lmers_pop$LS_F_Pasyukova_lmer_pop <- lmer(LSL_AgeAtDeath_days ~ Population + Wolbachia +(1|Line:Population) + (1|Line:ReplicateVial), data = filter(LSL, Censor == "0" & Supervisor.PI == "Pasyukova" & Sex == "F"))


## Males

# Parsch
LS_lmers_pop$LS_M_Parsch_lmer_pop <- lmer(LSL_AgeAtDeath_days ~ Population + Wolbachia +(1|Batch) + (1|Line:Population) + (1|Line:ReplicateVial), data = filter(LSL, Censor == "0" & Supervisor.PI == "Parsch" & Sex == "M"))

# Pasyukova, failed to converge, removed Batch as it explains the least, same output as when does not converge
LS_lmers_pop$LS_M_Pasyukova_lmer_pop <- lmer(LSL_AgeAtDeath_days ~ Population + Wolbachia +(1|Line:Population) + (1|Line:ReplicateVial), data = filter(LSL, Censor == "0" & Supervisor.PI == "Pasyukova" & Sex == "M"))


# save output list
saveRDS(LS_lmers_pop, file = file.path(lmer_dir, out_dir, "LS_lmers_pop.rds"))



############# COLD-SHOCK MORTALITY #############

# create output directory
out_dir <- "ColdShock"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F)

# initialize output list
CSM_lmers_pop <- list()

CSM <- droseu$csm %>%
  inner_join(Wolb.cons,by=c("Line"))


## Females

# Gonzalez
CSM_lmers_pop$CSM_F_Gonzalez_lmer_pop <- lmer(CSM_PropDead_ED_asin ~ Population + Wolbachia +(1|Line:Population) + (1|Batch), data = filter(CSM, Supervisor.PI == "Gonzalez" & Sex == "F"))

# Kozeretska
CSM_lmers_pop$CSM_F_Kozeretska_lmer_pop <- lmer(CSM_PropDead_ED_asin ~ Population + Wolbachia +(1|Line:Population) + (1|Batch), data = filter(CSM, Supervisor.PI == "Kozeretska" & Sex == "F"))

# Vieira
CSM_lmers_pop$CSM_F_Vieira_lmer_pop <- lmer(CSM_PropDead_ED_asin ~ Population + Wolbachia +(1|Line:Population) + (1|Batch), data = filter(CSM, Supervisor.PI == "Vieira" & Sex == "F"))


## Males

# Gonzalez
CSM_lmers_pop$CSM_M_Gonzalez_lmer_pop <- lmer(CSM_PropDead_ED_asin ~ Population + Wolbachia +(1|Line:Population) + (1|Batch), data = filter(CSM, Supervisor.PI == "Gonzalez" & Sex == "M"))

# Kozeretska
CSM_lmers_pop$CSM_M_Kozeretska_lmer_pop <- lmer(CSM_PropDead_ED_asin ~ Population + Wolbachia +(1|Line:Population) + (1|Batch), data = filter(CSM, Supervisor.PI == "Kozeretska" & Sex == "M"))

# Vieira
CSM_lmers_pop$CSM_M_Vieira_lmer_pop <- lmer(CSM_PropDead_ED_asin ~ Population + Wolbachia +(1|Line:Population) + (1|Batch), data = filter(CSM, Supervisor.PI == "Vieira" & Sex == "M"))


# save output list
saveRDS(CSM_lmers_pop, file = file.path(lmer_dir, out_dir, "CSM_lmers_pop.rds"))




############# CHILL-COMA RECOVERY TIME #############

# create output directory
out_dir <- "ChillComa"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F)

# initialize output list
CCRT_lmers_pop <- list()

CCRT <- droseu$ccrt %>%
  inner_join(Wolb.cons,by=c("Line"))


## Females

# Vieira
CCRT_lmers_pop$CCRT_F_Vieira_lmer_pop <- lmer(CCRT_seconds ~ Population + Wolbachia +(1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(CCRT, Censor == "0" & Supervisor.PI == "Vieira" & Sex == "F"))

# Mensch, singular fit, removed Rep
CCRT_lmers_pop$CCRT_F_Mensch_lmer_pop <- lmer(CCRT_seconds ~ Population + Wolbachia +(1|Batch) + (1|Line:Population), data = filter(CCRT, Censor == "0" & Supervisor.PI == "Mensch" & Sex == "F"))



## Males

# Vieira, singular fit, removed Batch
CCRT_lmers_pop$CCRT_M_Vieira_lmer_pop <- lmer(CCRT_seconds ~ Population + Wolbachia +(1|Line:Population) + (1|ReplicateVial:Line), data = filter(CCRT, Censor == "0" & Supervisor.PI == "Vieira" & Sex == "M"))

#  Mensch, singular fit, removed Rep
CCRT_lmers_pop$CCRT_M_Mensch_lmer_pop <- lmer(CCRT_seconds ~ Population + Wolbachia +(1|Batch) + (1|Line:Population), data = filter(CCRT, Censor == "0" & Supervisor.PI == "Mensch" & Sex == "M"))

# save output list
saveRDS(CCRT_lmers_pop, file = file.path(lmer_dir, out_dir, "CCRT_lmers_pop.rds"))





############# HEAT-SHOCK MORTALITY #############

# create output directory
out_dir <- "HeatShock"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F)

# initialize output list
HSM_lmers_pop <- list()

HSM <- droseu$hsm %>%
  inner_join(Wolb.cons,by=c("Line"))


## Females

# Parsch
HSM_lmers_pop$HSM_F_Parsch_lmer_pop <- lmer(TimeDeath_min ~ Population + Wolbachia +(1|Line:Population) + (1|Batch), data = filter(HSM, Censor == "0" & Supervisor.PI == "Parsch" & Sex == "F"))

# Vieira
HSM_lmers_pop$HSM_F_Vieira_lmer_pop <- lmer(TimeDeath_min ~ Population + Wolbachia +(1|Line:Population) + (1|Batch), data = filter(HSM, Censor == "0" & Supervisor.PI == "Vieira" & Sex == "F"))


## Males

# Parsch
HSM_lmers_pop$HSM_M_Parsch_lmer_pop <- lmer(TimeDeath_min ~ Population + Wolbachia +(1|Line:Population) + (1|Batch), data = filter(HSM, Censor == "0" & Supervisor.PI == "Parsch" & Sex == "M"))

# Vieira
HSM_lmers_pop$HSM_M_Vieira_lmer_pop <- lmer(TimeDeath_min ~ Population + Wolbachia +(1|Line:Population) + (1|Batch), data = filter(HSM, Censor == "0" & Supervisor.PI == "Vieira" & Sex == "M"))

# save output list
saveRDS(HSM_lmers_pop, file = file.path(lmer_dir, out_dir, "HSM_lmers_pop.rds"))







############# DIAPAUSE #############

# create output directory
out_dir <- "Diapause"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F)

# initialize output list
Dia_lmers_pop <- list()

DIA <- droseu$dia %>%
  inner_join(Wolb.cons,by=c("Line"))


## LMERs

# Bergland, singular fit, removed Batch
Dia_lmers_pop$Dia_Bergland_lmer_pop <- lmer(Prop_Max_Stage9_asin ~ Population + Wolbachia +(1|Line:Population), data = filter(DIA, Supervisor.PI == "Bergland"))

# Flatt
Dia_lmers_pop$Dia_Flatt_lm_pop <- lm(Prop_Max_Stage9_asin ~ Population, data = filter(DIA, Supervisor.PI == "Flatt"))

# Schlotterer
Dia_lmers_pop$Dia_Schlotterer_lmer_pop <- lmer(Prop_Max_Stage9_asin ~ Population + Wolbachia +(1|Line:Population) + (1|Batch), data = filter(DIA, Supervisor.PI == "Schlotterer"))


# save output list
saveRDS(Dia_lmers_pop, file = file.path(lmer_dir, out_dir, "Dia_lmers_pop.rds"))





## GLMERs

dia <- read.csv("Data/MasterSheets_May22_git/DIA_MasterSheet_Feb22.csv")

dia <- group_by(dia, Supervisor.PI, Diet, Batch, Population, Line) %>%
  summarise(Prop_Max_Stage9 = mean(MostAdvancedStage <= 9 & NumberOfEggs == 0),
            Prop_Max_Stage9_asin = asin(sqrt(Prop_Max_Stage9)),
            n = n())

dia <- dia %>%
  inner_join(Wolb.cons,by=c("Line"))


Dia_glmers_pop <- foreach(pi = unique(dia$Supervisor.PI)) %do% {
  m1 <- glmer(Prop_Max_Stage9 ~ Population + Wolbachia +(1|Line:Population), weights = n, family = binomial(), data = filter(dia, Supervisor.PI == pi))
  m2 <- glmer(Prop_Max_Stage9 ~ + (1|Line:Population), weights = n, family = binomial(), data = filter(dia, Supervisor.PI == pi))
  a <- anova(m1, m2)
  l <- list(m1, a)
  names(l) <- c(paste0("Dia_", pi, "_glmer_pop"), paste0("Dia_", pi, "_glmer_pop_anova"))
  l }

Dia_glmers_pop <- unlist(Dia_glmers_pop, recursive=FALSE)
Dia_glmers_pop_anova <- Dia_glmers_pop[grep("anova", names(Dia_glmers_pop))]
names(Dia_glmers_pop_anova) <- sub("_anova", "", names(Dia_glmers_pop_anova))
Dia_glmers_pop <- Dia_glmers_pop[grep("anova", names(Dia_glmers_pop), invert = T)]

saveRDS(Dia_glmers_pop, file = file.path(lmer_dir, out_dir, "Dia_glmers_pop.rds"))
saveRDS(Dia_glmers_pop_anova, file = file.path(lmer_dir, out_dir, "Dia_glmers_pop_anova.rds"))
capture.output(Dia_glmers_pop_anova, file = file.path(lmer_dir, out_dir, "Dia_glmers_pop_anova.txt"))


############# CIRCADIAN ECLOSION TIMING #############

# create output directory
out_dir <- "CircadianEclosion"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F)

# initialize output list
CET_lmers_pop <- list()

CET <- droseu$cet %>%
  inner_join(Wolb.cons,by=c("Line"))




############# LOCOMOTOR ACTIVITY #############

# create output directory
out_dir <- "Locomotor"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F)

# initialize output list
LA_lmers_pop <- list()

LA <- droseu$la %>%
  inner_join(Wolb.cons,by=c("Line"))


# Tauber
#LA_lmers_pop$LA_ND_Tauber_lmer_pop <- lmer(ND_log2 ~ Population + Wolbachia +(1|Line:Population) + (1|Batch), data = LA) # log2 transformation creates NA values, should they be removed?
# removed "-Inf", not converging, removed Batch that explains the least
LA_lmers_pop$LA_NDlog2_Tauber_lmer_pop <- lmer(ND_log2 ~ Population + Wolbachia +(1|Line:Population), data = filter(LA, ND_log2 != -Inf)) # log2 transformation creates -Inf values (2 cases), should they be removed? Removed for the time being

# untransformed, not converging, removed Batch that explains the least
#LA_lmers_pop$LA_ND_Tauber_lmer_pop <- lmer(ND ~ Population + Wolbachia +(1|Line:Population), data = LA)

LA_lmers_pop$LA_Period_Tauber_lmer_pop <- lmer(Period ~ Population + Wolbachia +(1|Line:Population), data = LA)

LA_lmers_pop$LA_CircPhase_Tauber_lmer_pop <- lmer(CircPhase ~ Population + Wolbachia +(1|Line:Population), data = LA)

# singular fit, removed Line
LA_lmers_pop$LA_AbsPhase_Tauber_lm_pop <- lm(AbsPhase ~ Population, data = LA)

LA_lmers_pop$LA_Activity_Tauber_lmer_pop <- lmer(Activity ~ Population + Wolbachia +(1|Line:Population), data = LA)


# save output list
saveRDS(LA_lmers_pop, file = file.path(lmer_dir, out_dir, "LA_lmers_pop.rds"))






############# STARVATION RESISTANCE #############

# create output directory
out_dir <- "Starvation"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F)

# initialize output list
SR_lmers_pop <- list()

SR <- droseu$sr %>%
  inner_join(Wolb.cons,by=c("Line"))


## Females

# Gonzalez
SR_lmers_pop$SR_F_Gonzalez_lmer_pop <- lmer(AgeAtDeath_hours ~ Population + Wolbachia +(1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line:Population), data = filter(SR, Supervisor.PI == "Gonzalez" & Sex == "F"))

# Onder, singular fit, removed Batch
SR_lmers_pop$SR_F_Onder_lmer_pop <- lmer(AgeAtDeath_hours ~ Population + Wolbachia +(1|Line:Population) + (1|ReplicateVial:Line:Population), data = filter(SR, Supervisor.PI == "Onder" & Sex == "F"))

# Pasyukova
SR_lmers_pop$SR_F_Pasyukova_lmer_pop <- lmer(AgeAtDeath_hours ~ Population + Wolbachia +(1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line:Population), data = filter(SR, Supervisor.PI == "Pasyukova" & Sex == "F"))

## Males

# Gonzalez
SR_lmers_pop$SR_M_Gonzalez_lmer_pop <- lmer(AgeAtDeath_hours ~ Population + Wolbachia +(1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line:Population), data = filter(SR, Supervisor.PI == "Gonzalez" & Sex == "M"))

# Onder, failed to converge, removed Batch (explain the least) to simplify model
SR_lmers_pop$SR_M_Onder_lmer_pop <- lmer(AgeAtDeath_hours ~ Population + Wolbachia +(1|Line:Population) + (1|ReplicateVial:Line:Population), data = filter(SR, Supervisor.PI == "Onder" & Sex == "M"))

# Pasyukova
SR_lmers_pop$SR_M_Pasyukova_lmer_pop <- lmer(AgeAtDeath_hours ~ Population + Wolbachia +(1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line:Population), data = filter(SR, Supervisor.PI == "Pasyukova" & Sex == "M"))

# save output list
saveRDS(SR_lmers_pop, file = file.path(lmer_dir, out_dir, "SR_lmers_pop.rds"))





############# PIGMENTATION #############

# create output directory
out_dir <- "Pigmentation"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F)

# initialize output list
Pgm_lmers_pop <- list()


PGM <- droseu$pgm %>%
  inner_join(Wolb.cons,by=c("Line"))


# Abbott
Pgm_lmers_pop$Pgm_T4_Abbott_lmer_pop <- lmer(PercT4_asin ~ Population + Wolbachia +(1|Line:Population), data = filter(PGM, Supervisor.PI == "Abbott"))

Pgm_lmers_pop$Pgm_T5_Abbott_lmer_pop <- lmer(PercT5_asin ~ Population + Wolbachia +(1|Line:Population), data = filter(PGM, Supervisor.PI == "Abbott"))

Pgm_lmers_pop$Pgm_T6_Abbott_lmer_pop <- lmer(PercT6_asin ~ Population + Wolbachia +(1|Line:Population), data = filter(PGM, Supervisor.PI == "Abbott"))

Pgm_lmers_pop$Pgm_Total_Abbott_lmer_pop <- lmer(TotalPerc_asin ~ Population + Wolbachia +(1|Line:Population), data = filter(PGM, Supervisor.PI == "Abbott"))


# Gibert
Pgm_lmers_pop$Pgm_T4_Gibert_lmer_pop <- lmer(PercT4_asin ~ Population + Wolbachia +(1|Line:Population), data = filter(PGM, Supervisor.PI == "Gibert"))

Pgm_lmers_pop$Pgm_T5_Gibert_lmer_pop <- lmer(PercT5_asin ~ Population + Wolbachia +(1|Line:Population), data = filter(PGM, Supervisor.PI == "Gibert"))

Pgm_lmers_pop$Pgm_T6_Gibert_lmer_pop <- lmer(PercT6_asin ~ Population + Wolbachia +(1|Line:Population), data = filter(PGM, Supervisor.PI == "Gibert"))

Pgm_lmers_pop$Pgm_Total_Gibert_lmer_pop <- lmer(TotalPerc_asin ~ Population + Wolbachia +(1|Line:Population), data = filter(PGM, Supervisor.PI == "Gibert"))


# Schmidt
Pgm_lmers_pop$Pgm_T4_Schmidt_lmer_pop <- lmer(PercT4_asin ~ Population + Wolbachia +(1|Line:Population), data = filter(PGM, Supervisor.PI == "Schmidt"))

Pgm_lmers_pop$Pgm_T5_Schmidt_lmer_pop <- lmer(PercT5_asin ~ Population + Wolbachia +(1|Line:Population), data = filter(PGM, Supervisor.PI == "Schmidt"))

Pgm_lmers_pop$Pgm_T6_Schmidt_lmer_pop <- lmer(PercT6_asin ~ Population + Wolbachia +(1|Line:Population), data = filter(PGM, Supervisor.PI == "Schmidt"))

Pgm_lmers_pop$Pgm_Total_Schmidt_lmer_pop <- lmer(TotalPerc_asin ~ Population + Wolbachia +(1|Line:Population), data = filter(PGM, Supervisor.PI == "Schmidt"))


# save output list
saveRDS(Pgm_lmers_pop, file = file.path(lmer_dir, out_dir, "Pgm_lmers_pop.rds"))





###############################################################################
###############################################################################
###############################################################################

######### combine all linear models into global lists, one for all LMERs and one for all GLMERs

### LMERs

# list with 2 levels where each element is a trait and sub element is a lab
all_lmers_pop <- rdsBatchReaderToList(path = lmer_dir, recursive = T, full.names = T, pattern = "_lmers_pop.rds")
# flatten the list and rename elements
all_lmers_pop <- unlist(all_lmers_pop, recursive=FALSE)
names(all_lmers_pop) <- str_split(names(all_lmers_pop), "\\.", simplify = T)[,2]
# output the list
saveRDS(all_lmers_pop, file = file.path(lmer_dir, "all_lmers_pop_list.rds"))

### GLMERs - currently for diapause only
all_glmers_pop <- rdsBatchReaderToList(path = lmer_dir, recursive = T, full.names = T, pattern = "_glmers_pop.rds")
all_glmers_pop <- unlist(all_glmers_pop, recursive=FALSE)
names(all_glmers_pop) <- str_split(names(all_glmers_pop), "\\.", simplify = T)[,2]
saveRDS(all_glmers_pop, file = file.path(lmer_dir, "all_glmers_pop_list.rds"))



######### output all model summaries, anovas and tukeys as global lists
######### anovas only for LMERs as they are already performed in the trait sections for GLMERs
######### tukeys only when applicable (eg for Population)


### LMERs

all_lmers_pop <- readRDS(file.path(lmer_dir, "all_lmers_pop_list.rds"))
# tukey
all_lmers_pop_tukey <- lapply(all_lmers_pop, compTukeyCLD)
saveRDS(all_lmers_pop_tukey, file = file.path(lmer_dir, "all_lmers_pop_tukey_list.rds"))
# summary
all_lmers_pop_summary <- lapply(all_lmers_pop, summary)
saveRDS(all_lmers_pop_summary, file = file.path(lmer_dir, "all_lmers_pop_summary_list.rds"))
# anova
all_lmers_pop_anova <- lapply(all_lmers_pop, anova)
saveRDS(all_lmers_pop_anova, file = "LinearModelsPop_Wolb/all_lmers_pop_anova_list.rds")


### GLMERs

all_glmers_pop <- readRDS(file.path(lmer_dir, "all_glmers_pop_list.rds"))
# tukey
all_glmers_pop_tukey <- lapply(all_glmers_pop, compTukeyCLD)
saveRDS(all_glmers_pop_tukey, file = file.path(lmer_dir, "all_glmers_pop_tukey_list.rds"))
# summary
all_glmers_pop_summary <- lapply(all_glmers_pop, summary)
saveRDS(all_glmers_pop_summary, file = file.path(lmer_dir, "all_glmers_pop_summary_list.rds"))
# anova
all_glmers_pop_anova <- rdsBatchReaderToList(path = lmer_dir, recursive = T, full.names = T, pattern = "_glmers_pop_anova.rds")
all_glmers_pop_anova <- unlist(all_glmers_pop_anova, recursive=FALSE)
names(all_glmers_pop_anova) <- str_split(names(all_glmers_pop_anova), "\\.", simplify = T)[,2]
saveRDS(all_glmers_pop_anova, file = file.path(lmer_dir, "all_glmers_pop_anova_list.rds"))





######### plot linear models residuals per trait

# list models outputs to keep the directory structure
lmers <- list.files(path = lmer_dir, recursive = T, full.names = T, pattern = "_lmers_pop.rds")
glmers <- list.files(path = lmer_dir, recursive = T, full.names = T, pattern = "_glmers_pop.rds")
models <- c(lmers, glmers)


# loop over models to get residuals
for (m in 1:length(models)){
  plotResiduals(models[m])
}



######### output all models summaries, anovas and tukeys by trait

# loop over models to get statistics
for (m in 1:length(models)){
  outputModelsStats(models[m])
}



######### output all models summaries, anovas and tukeys by trait and lab

# loop over models to get statistics
for (m in 1:length(models)){
  outputModelsStatsLab(models[m])
}


######### output all models estimates as a global list
######### estimates are the fitted Population values and their corresponding SE

all_models_estimates <- list()
for (trait in 1:length(models)){
  f <- models[trait]
  n <- str_match(f, '([^/]+)(?:/[^/]+){0}$')[,1]
  n <- sub("_lmers_pop.rds", "_lmer", n)
  n <- sub("_glmers_pop.rds", "_glmer", n)
  n <- tolower(n)
  m <- readRDS(f)
  e <- lapply(m, getEstSE)
  all_models_estimates[[trait]] <- combineEstSE(e)
  names(all_models_estimates)[trait] <- n
}

saveRDS(all_models_estimates, file = file.path(lmer_dir, "all_models_pop_estimates_list.rds"))
write.csv(bind_rows(all_models_estimates), file = file.path(lmer_dir, "all_models_pop_estimates_list.csv"), row.names = F)



######### output all models estimates by trait

for (i in 1:length(models)){
  f <- models[i]
  m <- readRDS(f)
  e <- lapply(m, getEstSE)
  e <- combineEstSE(e)
  e_out_rds <- sub(".rds", "_model_estimates.rds", f)
  e_out_txt <- sub(".rds", "_model_estimates.txt", f)
  saveRDS(e, file = e_out_rds)
  write.table(e, file = e_out_txt, row.names = F, quote = F, sep = "\t")
}


######### output all models P values

all_lmers_pop_anova <- readRDS(file.path(lmer_dir, "all_lmers_pop_anova_list.rds"))
all_glmers_pop_anova <- readRDS(file.path(lmer_dir, "all_glmers_pop_anova_list.rds"))
all_pop_anova <- c(all_lmers_pop_anova, all_glmers_pop_anova)

pvalues <- combinePValues(all_pop_anova)

saveRDS(pvalues, file = file.path(lmer_dir, "all_models_pop_pvalues.rds"))
write.csv(pvalues, file = file.path(lmer_dir, "all_models_pop_pvalues.csv"), row.names = F)
