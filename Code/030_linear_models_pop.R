

########################################################################### 
######################  RUN ALL LMMs FOR POPULATIONS ######################
###########################################################################




##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(lme4)
library(lsmeans)
library(afex)
library(multcomp)
library(foreach)

##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")


##### source functions
source("Functions/lab_correlations_functions.R")

##### load data
droseu <- readRDS("Data/droseu_master_list_2022-05-02.rds")

##### create output directory
lmer_dir <- "LinearModelsPop"
dir.create(lmer_dir, showWarnings = F) 





############# VIABILITY #############

# create output directory
out_dir <- "Viability"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
Via_lmers_pop <- list()

# Gibert
Via_lmers_pop$Via_Gibert_lmer_pop <- lmer(ProportionEggtoAdultSurvival_asin ~ Population + (1|Line:Population) + (1|Batch), data = filter(droseu$via, Supervisor.PI == "Gibert"))

# Grath, Batch is removed
Via_lmers_pop$Via_Grath_lmer_pop <- lmer(ProportionEggtoAdultSurvival_asin ~ Population + (1|Line:Population), data = filter(droseu$via, Supervisor.PI == "Grath"))

# Hoedjes, Batch is removed because of singularity warnings
Via_lmers_pop$Via_Hoedjes_lmer_pop <- lmer(ProportionEggtoAdultSurvival_asin ~ Population + (1|Line:Population), data = filter(droseu$via, Supervisor.PI == "Hoedjes"))

# Schmidt, LM because no Line replication
Via_lmers_pop$Via_Schmidt_lm_pop <- lm(ProportionEggtoAdultSurvival_asin ~ Population, data = filter(droseu$via, Supervisor.PI == "Schmidt"))

# StamenkovicRadak
Via_lmers_pop$Via_StamenkovicRadak_lmer_pop <- lmer(ProportionEggtoAdultSurvival_asin ~ Population + (1|Line:Population) + (1|Batch), data = filter(droseu$via, Supervisor.PI == "StamenkovicRadak"))

# Zwaan, Batch is removed because of singularity warnings
Via_lmers_pop$Via_Zwaan_lmer_pop <- lmer(ProportionEggtoAdultSurvival_asin ~ Population + (1|Line:Population), data = filter(droseu$via, Supervisor.PI == "Zwaan"))

# save output list
saveRDS(Via_lmers_pop, file = file.path(lmer_dir, out_dir, "Via_lmers_pop.rds"))







############# DEVELOPMENT TIME #############

# create output directory
out_dir <- "DevelopmentTime"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 


### Egg-to-Pupae

# initialize output list
DT_lmers_pop <- list()

# Schmidt
DT_lmers_pop$DT_P_Schmidt_lmer_pop <- lmer(DT_EggPupa ~ Population + (1|Line:Population), data = filter(droseu$dtp, Supervisor.PI == "Schmidt"))


### Egg-to-Adult

## Females

# Gibert
DT_lmers_pop$DT_A_F_Gibert_lmer_pop <- lmer(DT_EggAdult ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(droseu$dta, Supervisor.PI == "Gibert" & Sex == "F"))

# Grath, Batch and Rep removed because of singular fit
DT_lmers_pop$DT_A_F_Grath_lmer_pop <- lmer(DT_EggAdult ~ Population + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "Grath" & Sex == "F"))

# Hoedjes
DT_lmers_pop$DT_A_F_Hoedjes_lmer_pop <- lmer(DT_EggAdult ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(droseu$dta, Supervisor.PI == "Hoedjes" & Sex == "F"))

# Schmidt
DT_lmers_pop$DT_A_F_Schmidt_lmer_pop <- lmer(DT_EggAdult ~ Population + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "Schmidt" & Sex == "F"))

# StamenkovicRadak
DT_lmers_pop$DT_A_F_StamenkovicRadak_lmer_pop <- lmer(DT_EggAdult ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(droseu$dta, Supervisor.PI == "StamenkovicRadak" & Sex == "F"))

# Zwaan
DT_lmers_pop$DT_A_F_Zwaan_lmer_pop <- lmer(DT_EggAdult ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(droseu$dta, Supervisor.PI == "Zwaan" & Sex == "F"))


## Males

# Gibert
DT_lmers_pop$DT_A_M_Gibert_lmer_pop <- lmer(DT_EggAdult ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(droseu$dta, Supervisor.PI == "Gibert" & Sex == "M"))

# Grath, Batch and Rep removed because of singular fit
DT_lmers_pop$DT_A_M_Grath_lmer_pop <- lmer(DT_EggAdult ~ Population + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "Grath" & Sex == "M"))

# Hoedjes
DT_lmers_pop$DT_A_M_Hoedjes_lmer_pop <- lmer(DT_EggAdult ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(droseu$dta, Supervisor.PI == "Hoedjes" & Sex == "M"))

# Schmidt
DT_lmers_pop$DT_A_M_Schmidt_lmer_pop <- lmer(DT_EggAdult ~ Population + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "Schmidt" & Sex == "M"))

# StamenkovicRadak
DT_lmers_pop$DT_A_M_StamenkovicRadak_lmer_pop <- lmer(DT_EggAdult ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(droseu$dta, Supervisor.PI == "StamenkovicRadak" & Sex == "M"))

# Zwaan
DT_lmers_pop$DT_A_M_Zwaan_lmer_pop <- lmer(DT_EggAdult ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(droseu$dta, Supervisor.PI == "Zwaan" & Sex == "M"))

# save output list
saveRDS(DT_lmers_pop, file = file.path(lmer_dir, out_dir, "DT_lmers_pop.rds"))







############# DRY WEIGHT #############

# create output directory
out_dir <- "DryWeight"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
DW_lmers_pop <- list()

## Females

# Colinet
DW_lmers_pop$DW_F_Colinet_lmer_pop <- lmer(DW_micrograms ~ Population + (1|Line:Population) + (1|Batch), data = filter(droseu$dw, Supervisor.PI == "Colinet" & Sex == "F"))

# Hoedjes
DW_lmers_pop$DW_F_Hoedjes_lmer_pop <- lmer(DW_micrograms ~ Population + (1|Line:Population) + (1|Batch), data = filter(droseu$dw, Supervisor.PI == "Hoedjes" & Sex == "F"))

# Onder
DW_lmers_pop$DW_F_Onder_lmer_pop <- lmer(DW_micrograms ~ Population + (1|Line:Population) + (1|Batch), data = filter(droseu$dw, Supervisor.PI == "Onder" & Sex == "F"))

## Males

# Colinet, singular fit, removed Batch
DW_lmers_pop$DW_M_Colinet_lmer_pop <- lmer(DW_micrograms ~ Population + (1|Line:Population), data = filter(droseu$dw, Supervisor.PI == "Colinet" & Sex == "M"))

# Hoedjes
DW_lmers_pop$DW_M_Hoedjes_lmer_pop <- lmer(DW_micrograms ~ Population + (1|Line:Population) + (1|Batch), data = filter(droseu$dw, Supervisor.PI == "Hoedjes" & Sex == "M"))

# Onder
DW_lmers_pop$DW_M_Onder_lmer_pop <- lmer(DW_micrograms ~ Population + (1|Line:Population) + (1|Batch), data = filter(droseu$dw, Supervisor.PI == "Onder" & Sex == "M"))

# save output list
saveRDS(DW_lmers_pop, file = file.path(lmer_dir, out_dir, "DW_lmers_pop.rds"))







############# THORAX LENGTH #############

# create output directory
out_dir <- "ThoraxLength"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
TL_lmers_pop <- list()

## Females

# Kozeretska
TL_lmers_pop$TL_F_Kozeretska_lmer_pop <- lmer(TL_micrometers ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(droseu$tl, Supervisor.PI == 'Kozeretska' & Sex == "F"))

# Posnien, Batch and Rep removed
TL_lmers_pop$TL_F_Posnien_lmer_pop <- lmer(TL_micrometers ~ Population + (1|Line:Population), data = filter(droseu$tl, Supervisor.PI == 'Posnien' & Sex == "F"))

# Ritchie
TL_lmers_pop$TL_F_Ritchie_lmer_pop <- lmer(TL_micrometers ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(droseu$tl, Supervisor.PI == 'Ritchie' & Sex == "F"))

# Schmidt, Batch and Rep removed
TL_lmers_pop$TL_F_Schmidt_lmer_pop <- lmer(TL_micrometers ~ Population + (1|Line:Population), data = filter(droseu$tl, Supervisor.PI == 'Schmidt' & Sex == "F"))

## Males

# Kozeretska
TL_lmers_pop$TL_M_Kozeretska_lmer_pop <- lmer(TL_micrometers ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(droseu$tl, Supervisor.PI == 'Kozeretska' & Sex == "M"))

# Posnien, Batch and Rep removed
TL_lmers_pop$TL_M_Posnien_lmer_pop <- lmer(TL_micrometers ~ Population + (1|Line:Population), data = filter(droseu$tl, Supervisor.PI == 'Posnien' & Sex == "M"))

# Ritchie
TL_lmers_pop$TL_M_Ritchie_lmer_pop <- lmer(TL_micrometers ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(droseu$tl, Supervisor.PI == 'Ritchie' & Sex == "M"))

# save output list
saveRDS(TL_lmers_pop, file = file.path(lmer_dir, out_dir, "TL_lmers_pop.rds"))







############# WING AREA ############# 

# create output directory
out_dir <- "WingArea"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
WA_lmers_pop <- list()


## Females left

# Onder
WA_lmers_pop$WA_L_F_Onder_lmer_pop <- lmer(CentroidSizeLeft_micrometers ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(droseu$wa, Supervisor.PI == "Onder" & Sex == "F"))

# Posnien, removed Batch and Rep
WA_lmers_pop$WA_L_F_Posnien_lmer_pop <- lmer(CentroidSizeLeft_micrometers ~ Population + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Posnien" & Sex == "F"))

# Ritchie
WA_lmers_pop$WA_L_F_Ritchie_lmer_pop <- lmer(CentroidSizeLeft_micrometers ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line:Population), data = filter(droseu$wa, Supervisor.PI == "Ritchie" & Sex == "F"))

# StamenkovicRadak, warning nearly unidentifiable model, removed Rep because it explains the least, no more warnings, same output as when "nearly unidentifiable"
WA_lmers_pop$WA_L_F_StamenkovicRadak_lmer_pop <- lmer(CentroidSizeLeft_micrometers ~ Population + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "StamenkovicRadak" & Sex == "F"))


## Males left

# Onder, failed to converge, removing Batch or Rep fixes it. Removed batch as it is the variable that explains the least
WA_lmers_pop$WA_L_M_Onder_lmer_pop <- lmer(CentroidSizeLeft_micrometers ~ Population + (1|Line:Population) + (1|ReplicateVial:Line:Population), data = filter(droseu$wa, Supervisor.PI == "Onder" & Sex == "M"))

# Posnien, removed Batch and Rep
WA_lmers_pop$WA_L_M_Posnien_lmer_pop <- lmer(CentroidSizeLeft_micrometers ~ Population + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Posnien" & Sex == "M"))

# Ritchie
WA_lmers_pop$WA_L_M_Ritchie_lmer_pop <- lmer(CentroidSizeLeft_micrometers ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line:Population), data = filter(droseu$wa, Supervisor.PI == "Ritchie" & Sex == "M"))

# StamenkovicRadak, warning nearly unidentifiable model, removed Rep because it explains the least, no more warnings, same output as when "nearly unidentifiable"
WA_lmers_pop$WA_L_M_StamenkovicRadak_lmer_pop <- lmer(CentroidSizeLeft_micrometers ~ Population + (1|Line:Population) + (1|Batch), data = filter(droseu$wa, Supervisor.PI == "StamenkovicRadak" & Sex == "M"))


## Females right

# Onder
WA_lmers_pop$WA_R_F_Onder_lmer_pop <- lmer(CentroidSizeRight_micrometers ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line:Population), data = filter(droseu$wa, Supervisor.PI == "Onder" & Sex == "F"))

# Posnien, removed Batch and Rep
WA_lmers_pop$WA_R_F_Posnien_lmer_pop <- lmer(CentroidSizeRight_micrometers ~ Population + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Posnien" & Sex == "F"))

# Ritchie
WA_lmers_pop$WA_R_F_Ritchie_lmer_pop <- lmer(CentroidSizeRight_micrometers ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(droseu$wa, Supervisor.PI == "Ritchie" & Sex == "F"))

# StamenkovicRadak, warning nearly unidentifiable model, removed Rep because it explains the least, no more warnings, same output as when "nearly unidentifiable"
WA_lmers_pop$WA_R_F_StamenkovicRadak_lmer_pop <- lmer(CentroidSizeRight_micrometers ~ Population + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "StamenkovicRadak" & Sex == "F"))


## Males right

# Onder, failed to converge, removing Batch or Rep fixes it. Removed Batch as it is the variable that explains the least
WA_lmers_pop$WA_R_M_Onder_lmer_pop <- lmer(CentroidSizeRight_micrometers ~ Population + (1|Line:Population) + (1|ReplicateVial:Line), data = filter(droseu$wa, Supervisor.PI == "Onder" & Sex == "M"))

# Posnien, removed Batch and Rep
WA_lmers_pop$WA_R_M_Posnien_lmer_pop <- lmer(CentroidSizeRight_micrometers ~ Population + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Posnien" & Sex == "M"))

# Ritchie
WA_lmers_pop$WA_R_M_Ritchie_lmer_pop <- lmer(CentroidSizeRight_micrometers ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line:Population), data = filter(droseu$wa, Supervisor.PI == "Ritchie" & Sex == "M"))

# StamenkovicRadak, warning nearly unidentifiable model, removed Rep because it explains the least, no more warnings, same output as when "nearly unidentifiable"
WA_lmers_pop$WA_R_M_StamenkovicRadak_lmer_pop <- lmer(CentroidSizeRight_micrometers ~ Population + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "StamenkovicRadak" & Sex == "M"))


# save output list
saveRDS(WA_lmers_pop, file = file.path(lmer_dir, out_dir, "WA_lmers_pop.rds"))





############# FECUNDITY ############# 

# create output directory
out_dir <- "Fecundity"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
Fec_lmers_pop <- list()


# Billeter, no Batch
Fec_lmers_pop$Fec_Billeter_lmer_pop <- lmer(NumberOfAdultsEclosed ~ Population + (1|Line:Population), data = filter(droseu$fec, Supervisor.PI == "Billeter"))

# Fricke
Fec_lmers_pop$Fec_Fricke_lmer_pop <- lmer(NumberOfAdultsEclosed ~ Population + (1|Line:Population) + (1|Batch), data = filter(droseu$fec, Supervisor.PI == "Fricke"))


# save output list
saveRDS(Fec_lmers_pop, file = file.path(lmer_dir, out_dir, "Fec_lmers_pop.rds"))






############# LIFESPAN ############# 

# create output directory
out_dir <- "Lifespan"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
LS_lmers_pop <- list()


## Females

# Flatt
LS_lmers_pop$LS_F_Flatt_lmer_pop <- lmer(LSP_AgeAtDeath_days ~ Population + (1|Population:ReplicateCage), data = filter(droseu$lsp, Censor == "0" & Supervisor.PI == "Flatt" & Sex == "F"))

# Parsch
LS_lmers_pop$LS_F_Parsch_lmer_pop <- lmer(LSL_AgeAtDeath_days ~ Population + (1|Batch) + (1|Line:Population) + (1|Line:ReplicateVial), data = filter(droseu$lsl, Censor == "0" & Supervisor.PI == "Parsch" & Sex == "F"))

# Pasyukova, failed to converge, removed Batch as it explains the least, same output as when does not converge
LS_lmers_pop$LS_F_Pasyukova_lmer_pop <- lmer(LSL_AgeAtDeath_days ~ Population + (1|Line:Population) + (1|Line:ReplicateVial), data = filter(droseu$lsl, Censor == "0" & Supervisor.PI == "Pasyukova" & Sex == "F"))


## Males

# Flatt
LS_lmers_pop$LS_M_Flatt_lmer_pop <- lmer(LSP_AgeAtDeath_days ~ Population + (1|Population:ReplicateCage), data = filter(droseu$lsp, Censor == "0" & Supervisor.PI == "Flatt" & Sex == "M"))

# Parsch
LS_lmers_pop$LS_M_Parsch_lmer_pop <- lmer(LSL_AgeAtDeath_days ~ Population + (1|Batch) + (1|Line:Population) + (1|Line:ReplicateVial), data = filter(droseu$lsl, Censor == "0" & Supervisor.PI == "Parsch" & Sex == "M"))

# Pasyukova, failed to converge, removed Batch as it explains the least, same output as when does not converge
LS_lmers_pop$LS_M_Pasyukova_lmer_pop <- lmer(LSL_AgeAtDeath_days ~ Population + (1|Line:Population) + (1|Line:ReplicateVial), data = filter(droseu$lsl, Censor == "0" & Supervisor.PI == "Pasyukova" & Sex == "M"))


# save output list
saveRDS(LS_lmers_pop, file = file.path(lmer_dir, out_dir, "LS_lmers_pop.rds"))



############# COLD-SHOCK MORTALITY ############# 

# create output directory
out_dir <- "ColdShock"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
CSM_lmers_pop <- list()


## Females

# Gonzalez
CSM_lmers_pop$CSM_F_Gonzalez_lmer_pop <- lmer(CSM_PropDead_ED_asin ~ Population + (1|Line:Population) + (1|Batch), data = filter(droseu$csm, Supervisor.PI == "Gonzalez" & Sex == "F"))

# Kozeretska
CSM_lmers_pop$CSM_F_Kozeretska_lmer_pop <- lmer(CSM_PropDead_ED_asin ~ Population + (1|Line:Population) + (1|Batch), data = filter(droseu$csm, Supervisor.PI == "Kozeretska" & Sex == "F"))

# Vieira
CSM_lmers_pop$CSM_F_Vieira_lmer_pop <- lmer(CSM_PropDead_ED_asin ~ Population + (1|Line:Population) + (1|Batch), data = filter(droseu$csm, Supervisor.PI == "Vieira" & Sex == "F"))


## Males

# Gonzalez
CSM_lmers_pop$CSM_M_Gonzalez_lmer_pop <- lmer(CSM_PropDead_ED_asin ~ Population + (1|Line:Population) + (1|Batch), data = filter(droseu$csm, Supervisor.PI == "Gonzalez" & Sex == "M"))

# Kozeretska
CSM_lmers_pop$CSM_M_Kozeretska_lmer_pop <- lmer(CSM_PropDead_ED_asin ~ Population + (1|Line:Population) + (1|Batch), data = filter(droseu$csm, Supervisor.PI == "Kozeretska" & Sex == "M"))

# Vieira
CSM_lmers_pop$CSM_M_Vieira_lmer_pop <- lmer(CSM_PropDead_ED_asin ~ Population + (1|Line:Population) + (1|Batch), data = filter(droseu$csm, Supervisor.PI == "Vieira" & Sex == "M"))


# save output list
saveRDS(CSM_lmers_pop, file = file.path(lmer_dir, out_dir, "CSM_lmers_pop.rds"))




############# CHILL-COMA RECOVERY TIME ############# 

# create output directory
out_dir <- "ChillComa"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
CCRT_lmers_pop <- list()


## Females

# Vieira
CCRT_lmers_pop$CCRT_F_Vieira_lmer_pop <- lmer(CCRT_seconds ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(droseu$ccrt, Censor == "0" & Supervisor.PI == "Vieira" & Sex == "F"))

# Mensch, singular fit, removed Rep
CCRT_lmers_pop$CCRT_F_Mensch_lmer_pop <- lmer(CCRT_seconds ~ Population + (1|Batch) + (1|Line:Population), data = filter(droseu$ccrt, Censor == "0" & Supervisor.PI == "Mensch" & Sex == "F"))



## Males

# Vieira, singular fit, removed Batch
CCRT_lmers_pop$CCRT_M_Vieira_lmer_pop <- lmer(CCRT_seconds ~ Population + (1|Line:Population) + (1|ReplicateVial:Line), data = filter(droseu$ccrt, Censor == "0" & Supervisor.PI == "Vieira" & Sex == "M"))

#  Mensch, singular fit, removed Rep
CCRT_lmers_pop$CCRT_M_Mensch_lmer_pop <- lmer(CCRT_seconds ~ Population + (1|Batch) + (1|Line:Population), data = filter(droseu$ccrt, Censor == "0" & Supervisor.PI == "Mensch" & Sex == "M"))

# save output list
saveRDS(CCRT_lmers_pop, file = file.path(lmer_dir, out_dir, "CCRT_lmers_pop.rds"))





############# HEAT-SHOCK MORTALITY ############# 

# create output directory
out_dir <- "HeatShock"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
HSM_lmers_pop <- list()

## Females

# Parsch
HSM_lmers_pop$HSM_F_Parsch_lmer_pop <- lmer(TimeDeath_min ~ Population + (1|Line:Population) + (1|Batch), data = filter(droseu$hsm, Censor == "0" & Supervisor.PI == "Parsch" & Sex == "F"))

# Vieira
HSM_lmers_pop$HSM_F_Vieira_lmer_pop <- lmer(TimeDeath_min ~ Population + (1|Line:Population) + (1|Batch), data = filter(droseu$hsm, Censor == "0" & Supervisor.PI == "Vieira" & Sex == "F"))


## Males

# Parsch
HSM_lmers_pop$HSM_M_Parsch_lmer_pop <- lmer(TimeDeath_min ~ Population + (1|Line:Population) + (1|Batch), data = filter(droseu$hsm, Censor == "0" & Supervisor.PI == "Parsch" & Sex == "M"))

# Vieira
HSM_lmers_pop$HSM_M_Vieira_lmer_pop <- lmer(TimeDeath_min ~ Population + (1|Line:Population) + (1|Batch), data = filter(droseu$hsm, Censor == "0" & Supervisor.PI == "Vieira" & Sex == "M"))

# save output list
saveRDS(HSM_lmers_pop, file = file.path(lmer_dir, out_dir, "HSM_lmers_pop.rds"))







############# DIAPAUSE ############# 

# create output directory
out_dir <- "Diapause"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
Dia_lmers_pop <- list()

## LMERs

# Bergland, singular fit, removed Batch
Dia_lmers_pop$Dia_Bergland_lmer_pop <- lmer(Prop_Max_Stage9_asin ~ Population + (1|Line:Population), data = filter(droseu$dia, Supervisor.PI == "Bergland"))

# Flatt
Dia_lmers_pop$Dia_Flatt_lm_pop <- lm(Prop_Max_Stage9_asin ~ Population, data = filter(droseu$dia, Supervisor.PI == "Flatt"))

# Schlotterer
Dia_lmers_pop$Dia_Schlotterer_lmer_pop <- lmer(Prop_Max_Stage9_asin ~ Population + (1|Line:Population) + (1|Batch), data = filter(droseu$dia, Supervisor.PI == "Schlotterer"))


# save output list
saveRDS(Dia_lmers_pop, file = file.path(lmer_dir, out_dir, "Dia_lmers_pop.rds"))





## GLMERs

dia <- read.csv("Data/MasterSheets_May22_git/DIA_MasterSheet_Feb22.csv")

dia <- group_by(dia, Supervisor.PI, Diet, Batch, Population, Line) %>%
  summarise(Prop_Max_Stage9 = mean(MostAdvancedStage <= 9 & NumberOfEggs == 0), 
            Prop_Max_Stage9_asin = asin(sqrt(Prop_Max_Stage9)),          
            n = n())

Dia_glmers_pop <- foreach(pi = unique(dia$Supervisor.PI)) %do% {
  m1 <- glmer(Prop_Max_Stage9 ~ Population + (1|Line:Population), weights = n, family = binomial(), data = filter(dia, Supervisor.PI == pi)) 
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





############# LOCOMOTOR ACTIVITY ############# 

# create output directory
out_dir <- "Locomotor"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
LA_lmers_pop <- list()

# Tauber
#LA_lmers_pop$LA_ND_Tauber_lmer_pop <- lmer(ND_log2 ~ Population + (1|Line:Population) + (1|Batch), data = droseu$la) # log2 transformation creates NA values, should they be removed?
# removed "-Inf", not converging, removed Batch that explains the least
LA_lmers_pop$LA_NDlog2_Tauber_lmer_pop <- lmer(ND_log2 ~ Population + (1|Line:Population), data = filter(droseu$la, ND_log2 != -Inf)) # log2 transformation creates -Inf values (2 cases), should they be removed? Removed for the time being

# untransformed, not converging, removed Batch that explains the least
#LA_lmers_pop$LA_ND_Tauber_lmer_pop <- lmer(ND ~ Population + (1|Line:Population), data = droseu$la)

LA_lmers_pop$LA_Period_Tauber_lmer_pop <- lmer(Period ~ Population + (1|Line:Population), data = droseu$la)

LA_lmers_pop$LA_CircPhase_Tauber_lmer_pop <- lmer(CircPhase ~ Population + (1|Line:Population), data = droseu$la)

# singular fit, removed Line
LA_lmers_pop$LA_AbsPhase_Tauber_lm_pop <- lm(AbsPhase ~ Population, data = droseu$la)

LA_lmers_pop$LA_Activity_Tauber_lmer_pop <- lmer(Activity ~ Population + (1|Line:Population), data = droseu$la)


# save output list
saveRDS(LA_lmers_pop, file = file.path(lmer_dir, out_dir, "LA_lmers_pop.rds"))






############# STARVATION RESISTANCE ############# 

# create output directory
out_dir <- "Starvation"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
SR_lmers_pop <- list()

## Females

# Gonzalez
SR_lmers_pop$SR_F_Gonzalez_lmer_pop <- lmer(AgeAtDeath_hours ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line:Population), data = filter(droseu$sr, Supervisor.PI == "Gonzalez" & Sex == "F"))
                                          
# Onder, singular fit, removed Batch
SR_lmers_pop$SR_F_Onder_lmer_pop <- lmer(AgeAtDeath_hours ~ Population + (1|Line:Population) + (1|ReplicateVial:Line:Population), data = filter(droseu$sr, Supervisor.PI == "Onder" & Sex == "F"))

# Pasyukova
SR_lmers_pop$SR_F_Pasyukova_lmer_pop <- lmer(AgeAtDeath_hours ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line:Population), data = filter(droseu$sr, Supervisor.PI == "Pasyukova" & Sex == "F"))

## Males

# Gonzalez
SR_lmers_pop$SR_M_Gonzalez_lmer_pop <- lmer(AgeAtDeath_hours ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line:Population), data = filter(droseu$sr, Supervisor.PI == "Gonzalez" & Sex == "M"))

# Onder, failed to converge, removed Batch (explain the least) to simplify model
SR_lmers_pop$SR_M_Onder_lmer_pop <- lmer(AgeAtDeath_hours ~ Population + (1|Line:Population) + (1|ReplicateVial:Line:Population), data = filter(droseu$sr, Supervisor.PI == "Onder" & Sex == "M"))

# Pasyukova
SR_lmers_pop$SR_M_Pasyukova_lmer_pop <- lmer(AgeAtDeath_hours ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line:Population), data = filter(droseu$sr, Supervisor.PI == "Pasyukova" & Sex == "M"))

# save output list
saveRDS(SR_lmers_pop, file = file.path(lmer_dir, out_dir, "SR_lmers_pop.rds"))





############# PIGMENTATION ############# 

# create output directory
out_dir <- "Pigmentation"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
Pgm_lmers_pop <- list()


# Abbott
Pgm_lmers_pop$Pgm_T4_Abbott_lmer_pop <- lmer(PercT4_asin ~ Population + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Abbott"))

Pgm_lmers_pop$Pgm_T5_Abbott_lmer_pop <- lmer(PercT5_asin ~ Population + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Abbott"))

Pgm_lmers_pop$Pgm_T6_Abbott_lmer_pop <- lmer(PercT6_asin ~ Population + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Abbott"))

Pgm_lmers_pop$Pgm_Total_Abbott_lmer_pop <- lmer(TotalPerc_asin ~ Population + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Abbott"))


# Gibert
Pgm_lmers_pop$Pgm_T4_Gibert_lmer_pop <- lmer(PercT4_asin ~ Population + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Gibert"))

Pgm_lmers_pop$Pgm_T5_Gibert_lmer_pop <- lmer(PercT5_asin ~ Population + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Gibert"))

Pgm_lmers_pop$Pgm_T6_Gibert_lmer_pop <- lmer(PercT6_asin ~ Population + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Gibert"))

Pgm_lmers_pop$Pgm_Total_Gibert_lmer_pop <- lmer(TotalPerc_asin ~ Population + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Gibert"))


# Schmidt
Pgm_lmers_pop$Pgm_T4_Schmidt_lmer_pop <- lmer(PercT4_asin ~ Population + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Schmidt"))

Pgm_lmers_pop$Pgm_T5_Schmidt_lmer_pop <- lmer(PercT5_asin ~ Population + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Schmidt"))

Pgm_lmers_pop$Pgm_T6_Schmidt_lmer_pop <- lmer(PercT6_asin ~ Population + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Schmidt"))

Pgm_lmers_pop$Pgm_Total_Schmidt_lmer_pop <- lmer(TotalPerc_asin ~ Population + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Schmidt"))


# save output list
saveRDS(Pgm_lmers_pop, file = file.path(lmer_dir, out_dir, "Pgm_lmers_pop.rds"))





############# output all lmers as Rdata ############# 

#all_lmers_pop <- ls()[grep("_lmers_", ls())]
#save(all_lmers_pop, file = "LinearModelsPop/all_lmers_pop.Rdata")







######### combine all linear models into one list ###########

rdsBatchReaderToList <- function(...) {
  temp = list.files(...)
  tnames <- str_split(temp, "/", simplify = T)[,3]
  tnames <- str_replace(tnames, ".rds", "")
  tlist <- lapply(temp, readRDS)
  names(tlist) <- tnames
  return(tlist)
}

### for lmers

# this is a 2 level list where each element is a trait and sub element is a lab / sex
all_lmers_pop <- rdsBatchReaderToList(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "_lmers_pop.rds")

# flatten the list
all_lmers_pop <- unlist(all_lmers_pop, recursive=FALSE)
names(all_lmers_pop) <- str_split(names(all_lmers_pop), "\\.", simplify = T)[,2]

# output the list
saveRDS(all_lmers_pop, file = "LinearModelsPop/all_lmers_list_pop.rds")

### for glmers - currently for diapause only

all_glmers_pop <- rdsBatchReaderToList(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "_glmers_pop.rds")

# flatten the list
all_glmers_pop <- unlist(all_glmers_pop, recursive=FALSE)
names(all_glmers_pop) <- str_split(names(all_glmers_pop), "\\.", simplify = T)[,2]

# output the list
saveRDS(all_glmers_pop, file = "LinearModelsPop/all_glmers_list_pop.rds")





############# check linear models residuals per trait ############# 


# list models outputs to keep the directory structure
lmers <- list.files(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "_lmers_pop.rds")
glmers <- list.files(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "_glmers_pop.rds")
models <- c(lmers, glmers)

# loop over lmers to get residuals
for (i in 1:length(models)){
  f <- models[i]
  p <- str_match(f, '(.*[^/]+)(?:/[^/]+){1}$')[,2]
  dir.create(file.path(p, "by_lab_models_residuals"), showWarnings = F) 
  m <- readRDS(f)
  n <- names(m)
  f_out <- file.path(p, "by_lab_models_residuals", n)
  qq_out_png <- paste0(f_out, "_qq_plot_residuals.png")
  hist_out_png <- paste0(f_out, "_hist_residuals.png")
  for (j in 1:length(n)){
    png(filename = qq_out_png[j], height = 2100, width = 2100, res = 300)
    qqnorm(resid(m[[j]]), main = n[j])
    qqline(resid(m[[j]]))
    dev.off()
    png(filename = hist_out_png[j], height = 2100, width = 2100, res = 300)
    hist(resid(m[[j]]), main = n[j], xlab = "Residuals")
    dev.off()
  } 
}
 





############# output all lmers summaries, anovas and tukeys as global lists ############# 


compTukeyCLD <- function(x) {
  em <- emmeans(x, pairwise ~ Population, mode = "asymp", adjust = "tukey")
  let <- cld(em, Letters = letters, alpha = 0.05) %>% 
    as.data.frame %>% mutate(Population = as.factor(Population)) %>% 
    arrange(Population) %>% dplyr::rename(cld = .group)
  list(contrasts = em, letters = let)
}

# for lmers

all_lmers_pop <- readRDS("LinearModelsPop/all_lmers_list_pop.rds")

all_lmers_pop_tukey <- lapply(all_lmers_pop, compTukeyCLD)
saveRDS(all_lmers_pop_tukey, file = "LinearModelsPop/all_lmers_list_pop_tukey.rds")

all_lmers_pop_summary <- lapply(all_lmers_pop, summary)
saveRDS(all_lmers_pop_summary, file = "LinearModelsPop/all_lmers_list_pop_summary.rds")

all_lmers_pop_anova <- lapply(all_lmers_pop, anova)
saveRDS(all_lmers_pop_anova, file = "LinearModelsPop/all_lmers_list_pop_anova.rds")


# for glmers

all_glmers_pop <- readRDS("LinearModelsPop/all_glmers_list_pop.rds")

all_glmers_pop_tukey <- lapply(all_glmers_pop, compTukeyCLD)
saveRDS(all_glmers_pop_tukey, file = "LinearModelsPop/all_glmers_list_pop_tukey.rds")

all_glmers_pop_summary <- lapply(all_glmers_pop, summary)
saveRDS(all_glmers_pop_summary, file = "LinearModelsPop/all_glmers_list_pop_summary.rds")


# anovas have already been generated at the trait level, they just need to be combined

all_glmers_pop_anova <- rdsBatchReaderToList(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "_glmers_pop_anova.rds")

all_glmers_pop_anova <- unlist(all_glmers_pop_anova, recursive=FALSE)
names(all_glmers_pop_anova) <- str_split(names(all_glmers_pop_anova), "\\.", simplify = T)[,2]

saveRDS(all_glmers_pop_anova, file = "LinearModelsPop/all_glmers_list_pop_anova.rds")






############# output all lmers summaries, anovas and tukeys by trait ############# 


# list models outputs to keep the directory structure
lmers <- list.files(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "_lmers_pop.rds")
glmers <- list.files(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "_glmers_pop.rds")
models <- c(lmers, glmers)


for (i in 1:length(models)){
  f <- models[i]
  s_out_rds <- sub(".rds", "_summary.rds", f)
  s_out_txt <- sub(".rds", "_summary.txt", f)
  a_out_rds <- sub(".rds", "_anova.rds", f)
  a_out_txt <- sub(".rds", "_anova.txt", f)
  t_out_rds <- sub(".rds", "_tukey.rds", f)
  t_out_txt <- sub(".rds", "_tukey.txt", f)
  m <- readRDS(f)
  s <- lapply(m, summary)
  saveRDS(s, file = s_out_rds)
  capture.output(s, file = s_out_txt)
  if (class(m[[1]]) != "glmerMod") {
    a <- lapply(m, anova)
    saveRDS(a, file = a_out_rds) 
    capture.output(a, file = a_out_txt) }
  t <- lapply(m, compTukeyCLD)
  saveRDS(t, file = t_out_rds)
  capture.output(t, file = t_out_txt)
}
  



############# output all lmers summaries, anovas and tukeys by trait and lab ############# 

# list models outputs to keep the directory structure
lmers <- list.files(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "_lmers_pop.rds")
glmers <- list.files(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "_glmers_pop.rds")
models <- c(lmers, glmers)

for (i in 1:length(models)){
  f <- models[i]
  p <- str_match(f, '(.*[^/]+)(?:/[^/]+){1}$')[,2]
  dir.create(file.path(p, "by_lab_txt_output"), showWarnings = F) 
  dir.create(file.path(p, "by_lab_rds_output"), showWarnings = F) 
  m <- readRDS(f)
  n <- names(m)
  f_out_txt <- file.path(p, "by_lab_txt_output", n)
  f_out_rds <- file.path(p, "by_lab_rds_output", n)
  s_out_rds <- paste0(f_out_rds, "_summary.rds")
  s_out_txt <- paste0(f_out_txt, "_summary.txt")
  a_out_rds <- paste0(f_out_rds, "_anova.rds")
  a_out_txt <- paste0(f_out_txt, "_anova.txt")
  t_out_rds <- paste0(f_out_rds, "_tukey.rds")
  t_out_txt <- paste0(f_out_txt, "_tukey.txt")
  for (j in 1:length(n)){
    s <- summary(m[[j]])
    saveRDS(s, file = s_out_rds[j])
    capture.output(s, file = s_out_txt[j])
    if (class(m[[j]]) != "glmerMod") {
      a <- anova(m[[j]])
      saveRDS(a, file = a_out_rds[j])
      capture.output(a, file = a_out_txt[j]) }
    t <- compTukeyCLD(m[[j]])
    saveRDS(t, file = t_out_rds[j])
    capture.output(t, file = t_out_txt[j])
  } 
}




############# output all models estimates as a global list ############# 

getEstSE <- function(x, groups = "Population") {
  require(emmeans)
  df <- as.data.frame(emmeans(x, groups, mode = "asymp"))[,1:3]
  colnames(df)[2] <- "Estimate"
  df
}

combineEstSE <- function(x) {
  for (lab in 1:length(x)) {
    info <- str_split(names(x)[lab], "_") %>% unlist
    if (length(info) >= 5) info[1] <- paste(info[1], info[2], sep = "_")
    info[1] <- sub("_F", "", info[1])
    info[1] <- sub("_M", "", info[1])
    es <- x[[lab]] %>% 
      mutate(Model = paste(info[length(info)-1], info[length(info)], sep = "_"), Trait = info[1], Lab = info[length(info)-2], Sex = info[length(info)-3]) %>% relocate(Model, Trait, Lab, Sex) 
    if (!unique(es$Sex) %in% c("F", "M")) es$Sex <- "NA" 
    if (length(unique(es$Trait)) == 1 & unique(es$Trait) %in% c("Dia", "Fec")) es$Sex <- "F"
    if (length(unique(es$Trait)) == 1 & grepl("Pgm_", unique(es$Trait))) es$Sex <- "F"
    if (length(unique(es$Trait)) == 1 & grepl("LA_", unique(es$Trait))) es$Sex <- "B" 
    x[[lab]] <- es }
  bind_rows(x) }





# list models outputs to keep the directory structure
lmers <- list.files(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "_lmers_pop.rds")
glmers <- list.files(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "_glmers_pop.rds")
models <- c(lmers, glmers)

all_model_estimates <- list()
for (trait in 1:length(models)){
  f <- models[trait]
  n <- str_match(f, '([^/]+)(?:/[^/]+){0}$')[,1]
  n <- sub("_lmers_pop.rds", "_lmer", n)
  n <- sub("_glmers_pop.rds", "_glmer", n)
  n <- tolower(n)
  m <- readRDS(f)
  e <- lapply(m, getEstSE)
  all_model_estimates[[trait]] <- combineEstSE(e)
  names(all_model_estimates)[trait] <- n
}
saveRDS(all_model_estimates, file = "LinearModelsPop/all_models_list_pop_estimates.rds")

write.csv(bind_rows(all_model_estimates), file = "LinearModelsPop/all_models_pop_estimates.csv", row.names = F)



############# output all models estimates by trait ############# 

lmers <- list.files(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "_lmers_pop.rds")
glmers <- list.files(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "_glmers_pop.rds")
models <- c(lmers, glmers)

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


############ get all lmers pvalues

all_lmers_pop_anova <- readRDS("LinearModelsPop/all_lmers_list_pop_anova.rds")
all_glmers_pop_anova <- readRDS("LinearModelsPop/all_glmers_list_pop_anova.rds")
all_pop_anova <- c(all_lmers_pop_anova, all_glmers_pop_anova)


combinePValues <- function(x) {
  pvals <- list()
  for (a in 1:length(x)) {
    info <- str_split(names(x)[a], "_") %>% unlist
    if (length(info) >= 5) info[1] <- paste(info[1], info[2], sep = "_")
    info[1] <- sub("_F", "", info[1])
    info[1] <- sub("_M", "", info[1])
    info <- data.frame(Lab = info[length(info)-2], 
                       Trait = info[1], 
                       Sex = info[length(info)-3],
                       Model = paste(info[length(info)-1], info[length(info)], sep = "_"))
    info$P <- ifelse(grepl("glmer", info$Model), x[[a]]$P[2], x[[a]]$P[1])
    if (!info$Sex %in% c("F", "M")) info$Sex <- "NA" 
    if (info$Trait %in% c("Dia", "Fec")) info$Sex <- "F"
    if (grepl("Pgm_", info$Trait)) info$Sex <- "F"
    if (grepl("LA_", info$Trait)) info$Sex <- "B" 
    pvals[[a]] <- info }
  bind_rows(pvals) }
    

pop_pvalues <- combinePValues(all_pop_anova)

saveRDS(pop_pvalues, file = "LinearModelsPop/all_models_pop_pvalues.rds")
write.csv(pop_pvalues, "LinearModelsPop/all_models_pop_pvalues.csv", row.names = F)










