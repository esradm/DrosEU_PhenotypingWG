

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


##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")


##### source functions
source("Functions/lab_correlations_functions.R")
source("Functions/meta_analysis_functions.R")

##### load data
droseu <- readRDS("Data/droseu_master_list_2022-03-25.rds")

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
DT_lmers_pop$DT_P_Schmidt_lmer_pop <- lmer(DT_EggPupa ~ Population + (1|Population:Line), data = filter(droseu$dtp, Supervisor.PI == "Schmidt"))

# save output list
#saveRDS(DT_P_lmers_pop, file = file.path(lmer_dir, out_dir, "DT_P_lmers_pop.rds"))


### Egg-to-Adult

# initialize output list
#DT_A_lmers_pop <- list()

## Females

# Gibert
DT_lmers_pop$DT_A_F_Gibert_lmer_pop <- lmer(DT_EggAdult ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(droseu$dta, Supervisor.PI == "Gibert" & Sex == "F"))

# Grath, Batch is removed. Singular fit: Line explains nothing, removed it
DT_lmers_pop$DT_A_F_Grath_lmer_pop <- lmer(DT_EggAdult ~ Population + (1|ReplicateVial:Line), data = filter(droseu$dta, Supervisor.PI == "Grath" & Sex == "F"))

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

# Grath, Batch is removed. Singular fit: Line explains nothing, removed it
DT_lmers_pop$DT_A_M_Grath_lmer_pop <- lmer(DT_EggAdult ~ Population + (1|ReplicateVial:Line), data = filter(droseu$dta, Supervisor.PI == "Grath" & Sex == "M"))

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
DW_lmers_pop$DW_F_Colinet_lmer_pop <- lmer(DW_micrograms ~ Population + (1|Population:Line) + (1|Batch), data = filter(droseu$dw, Supervisor.PI == "Colinet" & Sex == "F"))

# Hoedjes
DW_lmers_pop$DW_F_Hoedjes_lmer_pop <- lmer(DW_micrograms ~ Population + (1|Population:Line) + (1|Batch), data = filter(droseu$dw, Supervisor.PI == "Hoedjes" & Sex == "F"))

# Onder
DW_lmers_pop$DW_F_Onder_lmer_pop <- lmer(DW_micrograms ~ Population + (1|Population:Line) + (1|Batch), data = filter(droseu$dw, Supervisor.PI == "Onder" & Sex == "F"))

## Males

# Colinet, singular fit, removed batch
DW_lmers_pop$DW_M_Colinet_lmer_pop <- lmer(DW_micrograms ~ Population + (1|Population:Line), data = filter(droseu$dw, Supervisor.PI == "Colinet" & Sex == "M"))

# Hoedjes
DW_lmers_pop$DW_M_Hoedjes_lmer_pop <- lmer(DW_micrograms ~ Population + (1|Population:Line) + (1|Batch), data = filter(droseu$dw, Supervisor.PI == "Hoedjes" & Sex == "M"))

# Onder
DW_lmers_pop$DW_M_Onder_lmer_pop <- lmer(DW_micrograms ~ Population + (1|Population:Line) + (1|Batch), data = filter(droseu$dw, Supervisor.PI == "Onder" & Sex == "M"))

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

# Posnien, Batch is removed
TL_lmers_pop$TL_F_Posnien_lmer_pop <- lmer(TL_micrometers ~ Population + (1|Line:Population) + (1|ReplicateVial:Line), data = filter(droseu$tl, Supervisor.PI == 'Posnien' & Sex == "F"))

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

# Posnien, does not converge, removed Rep, no more warnings, same output as when "does not converge". Also, no Batch
WA_lmers_pop$WA_L_F_Posnien_lmer_pop <- lmer(CentroidSizeLeft_micrometers ~ Population + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Posnien" & Sex == "F"))

# Ritchie
WA_lmers_pop$WA_L_F_Ritchie_lmer_pop <- lmer(CentroidSizeLeft_micrometers ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line:Population), data = filter(droseu$wa, Supervisor.PI == "Ritchie" & Sex == "F"))

# StamenkovicRadak, singular fit warning, removed Batch. Nearly unidentifiable model, removed Rep, no more warnings, same output as when "nearly unidentifiable".
WA_lmers_pop$WA_L_F_StamenkovicRadak_lmer_pop <- lmer(CentroidSizeLeft_micrometers ~ Population + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "StamenkovicRadak" & Sex == "F"))


## Males left

# Onder
WA_lmers_pop$WA_L_M_Onder_lmer_pop <- lmer(CentroidSizeLeft_micrometers ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line:Population), data = filter(droseu$wa, Supervisor.PI == "Onder" & Sex == "M"))

# Posnien, warning nearly unidentifiable model, removed Rep, no more warnings, same output as when "nearly unidentifiable". Also, no Batch
WA_lmers_pop$WA_L_M_Posnien_lmer_pop <- lmer(CentroidSizeLeft_micrometers ~ Population + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Posnien" & Sex == "M"))

# Ritchie
WA_lmers_pop$WA_L_M_Ritchie_lmer_pop <- lmer(CentroidSizeLeft_micrometers ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line:Population), data = filter(droseu$wa, Supervisor.PI == "Ritchie" & Sex == "M"))

# StamenkovicRadak, warning nearly unidentifiable model, removed Rep, no more warnings, same output as when "nearly unidentifiable"
WA_lmers_pop$WA_L_M_StamenkovicRadak_lmer_pop <- lmer(CentroidSizeLeft_micrometers ~ Population + (1|Line:Population) + (1|Batch), data = filter(droseu$wa, Supervisor.PI == "StamenkovicRadak" & Sex == "M"))


## Females right

# Onder
WA_lmers_pop$WA_R_F_Onder_lmer_pop <- lmer(CentroidSizeRight_micrometers ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(droseu$wa, Supervisor.PI == "Onder" & Sex == "F"))

# Posnien, does not converge, removed Rep, no more warnings, same output as when "does not converge". Also, no Batch
WA_lmers_pop$WA_R_F_Posnien_lmer_pop <- lmer(CentroidSizeRight_micrometers ~ Population + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Posnien" & Sex == "F"))

# Ritchie
WA_lmers_pop$WA_R_F_Ritchie_lmer_pop <- lmer(CentroidSizeRight_micrometers ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line), data = filter(droseu$wa, Supervisor.PI == "Ritchie" & Sex == "F"))

# StamenkovicRadak, singular fit warning, removed Batch. Failed to converge, removed Rep, no more warnings, same output as when "does not converge"
WA_lmers_pop$WA_R_F_StamenkovicRadak_lmer_pop <- lmer(CentroidSizeRight_micrometers ~ Population + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "StamenkovicRadak" & Sex == "F"))


## Males right

# Onder, failed to converge, removing Batch or Rep fixes it. removed batch as it is the variable that explains the least
WA_lmers_pop$WA_R_M_Onder_lmer_pop <- lmer(CentroidSizeRight_micrometers ~ Population + (1|Line:Population) + (1|ReplicateVial:Line), data = filter(droseu$wa, Supervisor.PI == "Onder" & Sex == "M"))

# Posnien, does not converge, removed Rep, no more warnings, same output as when "does not converge". Also, no Batch
WA_lmers_pop$WA_R_M_Posnien_lmer_pop <- lmer(CentroidSizeRight_micrometers ~ Population + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Posnien" & Sex == "M"))

# Ritchie
WA_lmers_pop$WA_R_M_Ritchie_lmer_pop <- lmer(CentroidSizeRight_micrometers ~ Population + (1|Line:Population) + (1|Batch) + (1|ReplicateVial:Line:Population), data = filter(droseu$wa, Supervisor.PI == "Ritchie" & Sex == "M"))

# StamenkovicRadak, singular fit, removed Batch. Does not converge, removed Rep, no more warnings, same output as when "does not converge"
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
LS_lmers_pop$LS_P_F_Flatt_lmer_pop <- lmer(LSP_AgeAtDeath_days ~ Population + (1|Population:ReplicateCage), data = filter(droseu$lsp, Censor == "0" & Supervisor.PI == "Flatt" & Sex == "F"))

# Parsch
LS_lmers_pop$LS_L_F_Parsch_lmer_pop <- lmer(LSL_AgeAtDeath_days ~ Population + (1|Batch) + (1|Population:Line) + (1|Line:ReplicateVial), data = filter(droseu$lsl, Censor == "0" & Supervisor.PI == "Parsch" & Sex == "F"))

# Pasyukova
LS_lmers_pop$LS_L_F_Pasyukova_lmer_pop <- lmer(LSL_AgeAtDeath_days ~ Population + (1|Batch) + (1|Population:Line) + (1|Line:ReplicateVial), data = filter(droseu$lsl, Censor == "0" & Supervisor.PI == "Pasyukova" & Sex == "F"))


## Males

# Flatt
LS_lmers_pop$LS_P_M_Flatt_lmer_pop <- lmer(LSP_AgeAtDeath_days ~ Population + (1|Population:ReplicateCage), data = filter(droseu$lsp, Censor == "0" & Supervisor.PI == "Flatt" & Sex == "M"))

# Parsch
LS_lmers_pop$LS_L_M_Parsch_lmer_pop <- lmer(LSL_AgeAtDeath_days ~ Population + (1|Batch) + (1|Population:Line) + (1|Line:ReplicateVial), data = filter(droseu$lsl, Censor == "0" & Supervisor.PI == "Parsch" & Sex == "M"))

# Pasyukova, failed to coverge, removed Batch as it explains the least, same output as when does not converge
LS_lmers_pop$LS_L_M_Pasyukova_lmer_pop <- lmer(LSL_AgeAtDeath_days ~ Population + (1|Population:Line) + (1|Line:ReplicateVial), data = filter(droseu$lsl, Censor == "0" & Supervisor.PI == "Pasyukova" & Sex == "M"))


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
CSM_lmers_pop$CSM_F_Gonzalez_lmer_pop <- lmer(CSM_PropDead_ED_asin ~ Population + (1|Line: Population) + (1|Batch), data = filter(droseu$csm, Supervisor.PI == "Gonzalez" & Sex == "F"))

# Kozeretska
CSM_lmers_pop$CSM_F_Kozeretska_lmer_pop <- lmer(CSM_PropDead_ED_asin ~ Population + (1|Line: Population) + (1|Batch), data = filter(droseu$csm, Supervisor.PI == "Kozeretska" & Sex == "F"))

# Vieira
CSM_lmers_pop$CSM_F_Vieira_lmer_pop <- lmer(CSM_PropDead_ED_asin ~ Population + (1|Line: Population) + (1|Batch), data = filter(droseu$csm, Supervisor.PI == "Vieira" & Sex == "F"))


## Males

# Gonzalez
CSM_lmers_pop$CSM_M_Gonzalez_lmer_pop <- lmer(CSM_PropDead_ED_asin ~ Population + (1|Line: Population) + (1|Batch), data = filter(droseu$csm, Supervisor.PI == "Gonzalez" & Sex == "M"))

# Kozeretska
CSM_lmers_pop$CSM_M_Kozeretska_lmer_pop <- lmer(CSM_PropDead_ED_asin ~ Population + (1|Line: Population) + (1|Batch), data = filter(droseu$csm, Supervisor.PI == "Kozeretska" & Sex == "M"))

# Vieira
CSM_lmers_pop$CSM_M_Vieira_lmer_pop <- lmer(CSM_PropDead_ED_asin ~ Population + (1|Line: Population) + (1|Batch), data = filter(droseu$csm, Supervisor.PI == "Vieira" & Sex == "M"))


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
CCRT_lmers_pop$CCRT_F_Vieira_lmer_pop <- lmer(CCRT_seconds ~ Population + (1|Population:Line) + (1|Batch) + (1|ReplicateVial:Line), data = filter(droseu$ccrt, Censor == "0" & Supervisor.PI == "Vieira" & Sex == "F"))

# Mensh
#CCRT_lmers_pop$CCRT_F_Mensh_lmer_pop <- lmer(CCRT_seconds ~ Population + (1|Population:Line) + (1|Batch) + (1|ReplicateVial:Line), data = filter(droseu$ccrt, Censor == "0" & Supervisor.PI == "Mensh" & Sex == "F"))



## Males

# Vieira, singular fit, removed Batch
CCRT_lmers_pop$CCRT_M_Vieira_lmer_pop <- lmer(CCRT_seconds ~ Population + (1|Population:Line) + (1|ReplicateVial:Line), data = filter(droseu$ccrt, Censor == "0" & Supervisor.PI == "Vieira" & Sex == "M"))

# Mensh
#CCRT_lmers_pop$CCRT_M_Mensh_lmer_pop <- lmer(CCRT_seconds ~ Population + (1|Population:Line) + (1|Batch) + (1|ReplicateVial:Line), data = filter(droseu$ccrt, Censor == "0" & Supervisor.PI == "Mensh" & Sex == "M"))

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

# Bergland, singular fit, removed Batch
Dia_lmers_pop$Dia_Bergland_lmer_pop <- lmer(Prop_Max_Stage9_asin ~ Population + (1|Population:Line), data = filter(droseu$dia, Supervisor.PI == "Bergland"))

# Flatt
Dia_lmers_pop$Dia_Flatt_lm_pop <- lm(Prop_Max_Stage9_asin ~ Population, data = filter(droseu$dia, Supervisor.PI == "Flatt"))

# Schlotterer, singular fit, removed Batch
Dia_lmers_pop$Dia_Schlotterer_lmer_pop <- lmer(Prop_Max_Stage9_asin ~ Population + (1|Population:Line) + (1|Batch), data = filter(droseu$dia, Supervisor.PI == "Schlotterer"))

# save output list
saveRDS(Dia_lmers_pop, file = file.path(lmer_dir, out_dir, "Dia_lmers_pop.rds"))




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
LA_lmers_pop$LA_NDlog2_Tauber_lmer_pop <- lmer(ND_log2 ~ Population + (1|Line:Population), data = filter(droseu$la, ND_log2 != -Inf)) # log2 transformation creates -Inf values (2 cases), should they be removed?

# untransformed, not converging, removed Batch that explains the least
LA_lmers_pop$LA_ND_Tauber_lmer_pop <- lmer(ND ~ Population + (1|Line:Population), data = droseu$la)

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
Pgm_lmers_pop$Pgm_T8_Schmidt_lmer_pop <- lmer(PercT8 ~ Population + (1|Line:Population), data = filter(droseu$pgm2, Supervisor.PI == "Schmidt"))

Pgm_lmers_pop$Pgm_T9_Schmidt_lmer_pop <- lmer(PercT9 ~ Population + (1|Line:Population), data = filter(droseu$pgm2, Supervisor.PI == "Schmidt"))

Pgm_lmers_pop$Pgm_T10_Schmidt_lmer_pop <- lmer(PercT10 ~ Population + (1|Line:Population), data = filter(droseu$pgm2, Supervisor.PI == "Schmidt"))

Pgm_lmers_pop$Pgm_Total_Schmidt_lmer_pop <- lmer(TotalPerc ~ Population + (1|Line:Population), data = filter(droseu$pgm2, Supervisor.PI == "Schmidt"))


# save output list
saveRDS(Pgm_lmers_pop, file = file.path(lmer_dir, out_dir, "Pgm_lmers_pop.rds"))





############# ouput all lmers ############# 

all_lmers_pop <- ls()[grep("lmers", ls())]
saveRDS(all_lmers_pop, file = "LinearModelsPop/all_lmers_pop.rds")





############# check linear models residuals ############# 



lmers <- list.files(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "lmers_pop.rds")[-1]


for (i in 1:length(lmers)){
  f <- lmers[i]
  p <- str_match(f, '(.*[^/]+)(?:/[^/]+){1}$')[,2]
  dir.create(file.path(p, "by_lab_lmer_residuals"), showWarnings = F) 
  m <- readRDS(f)
  n <- names(m)
  f_out <- file.path(p, "by_lab_lmer_residuals", n)
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
 


############# ouput all lmers summaries and anovas ############# 





compTukeyCLD <- function(x) {
  em <- emmeans(x, pairwise ~ Population, mode = "asymp", adjust = "tukey")
  let <- cld(em, Letters = letters, alpha = 0.05) %>% 
    as.data.frame %>% mutate(Population = as.factor(Population)) %>% 
    arrange(Population) %>% dplyr::rename(cld = .group)
  list(contrasts = em, letters = let)
}





lmers <- list.files(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "lmers_pop.rds")[-1]




for (i in 1:length(lmers)){
  f <- lmers[i]
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
  a <- lapply(m, anova)
  saveRDS(a, file = a_out_rds)
  capture.output(a, file = a_out_txt)
  t <- lapply(m, compTukeyCLD)
  saveRDS(t, file = t_out_rds)
  capture.output(t, file = t_out_txt)
}
  

##### by trait and by lab

for (i in 1:length(lmers)){
  f <- lmers[i]
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
    a <- anova(m[[j]])
    saveRDS(a, file = a_out_rds[j])
    capture.output(a, file = a_out_txt[j])
    t <- compTukeyCLD(m[[j]])
    saveRDS(t, file = t_out_rds[j])
    capture.output(t, file = t_out_txt[j])
  } 
}


############# ouput all models estimates ############# 


all_model_estimates <- list()
for (i in 1:length(lmers)){
  f <- lmers[i]
  n <- str_match(f, '([^/]+)(?:/[^/]+){0}$')[,1]
  n <- sub("_lmers_pop.rds", "", n)
  n <- tolower(n)
  m <- readRDS(f)
  e <- lapply(m, getEstSE)
  all_model_estimates[[i]] <- combineEst3(e)
  names(all_model_estimates)[i] <- n
}
saveRDS(all_model_estimates, file = "LinearModelsPop/all_model_estimates.rds")


##### by trait

for (i in 1:length(lmers)){
  f <- lmers[i]
  m <- readRDS(f)
  e <- lapply(m, getEstSE)
  e <- combineEst3(e)
  e_out_rds <- sub(".rds", "_model_estimates.rds", f)
  e_out_txt <- sub(".rds", "_model_estimates.txt", f)
  saveRDS(e, file = e_out_rds)
  write.table(e, file = e_out_txt, row.names = F, quote = F, sep = "\t")
}










