


########################################################################
######################  RUN ALL LMMs FOR LATITUDE ######################
########################################################################




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
#source("Functions/lab_correlations_functions.R")

##### load data
droseu <- readRDS("Data/droseu_master_list_2022-04-05.rds")

##### create output directory
lmer_dir <- "LinearModelsLat"
dir.create(lmer_dir, showWarnings = F) 





############# VIABILITY #############

# create output directory
out_dir <- "Viability"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
Via_lmers_lat <- list()

#### Gibert Lab
Via_lmers_lat$Via_Gibert_lmer_lat <- lmer(ProportionEggtoAdultSurvival ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$via, Supervisor.PI == "Gibert"))

#### Grath Lab
Via_lmers_lat$Via_Grath_lmer_lat <- lmer(ProportionEggtoAdultSurvival ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$via, Supervisor.PI == "Grath"))

#### Hoedjes Lab
Via_lmers_lat$Via_Hoedjes_lmer_lat <- lmer(ProportionEggtoAdultSurvival ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$via, Supervisor.PI == "Hoedjes"))

#### Schmidt Lab
Via_lmers_lat$Via_Schmidt_lmer_lat <- lmer(ProportionEggtoAdultSurvival ~ Latitude + (1|Population), data = filter(droseu$via, Supervisor.PI == "Schmidt"))

#### Stamenkovic-Radak Lab
Via_lmers_lat$Via_StamenkovicRadak_lmer_lat <- lmer(ProportionEggtoAdultSurvival ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$via, Supervisor.PI == "StamenkovicRadak"))

#### Zwaan Lab
Via_lmers_lat$Via_Zwaan_lmer_lat <- lmer(ProportionEggtoAdultSurvival ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$via, Supervisor.PI == "Zwaan"))

# save output list
saveRDS(Via_lmers_lat, file = file.path(lmer_dir, out_dir, "Via_lmers_lat.rds"))








############# DEVELOPMENT TIME #############

# create output directory
out_dir <- "DevelopmentTime"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 


### Egg-to-Pupae

# initialize output list
DT_lmers_lat <- list()

##### Schmidt Lab
DT_lmers_lat$DT_P_Schmidt_lmer_lat <- lmer(DT_EggPupa ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$dtp, Supervisor.PI == "Schmidt"))

### Egg-to-adult developmental time

### Females

##### Gibert Lab
DT_lmers_lat$DT_A_F_Gibert_lmer_lat <- lmer(DT_EggAdult ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "Gibert" & Sex == "F"))

##### Grath Lab
DT_lmers_lat$DT_A_F_Grath_lmer_lat <- lmer(DT_EggAdult ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "Grath" & Sex == "F"))

##### Hoedjes Lab
DT_lmers_lat$DT_A_F_Hoedjes_lmer_lat <- lmer(DT_EggAdult ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "Hoedjes" & Sex == "F"))

##### Schmidt Lab
DT_lmers_lat$DT_A_F_Schmidt_lmer_lat <- lmer(DT_EggAdult ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "Schmidt" & Sex == "F"))

##### Stamenkovic-Radak Lab
DT_lmers_lat$DT_A_F_StamenkovicRadak_lmer_lat <- lmer(DT_EggAdult ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "StamenkovicRadak" & Sex == "F"))

##### Zwaan Lab
DT_lmers_lat$DT_A_F_Zwaan_lmer_lat <- lmer(DT_EggAdult ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "Zwaan" & Sex == "F"))


### Males

##### Gibert Lab
DT_lmers_lat$DT_A_M_Gibert_lmer_lat <- lmer(DT_EggAdult ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "Gibert" & Sex == "M"))

##### Grath Lab
DT_lmers_lat$DT_A_M_Grath_lmer_lat <- lmer(DT_EggAdult ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "Grath" & Sex == "M"))

##### Hoedjes Lab
DT_lmers_lat$DT_A_M_Hoedjes_lmer_lat <- lmer(DT_EggAdult ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "Hoedjes" & Sex == "M"))

##### Schmidt Lab
DT_lmers_lat$DT_A_M_Schmidt_lmer_lat <- lmer(DT_EggAdult ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "Schmidt" & Sex == "M"))

##### Stamenkovic-Radak Lab
DT_lmers_lat$DT_A_M_StamenkovicRadak_lmer_lat <- lmer(DT_EggAdult ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "StamenkovicRadak" & Sex == "M"))

##### Zwaan Lab
DT_lmers_lat$DT_A_M_Zwaan_lmer_lat <- lmer(DT_EggAdult ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "Zwaan" & Sex == "M"))

# save output list
saveRDS(DT_lmers_lat, file = file.path(lmer_dir, out_dir, "DT_lmers_lat.rds"))









############# DRY WEIGHT #############

# create output directory
out_dir <- "DryWeight"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
DW_lmers_lat <- list()

## Females

##### Colinet Lab
DW_lmers_lat$DW_F_Colinet_lmer_lat <- lmer(DW_micrograms ~ Latitude + (1|Population) + (1|Population:Line), data = filter(droseu$dw, Supervisor.PI == "Colinet" & Sex == "F"))

##### Hoedjes Lab, not converging, removed Population
DW_lmers_lat$DW_F_Hoedjes_lmer_lat <- lmer(DW_micrograms ~ Latitude + (1|Population:Line), data = filter(droseu$dw, Supervisor.PI == "Hoedjes" & Sex == "F"))

##### Onder Lab
DW_lmers_lat$DW_F_Onder_lmer_lat <- lmer(DW_micrograms ~ Latitude + (1|Population) + (1|Population:Line), data = filter(droseu$dw, Supervisor.PI == "Onder" & Sex == "F"))


## Males

##### Colinet Lab
DW_lmers_lat$DW_M_Colinet_lmer_lat <- lmer(DW_micrograms ~ Latitude + (1|Population) + (1|Population:Line), data = filter(droseu$dw, Supervisor.PI == "Colinet" & Sex == "M"))

##### Hoedjes Lab, singular fit removed Population
DW_lmers_lat$DW_M_Hoedjes_lmer_lat <- lmer(DW_micrograms ~ Latitude + (1|Population:Line), data = filter(droseu$dw, Supervisor.PI == "Hoedjes" & Sex == "M"))

##### Onder Lab
DW_lmers_lat$DW_M_Onder_lmer_lat <- lmer(DW_micrograms ~ Latitude + (1|Population) + (1|Population:Line), data = filter(droseu$dw, Supervisor.PI == "Onder" & Sex == "M"))

# save output list
saveRDS(DW_lmers_lat, file = file.path(lmer_dir, out_dir, "DW_lmers_lat.rds"))











############# THORAX LENGTH #############

# create output directory
out_dir <- "ThoraxLength"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
TL_lmers_lat <- list()

## Females

##### Kozeretska Lab
TL_lmers_lat$TL_F_Kozeretska_lmer_lat <- lmer(TL_micrometers ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$tl, Supervisor.PI == 'Kozeretska' & Sex == "F"))

##### Posnien Lab
TL_lmers_lat$TL_F_Posnien_lmer_lat <- lmer(TL_micrometers ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$tl, Supervisor.PI == 'Posnien' & Sex == "F"))

##### Ritchie Lab, singular fit, removed Population
TL_lmers_lat$TL_F_Ritchie_lmer_lat <- lmer(TL_micrometers ~ Latitude + (1|Line:Population), data = filter(droseu$tl, Supervisor.PI == 'Ritchie' & Sex == "F"))

##### Schmidt Lab
TL_lmers_lat$TL_F_Schmidt_lmer_lat <- lmer(TL_micrometers ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$tl, Supervisor.PI == 'Schmidt' & Sex == "F"))


## Males

##### Kozeretska Lab
TL_lmers_lat$TL_M_Kozeretska_lmer_lat <- lmer(TL_micrometers ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$tl, Supervisor.PI == 'Kozeretska' & Sex == "M"))

##### Posnien Lab
TL_lmers_lat$TL_M_Posnien_lmer_lat <- lmer(TL_micrometers ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$tl, Supervisor.PI == 'Posnien' & Sex == "M"))

##### Ritchie Lab, singular fit, removed Population
TL_lmers_lat$TL_M_Ritchie_lmer_lat <- lmer(TL_micrometers ~ Latitude + (1|Line:Population), data = filter(droseu$tl, Supervisor.PI == 'Ritchie' & Sex == "M"))

# save output list
saveRDS(TL_lmers_lat, file = file.path(lmer_dir, out_dir, "TL_lmers_lat.rds"))







############# WING AREA ############# 

# create output directory
out_dir <- "WingArea"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
WA_lmers_lat <- list()


## Females left

##### Onder Lab
WA_lmers_lat$WA_L_F_Onder_lmer_lat <- lmer(CentroidSizeLeft_micrometers ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Onder" & Sex == "F"))

##### Posnien Lab, singular fit, removed Population
WA_lmers_lat$WA_L_F_Posnien_lmer_lat <- lmer(CentroidSizeLeft_micrometers ~ Latitude + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Posnien" & Sex == "F"))

##### Ritchie Lab
WA_lmers_lat$WA_L_F_Ritchie_lmer_lat <- lmer(CentroidSizeLeft_micrometers ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Ritchie" & Sex == "F"))

##### Onder Lab
WA_lmers_lat$WA_L_F_StamenkovicRadak_lmer_lat <- lmer(CentroidSizeLeft_micrometers ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "StamenkovicRadak" & Sex == "F"))


## Males left

##### Onder Lab
WA_lmers_lat$WA_L_M_Onder_lmer_lat <- lmer(CentroidSizeLeft_micrometers ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Onder" & Sex == "M"))

##### Posnien Lab, singular fit, removed Population
WA_lmers_lat$WA_L_M_Posnien_lmer_lat <- lmer(CentroidSizeLeft_micrometers ~ Latitude + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Posnien" & Sex == "M"))

##### Ritchie Lab
WA_lmers_lat$WA_L_M_Ritchie_lmer_lat <- lmer(CentroidSizeLeft_micrometers ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Ritchie" & Sex == "M"))

##### Onder Lab
WA_lmers_lat$WA_L_M_StamenkovicRadak_lmer_lat <- lmer(CentroidSizeLeft_micrometers ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "StamenkovicRadak" & Sex == "M"))



## Females right

##### Onder Lab
WA_lmers_lat$WA_R_F_Onder_lmer_lat <- lmer(CentroidSizeRight_micrometers ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Onder" & Sex == "F"))

##### Posnien Lab, singular fit, removed Population
WA_lmers_lat$WA_R_F_Posnien_lmer_lat <- lmer(CentroidSizeRight_micrometers ~ Latitude + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Posnien" & Sex == "F"))

##### Ritchie Lab
WA_lmers_lat$WA_R_F_Ritchie_lmer_lat <- lmer(CentroidSizeRight_micrometers ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Ritchie" & Sex == "F"))

##### Onder Lab
WA_lmers_lat$WA_R_F_StamenkovicRadak_lmer_lat <- lmer(CentroidSizeRight_micrometers ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "StamenkovicRadak" & Sex == "F"))


## Mmales right

##### Onder Lab
WA_lmers_lat$WA_R_M_Onder_lmer_lat <- lmer(CentroidSizeRight_micrometers ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Onder" & Sex == "M"))

##### Posnien Lab, singular fit, removed Population
WA_lmers_lat$WA_R_M_Posnien_lmer_lat <- lmer(CentroidSizeRight_micrometers ~ Latitude + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Posnien" & Sex == "M"))

##### Ritchie Lab
WA_lmers_lat$WA_R_M_Ritchie_lmer_lat <- lmer(CentroidSizeRight_micrometers ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Ritchie" & Sex == "M"))

##### Onder Lab
WA_lmers_lat$WA_R_M_StamenkovicRadak_lmer_lat <- lmer(CentroidSizeRight_micrometers ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "StamenkovicRadak" & Sex == "M"))

# save output list
saveRDS(WA_lmers_lat, file = file.path(lmer_dir, out_dir, "WA_lmers_lat.rds"))






############# FECUNDITY ############# 

# create output directory
out_dir <- "Fecundity"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
Fec_lmers_lat <- list()

##### Billeter Lab
Fec_lmers_lat$Fec_Billeter_lmer_lat <- lmer(NumberOfAdultsEclosed ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$fec, Supervisor.PI == "Billeter"))

##### Fricke Lab, singular fit, removed Population
Fec_lmers_lat$Fec_Fricke_lmer_lat <- lmer(NumberOfAdultsEclosed ~ Latitude + (1|Line:Population), data = filter(droseu$fec, Supervisor.PI == "Fricke"))

# save output list
saveRDS(Fec_lmers_lat, file = file.path(lmer_dir, out_dir, "Fec_lmers_lat.rds"))







############# LIFESPAN ############# 

# create output directory
out_dir <- "Lifespan"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
LS_lmers_lat <- list()


## Females

#### Flatt Lab
LS_lmers_lat$LS_F_Flatt_lmer_lat <- lmer(LSP_AgeAtDeath_days ~ Latitude + (1|Population), data = filter(droseu$lsp, Censor == "0" & Supervisor.PI == "Flatt" & Sex == "F"))

#### Parsch Lab
LS_lmers_lat$LS_F_Parsch_lmer_lat <- lmer(LSL_AgeAtDeath_days ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$lsl, Censor == "0" & Supervisor.PI == "Parsch" & Sex == "F"))

#### Pasyukova Lab
LS_lmers_lat$LS_F_Pasyukova_lmer_lat <- lmer(LSL_AgeAtDeath_days ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$lsl, Censor == "0" & Supervisor.PI == "Pasyukova" & Sex == "F"))


## Males

#### Flatt Lab
LS_lmers_lat$LS_M_Flatt_lmer_lat <- lmer(LSP_AgeAtDeath_days ~ Latitude + (1|Population), data = filter(droseu$lsp, Censor == "0" & Supervisor.PI == "Flatt" & Sex == "M"))

#### Parsch Lab
LS_lmers_lat$LS_M_Parsch_lmer_lat <- lmer(LSL_AgeAtDeath_days ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$lsl, Censor == "0" & Supervisor.PI == "Parsch" & Sex == "M"))

#### Pasyukova Lab
LS_lmers_lat$LS_M_Pasyukova_lmer_lat <- lmer(LSL_AgeAtDeath_days ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$lsl, Censor == "0" & Supervisor.PI == "Pasyukova" & Sex == "M"))


# save output list
saveRDS(LS_lmers_lat, file = file.path(lmer_dir, out_dir, "LS_lmers_lat.rds"))







############# COLD-SHOCK MORTALITY ############# 

# create output directory
out_dir <- "ColdShock"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
CSM_lmers_lat <- list()


## Females

#### Gonzalez Lab
CSM_lmers_lat$CSM_F_Gonzalez_lmer_lat <- lmer(CSM_PropDead_ED_asin ~ Latitude + (1|Line:Population) + (1|Population), data = filter(droseu$csm, Supervisor.PI == "Gonzalez" & Sex == "F"))

#### Kozeretska Lab, singular fit, removed population
CSM_lmers_lat$CSM_F_Kozeretska_lmer_lat <- lmer(CSM_PropDead_ED_asin ~ Latitude + (1|Line:Population), data = filter(droseu$csm, Supervisor.PI == "Kozeretska" & Sex == "F"))

#### Vieira Lab
CSM_lmers_lat$CSM_F_Vieira_lmer_lat <- lmer(CSM_PropDead_ED_asin ~ Latitude + (1|Line:Population) + (1|Population), data = filter(droseu$csm, Supervisor.PI == "Vieira" & Sex == "F"))


## Males

#### Gonzalez Lab
CSM_lmers_lat$CSM_M_Gonzalez_lmer_lat <- lmer(CSM_PropDead_ED_asin ~ Latitude + (1|Line:Population) + (1|Population), data = filter(droseu$csm, Supervisor.PI == "Gonzalez" & Sex == "M"))

#### Kozeretska Lab, singular fit, removed population
CSM_lmers_lat$CSM_M_Kozeretska_lmer_lat <- lmer(CSM_PropDead_ED_asin ~ Latitude + (1|Line:Population), data = filter(droseu$csm, Supervisor.PI == "Kozeretska" & Sex == "M"))

#### Vieira Lab
CSM_lmers_lat$CSM_M_Vieira_lmer_lat <- lmer(CSM_PropDead_ED_asin ~ Latitude + (1|Line:Population) + (1|Population), data = filter(droseu$csm, Supervisor.PI == "Vieira" & Sex == "M"))

# save output list
saveRDS(CSM_lmers_lat, file = file.path(lmer_dir, out_dir, "CSM_lmers_lat.rds"))






############# CHILL-COMA RECOVERY TIME ############# 

# create output directory
out_dir <- "ChillComa"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
CCRT_lmers_lat <- list()


## Females

#### Vieira Lab
CCRT_lmers_lat$CCRT_F_Vieira_lmer_lat <- lmer(CCRT_seconds ~ Latitude + (1|Population) + (1|Population:Line), data = filter(droseu$ccrt, Censor == "0" & Supervisor.PI == "Vieira" & Sex == "F"))

#### Mensh Lab
#CCRT_lmers_lat$CCRT_F_Mensh_lmer_lat <- lmer(CCRT_seconds ~ Latitude + (1|Population) + (1|Population:Line), data = filter(droseu$ccrt, Censor == "0" & Supervisor.PI == "Mensh" & Sex == "F"))


## Males

#### Vieira Lab
CCRT_lmers_lat$CCRT_M_Vieira_lmer_lat <- lmer(CCRT_seconds ~ Latitude + (1|Population) + (1|Population:Line), data = filter(droseu$ccrt, Censor == "0" & Supervisor.PI == "Vieira" & Sex == "M"))

#### Mensh Lab
#CCRT_lmers_lat$CCRT_M_Mensh_lmer_lat <- lmer(CCRT_seconds ~ Latitude + (1|Population) + (1|Population:Line), data = filter(droseu$ccrt, Censor == "0" & Supervisor.PI == "Mensh" & Sex == "M"))

# save output list
saveRDS(CCRT_lmers_lat, file = file.path(lmer_dir, out_dir, "CCRT_lmers_lat.rds"))






############# HEAT-SHOCK MORTALITY ############# 

# create output directory
out_dir <- "HeatShock"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
HSM_lmers_lat <- list()

## Females

#### Parsch Lab
HSM_lmers_lat$HSM_F_Parsch_lmer_lat <- lmer(TimeDeath_min ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$hsm, Censor == "0" & Supervisor.PI == "Parsch" & Sex == "F"))

#### Vieira Lab
HSM_lmers_lat$HSM_F_Vieira_lmer_lat <- lmer(TimeDeath_min ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$hsm, Censor == "0" & Supervisor.PI == "Vieira" & Sex == "F"))


## Males

#### Parsch Lab
HSM_lmers_lat$HSM_M_Parsch_lmer_lat <- lmer(TimeDeath_min ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$hsm, Censor == "0" & Supervisor.PI == "Parsch" & Sex == "M"))

#### Vieira Lab
HSM_lmers_lat$HSM_M_Vieira_lmer_lat <- lmer(TimeDeath_min ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$hsm, Censor == "0" & Supervisor.PI == "Vieira" & Sex == "M"))

# save output list
saveRDS(HSM_lmers_lat, file = file.path(lmer_dir, out_dir, "HSM_lmers_lat.rds"))






############# DIAPAUSE ############# 

# create output directory
out_dir <- "Diapause"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
Dia_lmers_lat <- list()

#### Bergland Lab
Dia_lmers_lat$Dia_Bergland_lmer_lat <- lmer(Prop_Max_Stage9_asin ~ Latitude + (1|Population) + (1|Population:Line), data = filter(droseu$dia, Supervisor.PI == "Bergland"))

#### Flatt Lab, removed Population because of singular fit
Dia_lmers_lat$Dia_Flatt_lm_lat <- lm(Prop_Max_Stage9_asin ~ Latitude, data = filter(droseu$dia, Supervisor.PI == "Flatt"))

#### Schlotterer Lab
Dia_lmers_lat$Dia_Schlotterer_lmer_lat <- lmer(Prop_Max_Stage9_asin ~ Latitude + (1|Population) + (1|Population:Line), data = filter(droseu$dia, Supervisor.PI == "Schlotterer"))

# save output list
saveRDS(Dia_lmers_lat, file = file.path(lmer_dir, out_dir, "Dia_lmers_lat.rds"))








############# CIRCADIAN ECLOSION TIMING ############# 

# create output directory
out_dir <- "CircadianEclosion"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
CET_lmers_lat <- list()




############# LOCOMOTOR ACTIVITY ############# 

# create output directory
out_dir <- "Locomotor"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
LA_lmers_lat <- list()

#### Tauber Lab, removed "-Inf"
LA_lmers_lat$LA_NDlog2_Tauber_lmer_lat <- lmer(ND_log2 ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$la, ND_log2 != -Inf)) # log2 transformation creates -Inf values (2 cases), should they be removed?

#### Tauber Lab, untransformed ND
#LA_lmers_lat$LA_ND_Tauber_lmer_lat <- lmer(ND ~ Latitude + (1|Population) + (1|Line:Population), data = droseu$la)

LA_lmers_lat$LA_Period_Tauber_lmer_lat <- lmer(Period ~ Latitude + (1|Population) + (1|Line:Population), data = droseu$la)

LA_lmers_lat$LA_CircPhase_Tauber_lmer_lat <- lmer(CircPhase ~ Latitude + (1|Population) + (1|Line:Population), data = droseu$la)

# singular fit, removed Line
LA_lmers_lat$LA_AbsPhase_Tauber_lmer_lat <- lmer(AbsPhase ~ Latitude + (1|Population), data = droseu$la)

LA_lmers_lat$LA_Activity_Tauber_lmer_lat <- lmer(Activity ~ Latitude + (1|Population) + (1|Line:Population), data = droseu$la)



# save output list
saveRDS(LA_lmers_lat, file = file.path(lmer_dir, out_dir, "LA_lmers_lat.rds"))







############# STARVATION RESISTANCE ############# 

# create output directory
out_dir <- "Starvation"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
SR_lmers_lat <- list()

## Females

#### Gonzalez Lab
SR_lmers_lat$SR_F_Gonzalez_lmer_lat <- lmer(AgeAtDeath_hours ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$sr, Supervisor.PI == "Gonzalez" & Sex == "F"))

#### Onder Lab
SR_lmers_lat$SR_F_Onder_lmer_lat <- lmer(AgeAtDeath_hours ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$sr, Supervisor.PI == "Onder" & Sex == "F"))

#### Pasyukova Lab
SR_lmers_lat$SR_F_Pasyukova_lmer_lat <- lmer(AgeAtDeath_hours ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$sr, Supervisor.PI == "Pasyukova" & Sex == "F"))

## Males

#### Gonzalez Lab
SR_lmers_lat$SR_M_Gonzalez_lmer_lat <- lmer(AgeAtDeath_hours ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$sr, Supervisor.PI == "Gonzalez" & Sex == "M"))

#### Onder Lab
SR_lmers_lat$SR_M_Onder_lmer_lat <- lmer(AgeAtDeath_hours ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$sr, Supervisor.PI == "Onder" & Sex == "M"))

#### Pasyukova Lab
SR_lmers_lat$SR_M_Pasyukova_lmer_lat <- lmer(AgeAtDeath_hours ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$sr, Supervisor.PI == "Pasyukova" & Sex == "M"))

# save output list
saveRDS(SR_lmers_lat, file = file.path(lmer_dir, out_dir, "SR_lmers_lat.rds"))






############# PIGMENTATION ############# 

# create output directory
out_dir <- "Pigmentation"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
Pgm_lmers_lat <- list()

#### Abbott Lab
Pgm_lmers_lat$Pgm_T4_Abbott_lmer_lat <- lmer(PercT4_asin ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Abbott"))

Pgm_lmers_lat$Pgm_T5_Abbott_lmer_lat <- lmer(PercT5_asin ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Abbott"))

Pgm_lmers_lat$Pgm_T6_Abbott_lmer_lat <- lmer(PercT6_asin ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Abbott"))

Pgm_lmers_lat$Pgm_Total_Abbott_lmer_lat <- lmer(TotalPerc_asin ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Abbott"))

#### Gibert Lab
Pgm_lmers_lat$Pgm_T4_Gibert_lmer_lat <- lmer(PercT4_asin ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Gibert"))

Pgm_lmers_lat$Pgm_T5_Gibert_lmer_lat <- lmer(PercT5_asin ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Gibert"))

Pgm_lmers_lat$Pgm_T6_Gibert_lmer_lat <- lmer(PercT6_asin ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Gibert"))

# singular fit, removed population
Pgm_lmers_lat$Pgm_Total_Gibert_lmer_lat <- lmer(TotalPerc_asin ~ Latitude + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Gibert"))

#### Schmidt Lab
Pgm_lmers_lat$Pgm_T4_Schmidt_lmer_lat <- lmer(PercT4_asin ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Schmidt"))

Pgm_lmers_lat$Pgm_T5_Schmidt_lmer_lat <- lmer(PercT5_asin ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Schmidt"))

Pgm_lmers_lat$Pgm_T6_Schmidt_lmer_lat <- lmer(PercT6_asin ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Schmidt"))

Pgm_lmers_lat$Pgm_Total_Schmidt_lmer_lat <- lmer(TotalPerc_asin ~ Latitude + (1|Population) + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Schmidt"))


# save output list
saveRDS(Pgm_lmers_lat, file = file.path(lmer_dir, out_dir, "Pgm_lmers_lat.rds"))









############# output all lmers as Rdata ############# 

all_lmers_lat <- ls()[grep("lmers", ls())]
save(all_lmers_lat, file = "LinearModelsLat/all_lmers_lat.Rdata")



######### combine all lmers into one list ###########

rdsBatchReaderToList <- function(...) {
  temp = list.files(...)
  tnames <- str_split(temp, "/", simplify = T)[,3]
  tnames <- str_replace(tnames, ".rds", "")
  tlist <- lapply(temp, readRDS)
  names(tlist) <- tnames
  return(tlist)
}

all_lmers_lat <- rdsBatchReaderToList(path = "LinearModelsLat", recursive = T, full.names = T, pattern = "lmers_lat.rds")

all_lmers_lat <- unlist(all_lmers_lat, recursive=FALSE)
names(all_lmers_lat) <- str_split(names(all_lmers_lat), "\\.", simplify = T)[,2]

saveRDS(all_lmers_lat, file = "LinearModelsLat/all_lmers_list_lat.rds")



############# check linear models residuals ############# 


lmers <- list.files(path = "LinearModelsLat", recursive = T, full.names = T, pattern = "lmers_lat.rds")


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



############# output all lmers summaries and anovas as global lists ############# 


all_lmers_lat_summary <- lapply(all_lmers_lat, summary)
saveRDS(all_lmers_lat_summary, file = "LinearModelsLat/all_lmers_list_lat_summary.rds")

all_lmers_lat_anova <- lapply(all_lmers_lat, anova)
saveRDS(all_lmers_lat_anova, file = "LinearModelsLat/all_lmers_list_lat_anova.rds")



############# output all lmers summaries and anovas by trait ############# 


lmers <- list.files(path = "LinearModelsLat", recursive = T, full.names = T, pattern = "lmers_lat.rds")


for (i in 1:length(lmers)){
  f <- lmers[i]
  s_out_rds <- sub(".rds", "_summary.rds", f)
  s_out_txt <- sub(".rds", "_summary.txt", f)
  a_out_rds <- sub(".rds", "_anova.rds", f)
  a_out_txt <- sub(".rds", "_anova.txt", f)
  m <- readRDS(f)
  s <- lapply(m, summary)
  saveRDS(s, file = s_out_rds)
  capture.output(s, file = s_out_txt)
  a <- lapply(m, anova)
  saveRDS(a, file = a_out_rds)
  capture.output(a, file = a_out_txt)
}




############# output all lmers summaries and anovas by trait and lab ############# 

lmers <- list.files(path = "LinearModelsLat", recursive = T, full.names = T, pattern = "lmers_lat.rds")

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
  for (j in 1:length(n)){
    s <- summary(m[[j]])
    saveRDS(s, file = s_out_rds[j])
    capture.output(s, file = s_out_txt[j])
    a <- anova(m[[j]])
    saveRDS(a, file = a_out_rds[j])
    capture.output(a, file = a_out_txt[j])
  } 
}



############ get all lmers pvalues


all_lmers_lat_anova <- readRDS("LinearModelsLat/all_lmers_list_lat_anova.rds")

lat_pvalues <- bind_rows(Trait = names(all_lmers_lat_anova), P_lat = lapply(all_lmers_lat_anova, function(x) x$P[1]) %>% unlist())
lat_pvalues$Trait <- sub("_lmer_lat", "", lat_pvalues$Trait)
lat_pvalues$Trait <- sub("_lm_lat", "", lat_pvalues$Trait)


write.csv(lat_pvalues, "LinearModelsLat/all_lmers_lat_pvalues.csv", row.names = F)


















