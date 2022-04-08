
###############################################################
###################### SURVIVAL ANALYSES ######################
###############################################################



# Running survival analyses takes too long and clogs the computer / knitting. 
# Thus, all were done externally and output files are shown in .Rmd file. 




##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(ggpubr)
library(coxme)
library(survival)
library(survminer)




##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")


##### load data
droseu <- readRDS("Data/droseu_master_list_2022-04-05.rds")

##### create output directory
surv_dir <- "SurvivalAnalyses"
dir.create(surv_dir, showWarnings = F) 



### STARVATION ###

####

SR_F_coxme_Gonzalez <- coxme(Surv(AgeAtDeath_hours, Censor) ~ Population + (1|Batch) + (1|Population/Line) , data = filter(d_SR_surv, Supervisor.PI == "Gonzalez", Sex == "F"))




capture.output(summary(SR_F_coxme_Gonzalez), file = "SR_F_coxme_Gonzalez_sum.txt")
capture.output(anova(SR_F_coxme_Gonzalez), file = "SR_F_coxme_Gonzalez.txt")

SR_M_coxme_Gonzalez <- coxme(Surv(AgeAtDeath_hours, Censor) ~ Population + (1|Batch) + (1|Population/Line) , data = filter(d_SR_surv, Supervisor.PI == "Gonzalez", Sex == "M"))
capture.output(summary(SR_M_coxme_Gonzalez), file = "SR_M_coxme_Gonzalez_sum.txt")
capture.output(anova(SR_M_coxme_Gonzalez), file = "SR_M_coxme_Gonzalez.txt")

SR_F_coxme_Onder <- coxme(Surv(AgeAtDeath_hours, Censor) ~ Population + (1|Batch) + (1|Population/Line) , data = filter(d_SR_surv, Supervisor.PI == "Onder", Sex == "F"))
capture.output(summary(SR_F_coxme_Onder), file = "SR_F_coxme_Onder_sum.txt")
capture.output(anova(SR_F_coxme_Onder), file = "SR_F_coxme_Onder.txt")

SR_M_coxme_Onder <- coxme(Surv(AgeAtDeath_hours, Censor) ~ Population + (1|Batch) + (1|Population/Line) , data = filter(d_SR_surv, Supervisor.PI == "Onder", Sex == "M"))
capture.output(summary(SR_M_coxme_Onder), file = "SR_M_coxme_Onder_sum.txt")
capture.output(anova(SR_M_coxme_Onder), file = "SR_M_coxme_Onder.txt")

SR_F_coxme_Pasyukova <- coxme(Surv(AgeAtDeath_hours, Censor) ~ Population + (1|Batch) + (1|Population/Line) , data = filter(d_SR_surv, Supervisor.PI == "Pasyukova", Sex == "F"))
capture.output(summary(SR_F_coxme_Pasyukova), file = "SR_F_coxme_Pasyukova_sum.txt")
capture.output(anova(SR_F_coxme_Pasyukova), file = "SR_F_coxme_Pasyukova.txt")

SR_M_coxme_Pasyukova <- coxme(Surv(AgeAtDeath_hours, Censor) ~ Population + (1|Batch) + (1|Population/Line) , data = filter(d_SR_surv, Supervisor.PI == "Pasyukova", Sex == "M"))
capture.output(summary(SR_M_coxme_Pasyukova), file = "SR_M_coxme_Pasyukova_sum.txt")
capture.output(anova(SR_M_coxme_Pasyukova), file = "SR_M_coxme_Pasyukova.txt")

################################################################################

### LIFESPAN ###

d_LS_L <- read.csv("LSL_MasterSheet_Oct21.csv")

d_LS_L$Supervisor.PI <- as.factor(d_LS_L$Supervisor.PI)
d_LS_L$Diet <- as.factor(d_LS_L$Diet)
d_LS_L$Batch <- as.factor(d_LS_L$Batch)
d_LS_L$Population <- as.factor(d_LS_L$Population)
d_LS_L$Population_Lat <- factor(d_LS_L$Population, levels= c("YE","RE","GI","MU","MA","UM","KA","VA","AK"))
d_LS_L$Population_Lon <- factor(d_LS_L$Population, levels= c("RE","GI","KA","MU","MA","AK","UM","YE","VA"))
d_LS_L$Population_Alt <- factor(d_LS_L$Population, levels= c("KA","AK","GI","RE","UM","VA","MU","MA","YE"))
d_LS_L$Line <- as.factor(d_LS_L$Line)
d_LS_L$ReplicateVial <- as.factor(d_LS_L$ReplicateVial)
d_LS_L$LSL_AgeAtDeath_days <- as.numeric(d_LS_L$LSL_AgeAtDeath_days)
d_LS_L$Censor <- as.factor(d_LS_L$Censor)


d_LS_P <- read.csv("LSP_MasterSheet_Oct21.csv")

d_LS_P$Supervisor.PI <- as.factor(d_LS_P$Supervisor.PI)
d_LS_P$Diet <- as.factor(d_LS_P$Diet)
d_LS_P$Batch <- as.factor(d_LS_P$Batch)
d_LS_P$Population <- as.factor(d_LS_P$Population)
d_LS_P$Population_Lat <- factor(d_LS_P$Population, levels= c("YE","RE","GI","MU","MA","UM","KA","VA","AK"))
d_LS_P$Population_Lon <- factor(d_LS_P$Population, levels= c("RE","GI","KA","MU","MA","AK","UM","YE","VA"))
d_LS_P$Population_Alt <- factor(d_LS_P$Population, levels= c("KA","AK","GI","RE","UM","VA","MU","MA","YE"))
d_LS_P$ReplicateCage <- as.factor(d_LS_P$ReplicateCage)
d_LS_P$LSP_AgeAtDeath_days <- as.numeric(d_LS_P$LSP_AgeAtDeath_days)
d_LS_P$Censor <- as.factor(d_LS_P$Censor)

d_LS_P_surv <- d_LS_P %>% mutate(Censor = ifelse(Censor == 0, 1, 0))
d_LS_L_surv <- d_LS_L %>% mutate(Censor = ifelse(Censor == 0, 1, 0))

###

LS_P_F_coxme_Flatt <- coxme(Surv(LSP_AgeAtDeath_days, Censor) ~ Population + (1|Population/ReplicateCage), data = filter(d_LS_P_surv, Supervisor.PI == "Flatt", Sex == "F"))
capture.output(summary(LS_P_F_coxme_Flatt), file = "LS_P_F_coxme_Flatt_sum.txt")
capture.output(anova(LS_P_F_coxme_Flatt), file = "LS_P_F_coxme_Flatt.txt")
           
LS_P_M_coxme_Flatt <- coxme(Surv(LSP_AgeAtDeath_days, Censor) ~ Population + (1|Population/ReplicateCage), data = filter(d_LS_P_surv, Supervisor.PI == "Flatt", Sex == "M"))
capture.output(summary(LS_P_M_coxme_Flatt), file = "LS_P_M_coxme_Flatt_sum.txt")
capture.output(anova(LS_P_M_coxme_Flatt), file = "LS_P_M_coxme_Flatt.txt")
           
LS_L_F_coxme_Parsch <- coxme(Surv(LSL_AgeAtDeath_days, Censor) ~ Population + (1|Batch) + (1|Population/Line) , data = filter(d_LS_L_surv, Supervisor.PI == "Parsch", Sex == "F"))
capture.output(summary(LS_L_F_coxme_Parsch), file = "LS_L_F_coxme_Parsch_sum.txt")
capture.output(anova(LS_L_F_coxme_Parsch), file = "LS_L_F_coxme_Parsch.txt")
  
LS_L_M_coxme_Parsch <- coxme(Surv(LSL_AgeAtDeath_days, Censor) ~ Population + (1|Batch) + (1|Population/Line) , data = filter(d_LS_L_surv, Supervisor.PI == "Parsch", Sex == "M"))
capture.output(summary(LS_L_M_coxme_Parsch), file = "LS_L_M_coxme_Parsch_sum.txt")
capture.output(anova(LS_L_M_coxme_Parsch), file = "LS_L_M_coxme_Parsch.txt")
  
LS_L_F_coxme_Pasyukova <- coxme(Surv(LSL_AgeAtDeath_days, Censor) ~ Population + (1|Batch) + (1|Population/Line) , data = filter(d_LS_L_surv, Supervisor.PI == "Pasyukova", Sex == "F"))
capture.output(summary(LS_L_F_coxme_Pasyukova), file = "LS_L_F_coxme_Pasyukova_sum.txt")
capture.output(anova(LS_L_F_coxme_Pasyukova), file = "LS_L_F_coxme_Pasyukova.txt")

LS_L_M_coxme_Pasyukova <- coxme(Surv(LSL_AgeAtDeath_days, Censor) ~ Population + (1|Batch) + (1|Population/Line) , data = filter(d_LS_L_surv, Supervisor.PI == "Pasyukova", Sex == "M"))
capture.output(summary(LS_L_M_coxme_Pasyukova), file = "LS_L_M_coxme_Pasyukova_sum.txt")
capture.output(anova(LS_L_M_coxme_Pasyukova), file = "LS_L_M_coxme_Pasyukova.txt")
  
####

### HEATSHOCK ### 

d_HSM <- read.csv("HSM_MasterSheet_Oct21.csv")

d_HSM$Supervisor.PI <- as.factor(d_HSM$Supervisor.PI)
d_HSM$Diet <- as.factor(d_HSM$Diet)
d_HSM$Batch <- as.factor(d_HSM$Batch)
d_HSM$Population_Lat <- factor(d_HSM$Population, levels= c("YE","RE","GI","MU","MA","UM","KA","VA","AK"))
d_HSM$Population_Lon <- factor(d_HSM$Population, levels= c("RE","GI","KA","MU","MA","AK","UM","YE","VA"))
d_HSM$Population_Alt <- factor(d_HSM$Population, levels= c("KA","AK","GI","RE","UM","VA","MU","MA","YE"))
d_HSM$Line <- as.factor(d_HSM$Line)
d_HSM$Sex <- as.factor(d_HSM$Sex)
d_HSM$ReplicateVial <- as.factor(d_HSM$ReplicateVial)
d_HSM$TimeDeath_min <- as.numeric(d_HSM$TimeDeath_min)
d_HSM$Censor <- as.numeric(d_HSM$Censor)
str(d_HSM)

d_HSM_surv <- d_HSM %>% mutate(Censor = ifelse(Censor == 0, 1, 0))

HSM_F_coxme_Parsch <- coxme(Surv(TimeDeath_min, Censor) ~ Population + (1|Batch) + (1|Population/Line) , data = filter(d_HSM_surv, Supervisor.PI == "Parsch", Sex == "F"))
capture.output(summary(HSM_F_coxme_Parsch), file = "HSM_F_coxme_Parsch_sum.txt")
capture.output(anova(HSM_F_coxme_Parsch), file = "HSM_F_coxme_Parsch.txt")

HSM_M_coxme_Parsch <- coxme(Surv(TimeDeath_min, Censor) ~ Population + (1|Batch) + (1|Population/Line) , data = filter(d_HSM_surv, Supervisor.PI == "Parsch", Sex == "M"))
capture.output(summary(HSM_M_coxme_Parsch), file = "HSM_M_coxme_Parsch_sum.txt")
capture.output(anova(HSM_M_coxme_Parsch), file = "HSM_M_coxme_Parsch.txt")

HSM_F_coxme_Vieira <- coxme(Surv(TimeDeath_min, Censor) ~ Population + (1|Batch) + (1|Population/Line) , data = filter(d_HSM_surv, Supervisor.PI == "Vieira", Sex == "F"))
capture.output(summary(HSM_F_coxme_Vieira), file = "HSM_F_coxme_Vieira_sum.txt")
capture.output(anova(HSM_F_coxme_Vieira), file = "HSM_F_coxme_Vieira.txt")

HSM_M_coxme_Vieira <- coxme(Surv(TimeDeath_min, Censor) ~ Population + (1|Batch) + (1|Population/Line) , data = filter(d_HSM_surv, Supervisor.PI == "Vieira", Sex == "M"))
capture.output(summary(HSM_M_coxme_Vieira), file = "HSM_M_coxme_Vieira_sum.txt")
capture.output(anova(HSM_M_coxme_Vieira), file = "HSM_M_coxme_Vieira.txt")


### CCRT ###

d_CCRT <- read.csv("CCRT_MasterSheet_Oct21.csv")

d_CCRT$Supervisor.PI <- as.factor(d_CCRT$Supervisor.PI)
d_CCRT$Diet <- as.factor(d_CCRT$Diet)
d_CCRT$Batch <- as.factor(d_CCRT$Batch)
d_CCRT$Population_Lat <- factor(d_CCRT$Population, levels= c("YE","RE","GI","MU","MA","UM","KA","VA","AK"))
d_CCRT$Population_Lon <- factor(d_CCRT$Population, levels= c("RE","GI","KA","MU","MA","AK","UM","YE","VA"))
d_CCRT$Population_Alt <- factor(d_CCRT$Population, levels= c("KA","AK","GI","RE","UM","VA","MU","MA","YE"))
d_CCRT$Line <- as.factor(d_CCRT$Line)
d_CCRT$Sex <- as.factor(d_CCRT$Sex)
d_CCRT$ReplicateVial <- as.factor(d_CCRT$ReplicateVial)
d_CCRT$CCRT_seconds <- as.numeric(d_CCRT$CCRT_seconds)
d_CCRT$Censor <- as.numeric(d_CCRT$Censor)

d_CCRT_surv <- d_CCRT %>% mutate(Censor = ifelse(Censor == 0, 1, 0))

CCRT_F_coxme_Vieira <- coxme(Surv(CCRT_seconds, Censor) ~ Population + (1|Batch) + (1|Population/Line) , data = filter(d_CCRT_surv, Supervisor.PI == "Vieira", Sex == "F"))
capture.output(summary(CCRT_F_coxme_Vieira), file = "CCRT_F_coxme_Vieira_sum.txt")
capture.output(anova(CCRT_F_coxme_Vieira), file = "CCRT_F_coxme_Vieira.txt")

CCRT_M_coxme_Vieira <- coxme(Surv(CCRT_seconds, Censor) ~ Population + (1|Batch) + (1|Population/Line) , data = filter(d_CCRT_surv, Supervisor.PI == "Vieira", Sex == "M"))
capture.output(summary(CCRT_M_coxme_Vieira), file = "CCRT_M_coxme_Vieira_sum.txt")
capture.output(anova(CCRT_M_coxme_Vieira), file = "CCRT_M_coxme_Vieira.txt")



