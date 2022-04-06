

##################################################################### 
######################  DESCRIPTIVE STATISTICS ######################
#####################################################################




##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
#library(lme4)
#library(lsmeans)
#library(afex)
#library(multcomp)


##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")


##### source functions
#source("Functions/lab_correlations_functions.R")

##### load data
droseu <- readRDS("Data/droseu_master_list_2022-04-05.rds")

##### create output directory
desc_dir <- "DescriptiveStatistics"
dir.create(desc_dir, showWarnings = F) 


##### define functions

std_err <- function(x) sd(x)/sqrt(length(x))
coef_var <- function(x) sd(x)/mean(x)
estimate_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}



var_list <- c("CCRT_seconds", "ZT_hours_MESA", "ZT_hours_LSPR", "Period_MESA", "Period_LSPR", "Rhythmicity_LSPR_amplitude", "Rhythmicity_JTK_p_BH_corrected", "CSM_PropDead_ED", "Prop_Max_Stage7", "Prop_Max_Stage8", "Prop_Max_Stage9", "DT_EggAdult", "DT_EggPupa", "DW_micrograms", "NumberOfAdultsEclosed", "TimeDeath_min", "Period", "CircPhase", "AbsPhase", "ND", "Activity", "LSL_AgeAtDeath_days", "LSM_AgeAtDeath_days", "LSP_AgeAtDeath_days", "PercT4", "PercT5", "PercT6", "TotalPerc", "ScoreT4", "ScoreT5", "ScoreT6", "TotalScore", "AgeAtDeath_hours", "TL_micrometers", "ProportionEggtoAdultSurvival", "CentroidSizeLeft_micrometers", "CentroidSizeRight_micrometers")

com_var <- c("Supervisor.PI", "Batch", "Population", "Line", "Sex")

names(droseu) <- toupper(names(droseu))

names(droseu)[names(droseu) == "DIA"] <- "Dia"
names(droseu)[names(droseu) == "FEC"] <- "Fec"
names(droseu)[names(droseu) == "PGM"] <- "Pgm"
names(droseu)[names(droseu) == "PGM2"] <- "Pgm2"
names(droseu)[names(droseu) == "VIA"] <- "Via"

droseu_long <- droseu
names_d <- names(droseu_long)
for (i in 1:length(droseu_long)){
  d <- droseu_long[[i]]
  if (!"Line" %in% colnames(d)) d$Line = NA
  if (!"Batch" %in% colnames(d)) dd$Batch = NA
  if (!"Sex" %in% colnames(d)) d$Sex = NA
  d <- d[, colnames(d) %in% c(com_var, var_list)]
  d <- pivot_longer(d, cols = -all_of(com_var), names_to = "Trait", values_to = "Value")
  d$Trait_name <- names_d[i]
  d <- split(d, d$Trait)
  droseu_long[[i]] <- d
}

droseu_long <- bind_rows(unlist(droseu_long, recursive = F))


##### Line level with Batch

line_wbatch <- group_by(droseu_long, Trait_name, Trait, Supervisor.PI, Batch, Sex, Population, Line) %>% summarise_at(vars(Value), list(Mean = mean, SD = sd, Median = median, Min = min, Max = max, SE = std_err, CV = coef_var, Mode = estimate_mode)) %>% ungroup

n <- unique(paste0(line_wbatch$Trait_name, "_", line_wbatch$Trait))
n[grep("CCRT_", n)] <- "CCRT"
n[grep("CSM_", n)] <- "CSM"
n[grep("DTP_", n)] <- "DT_P"
n[grep("DTA_", n)] <- "DT_A"
n[grep("DW_", n)] <- "DW"
n[grep("Fec_", n)] <- "Fec"
n[grep("HSM_", n)] <- "HSM"
n[grep("LSL_", n)] <- "LS_L"
n[grep("LSP_", n)] <- "LS_P"
n[grep("LSM_", n)] <- "LS_M"
n[grep("TL_", n)] <- "TL"
n[grep("SR_", n)] <- "SR"
n[grep("Via_", n)] <- "Via"
n[grep("WA_CentroidSizeLeft", n)] <- "WA_Left"
n[grep("WA_CentroidSizeRight", n)] <- "WA_Right"
line_wbatch <- group_split(line_wbatch, Trait_name, Trait)
names(line_wbatch) <- n

for (i in 1:length(line_wbatch)){
  d <- line_wbatch[[i]]
  out <- paste0("table_", n[i], "_Line_wbatch.csv")
  d <- select(d, -c(Trait, Trait_name))
  write.csv(d, file.path(desc_dir, out), row.names = F)
}

saveRDS(line_wbatch, file.path(desc_dir, "all_table_Line_wbatch.rds"))


##### Line level without Batch

line_wobatch <- group_by(droseu_long, Trait_name, Trait, Supervisor.PI, Sex, Population, Line) %>% summarise_at(vars(Value), list(Mean = mean, SD = sd, Median = median, Min = min, Max = max, SE = std_err, CV = coef_var, Mode = estimate_mode)) %>% ungroup

#n <- unique(paste0(line_wobatch$Trait_name, "_", line_wobatch$Trait)) 
line_wobatch <- group_split(line_wobatch, Trait_name, Trait)
names(line_wobatch) <- n

for (i in 1:length(line_wobatch)){
  d <- line_wobatch[[i]]
  out <- paste0("table_", n[i], "_Line_wobatch.csv")
  d <- select(d, -c(Trait, Trait_name))
  write.csv(d, file.path(desc_dir, out), row.names = F)
}

saveRDS(line_wobatch, file.path(desc_dir, "all_table_Line_wobatch.rds"))



##### Population level with Batch

pop_wbatch <- group_by(droseu_long, Trait_name, Trait, Supervisor.PI, Batch, Sex, Population) %>% summarise_at(vars(Value), list(Mean = mean, SD = sd, Median = median, Min = min, Max = max, SE = std_err, CV = coef_var, Mode = estimate_mode)) %>% ungroup

#n <- unique(paste0(pop_wbatch$Trait_name, "_", pop_wbatch$Trait)) 
pop_wbatch <- group_split(pop_wbatch, Trait_name, Trait)
names(pop_wbatch) <- n

for (i in 1:length(pop_wbatch)){
  d <- pop_wbatch[[i]]
  out <- paste0("table_", n[i], "_Pop_wbatch.csv")
  d <- select(d, -c(Trait, Trait_name))
  write.csv(d, file.path(desc_dir, out), row.names = F)
}

saveRDS(pop_wbatch, file.path(desc_dir, "all_table_Pop_wbatch.rds"))



##### Population level without Batch

pop_wobatch <- group_by(droseu_long, Trait_name, Trait, Supervisor.PI, Sex, Population) %>% summarise_at(vars(Value), list(Mean = mean, SD = sd, Median = median, Min = min, Max = max, SE = std_err, CV = coef_var, Mode = estimate_mode)) %>% ungroup

#n <- unique(paste0(pop_wobatch$Trait_name, "_", pop_wobatch$Trait)) 
pop_wobatch <- group_split(pop_wobatch, Trait_name, Trait)
names(pop_wobatch) <- n

for (i in 1:length(pop_wobatch)){
  d <- pop_wobatch[[i]]
  out <- paste0("table_", n[i], "_Pop_wobatch.csv")
  d <- select(d, -c(Trait, Trait_name))
  write.csv(d, file.path(desc_dir, out), row.names = F)
}

saveRDS(pop_wobatch, file.path(desc_dir, "all_table_Pop_wobatch.rds"))



