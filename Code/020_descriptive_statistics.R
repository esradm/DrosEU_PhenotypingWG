

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


removeXtraCols <- function(x) {
  if (unique(x$Trait) == "CETS") x <- subset(x, select = -c(SD, Median, Min, Max, SE, CV, Mode))
  x <- x[,!colnames(x) %in% c("Trait", "Trait_name")]
  x <- relocate(x, Condition, .after = Supervisor.PI)
  if (is.na(unique(x$Condition))) x <- subset(x, select = -c(Condition))
  if (is.na(unique(x$Sex))[1]) x <- subset(x, select = -c(Sex))
  if ("Line" %in% colnames(x)) {
    if (is.na(unique(x$Line)[1]) | unique(x$Line)[1] == "mixed population") x <- subset(x, select = -c(Line)) }
  return(x)
}


##### change names to match what was done previously in the rmd file

names(droseu) <- toupper(names(droseu))
names(droseu)[names(droseu) == "DIA"] <- "Dia"
names(droseu)[names(droseu) == "FEC"] <- "Fec"
names(droseu)[names(droseu) == "PGM"] <- "Pgm"
names(droseu)[names(droseu) == "VIA"] <- "Via"


##### define trait and common variables

trait_var <- c("CCRT_seconds", "ZT_hours_MESA", "ZT_hours_LSPR", "Period_MESA", "Period_LSPR", "Rhythmicity_LSPR_amplitude", "Rhythmicity_JTK_p_BH_corrected", "CSM_PropDead_ED", "Prop_Max_Stage7", "Prop_Max_Stage8", "Prop_Max_Stage9", "DT_EggAdult", "DT_EggPupa", "DW_micrograms", "NumberOfAdultsEclosed", "TimeDeath_min", "Period", "CircPhase", "AbsPhase", "ND", "Activity", "LSL_AgeAtDeath_days", "LSM_AgeAtDeath_days", "LSP_AgeAtDeath_days", "PercT4", "PercT5", "PercT6", "TotalPerc", "AgeAtDeath_hours", "TL_micrometers", "ProportionEggtoAdultSurvival", "CentroidSizeLeft_micrometers", "CentroidSizeRight_micrometers")

com_var <- c("Supervisor.PI", "Batch", "Population", "Line", "Sex", "Condition")


##### add extra columns for easier batch processing and split data by single trait variables

droseu_long <- droseu
for (i in 1:length(droseu_long)){
  d <- droseu_long[[i]]
  if (!"Line" %in% colnames(d)) d$Line = NA
  #if (!"Batch" %in% colnames(d)) d$Batch = NA
  if (!"Sex" %in% colnames(d)) d$Sex = NA
  if (!"Condition" %in% colnames(d)) d$Condition = NA
  d <- d[, colnames(d) %in% c(com_var, trait_var)]
  d <- pivot_longer(d, cols = -all_of(com_var), names_to = "Trait_name", values_to = "Value")
  d$Trait <- names(droseu_long)[i]
  d <- split(d, d$Trait_name)
  droseu_long[[i]] <- d
}

##### combine all traits into one data.frame

droseu_long <- bind_rows(unlist(droseu_long, recursive = F)) 


#### define output names to match what was done previously in rmd file

n <- unique(paste(droseu_long$Trait, droseu_long$Condition, droseu_long$Trait_name, sep ="_"))
n <- sub("_NA", "", n)
n[grep("CCRT_", n)] <- "CCRT"
n[grep("CSM_", n)] <- "CSM"
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



##### Line level with Batch

line_wbatch <- group_by(droseu_long, Trait, Trait_name, Condition, Supervisor.PI, Batch, Sex, Population, Line) %>% summarise_at(vars(Value), list(Mean = mean, SD = sd, Median = median, Min = min, Max = max, SE = std_err, CV = coef_var, Mode = estimate_mode)) %>% ungroup

line_wbatch <- group_split(line_wbatch, Trait, Trait_name, Condition)
line_wbatch <- lapply(line_wbatch, removeXtraCols)
names(line_wbatch) <- n

saveRDS(line_wbatch, file.path(desc_dir, "all_table_Line_wbatch.rds"))

for (i in 1:length(line_wbatch)){
  out <- paste0("table_", n[i], "_Line_wbatch.csv")
  write.csv(line_wbatch[[i]], file.path(desc_dir, out), row.names = F)
}



##### Line level without Batch

line_wobatch <- group_by(droseu_long, Trait, Trait_name, Condition, Supervisor.PI, Sex, Population, Line) %>% summarise_at(vars(Value), list(Mean = mean, SD = sd, Median = median, Min = min, Max = max, SE = std_err, CV = coef_var, Mode = estimate_mode)) %>% ungroup

line_wobatch <- group_split(line_wobatch, Trait, Trait_name, Condition)
line_wobatch <- lapply(line_wobatch, removeXtraCols)
names(line_wobatch) <- n

saveRDS(line_wobatch, file.path(desc_dir, "all_table_Line_wobatch.rds"))

for (i in 1:length(line_wobatch)){
  out <- paste0("table_", n[i], "_Line_wobatch.csv")
  write.csv(line_wobatch[[i]], file.path(desc_dir, out), row.names = F)
}




##### Population level with Batch


pop_wbatch <- group_by(droseu_long, Trait, Trait_name, Condition, Supervisor.PI, Batch, Sex, Population) %>% summarise_at(vars(Value), list(Mean = mean, SD = sd, Median = median, Min = min, Max = max, SE = std_err, CV = coef_var, Mode = estimate_mode)) %>% ungroup

pop_wbatch <- group_split(pop_wbatch, Trait, Trait_name, Condition)
pop_wbatch <- lapply(pop_wbatch, removeXtraCols)
names(pop_wbatch) <- n

saveRDS(pop_wbatch, file.path(desc_dir, "all_table_Pop_wbatch.rds"))

for (i in 1:length(pop_wbatch)){
  out <- paste0("table_", n[i], "_Pop_wbatch.csv")
  write.csv(pop_wbatch[[i]], file.path(desc_dir, out), row.names = F)
}




##### Population level without Batch

pop_wobatch <- group_by(droseu_long, Trait, Trait_name, Condition, Supervisor.PI, Sex, Population) %>% summarise_at(vars(Value), list(Mean = mean, SD = sd, Median = median, Min = min, Max = max, SE = std_err, CV = coef_var, Mode = estimate_mode)) %>% ungroup

pop_wobatch <- group_split(pop_wobatch, Trait, Trait_name, Condition)
pop_wobatch <- lapply(pop_wobatch, removeXtraCols)
names(pop_wobatch) <- n

saveRDS(pop_wobatch, file.path(desc_dir, "all_table_Pop_wobatch.rds"))

for (i in 1:length(pop_wobatch)){
  out <- paste0("table_", n[i], "_Pop_wobatch.csv")
  write.csv(pop_wobatch[[i]], file.path(desc_dir, out), row.names = F)
}








