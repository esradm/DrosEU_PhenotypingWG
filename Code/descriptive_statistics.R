

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



var_list <- c("CCRT_seconds", "ZT_hours_MESA", "ZT_hours_LSPR", "Period_MESA", "Period_LSPR", "Rhythmicity_LSPR_amplitude", "Rhythmicity_JTK_p_BH_corrected", "CSM_PropDead_ED", "Prop_Max_Stage7", "Prop_Max_Stage8", "Prop_Max_Stage9", "DT_EggAdult", "DT_EggPupa", "DW_micrograms", "NumberOfAdultsEclosed", "TimeDeath_min", "Period", "CircPhase", "AbsPhase", "ND", "Activity", "LSL_AgeAtDeath_days", "LSM_AgeAtDeath_days", "LSP_AgeAtDeath_days", "PercT4", "PercT5", "PercT6", "TotalPerc", "ScoreT5", "ScoreT6", "ScoreT7", "TotalScore", "AgeAtDeath_hours", "TL_micrometers", "ProportionEggtoAdultSurvival", "CentroidSizeLeft_micrometers", "CentroidSizeRight_micrometers")

com_var <- c("Supervisor.PI", "Batch", "Population", "Line", "Sex")



d <- droseu$wa

d <- d[, colnames(d) %in% c(com_var, var_list)]



pivot_longer(d, cols = dplyr::select(-all_of(com_var)), names_to = "Trait", values_to = "Value")



#Supervisor.PI, Batch, Population, Line, Sex



droseu2 <- droseu
for (i in (1:length(droseu2))){
  if (!"Line" %in% colnames(droseu2[[i]])) droseu2[[i]]$Line = NA
  if (!"Batch" %in% colnames(droseu2[[i]])) droseu2[[i]]$Batch = NA
  if (!"Sex" %in% colnames(droseu2[[i]])) droseu2[[i]]$Sex = NA
}


for (i in var_list){
  var_sel <- c(com_var, i)
  for (j in 1:length(droseu2)){
    a <- droseu2[[j]][droseu2[[j]][,var_sel]
  }
}





for (i in (1:length(droseu2))){
  droseu2[[i]] <- group_by(droseu2[[i]], Supervisor.PI, Batch, Sex, Population, Line) %>% summarise_at(vars(-group_cols()), list(Mean = mean, SD = sd, Median = median, Min = min, Max = max, SE = std_err, CV = coef_var, Mode = estimate_mode))
}


                   
                   
  summarise(across(where(is.numeric), ~ mean(.x)))



  summarise_at(vars(-group_cols()), Mean = mean)
  
  
  summarise(across(where(is.numeric)), ~ list(Mean = mean, SD = sd, Median = median, Min = min, Max = max, SE = std_err, CV = coef_var, Mode = estimate_mode))
summarise_at(vars(-group_cols(), ...), myoperation)

group_by(Trait, Supervisor.PI, Diet, Sex, Population, Line) %>%
  summarise(Replicate_count = length(unique(ReplicateVial)),
            Replicate_ids = paste(sort(unique(ReplicateVial)), collapse = "; "),
            Observation_count = n(), 
            Batch_count = length(unique(Batch)),
            Batch_ids = paste(sort(unique(Batch)), collapse = "; "),
            .groups = "drop")



table_Via_Line_wbatch <- write.csv(d_Via %>% group_by(Supervisor.PI, Batch, Population, Line) %>% summarise_at(vars(ProportionEggtoAdultSurvival), list(Mean = mean, SD= sd, Median = median, Min = min, Max = max, SE = std_err, CV = coef_var, Mode = estimate_mode)), file = "Viability/table_Via_Line_wbatch.csv", row.names = T)






