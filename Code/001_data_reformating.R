

rm(list = ls())

setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/Data")

library(tidyverse)
library(nortest)
library(MetBrewer)

##### define functions

# Wrapper for read.csv function. Takes list.files arguments. Each .csv file is turned into a list element. List element names are extracted from .csv files names. StringsAsFactors = F and na.strings = c("", " ", "NA"). 
csvBatchReaderToList3 <- function(dir, ...) {
  f <- list.files(path = dir, ...)
  n <- sub("_MasterSheet_.*", "", f)
  n <- tolower(n)
  p <- paste(dir, f, sep = "/")
  flist <- lapply(p, function(x) read.csv(x, stringsAsFactors = F, na.strings = c("", " ", "NA")))
  names(flist) <- n
  return(flist)
}


##### batch load all the .csv files, make sure that the path matches the most up to date directory

droseu <- csvBatchReaderToList3(dir = "MasterSheets_Feb22_git", pattern = "*.csv", full.names = F)



##### change the variables that are common to all traits

droseu <- lapply(droseu, mutate, 
                 Supervisor.PI = as.factor(Supervisor.PI), 
                 Diet = as.factor(Diet),  
                 Batch = as.factor(Batch), 
                 Population = as.factor(Population),
                 Population_Lat = factor(Population, levels = c("YE","RE","GI","MU","MA","UM","KA","VA","AK")),
                 Population_Lon = factor(Population, levels = c("RE","GI","KA","MU","MA","AK","UM","YE","VA")),
                 Population_Alt = factor(Population, levels = c("KA","AK","GI","RE","UM","VA","MU","MA","YE")),
                 Country = as.factor(Country),
                 Latitude = as.numeric(Latitude),
                 Longitude = as.numeric(Longitude),
                 Altitude = as.numeric(Altitude)
                 )


##### change the variables that are not common to all traits

for (i in 1:length(droseu)) {
  if ("Line" %in% colnames(droseu[[i]])) {
    droseu[[i]]$Line <- as.factor(droseu[[i]]$Line) }
  if ("Sex" %in% colnames(droseu[[i]])) {
    droseu[[i]]$Sex <- as.factor(droseu[[i]]$Sex) }
  if ("ReplicateVial" %in% colnames(droseu[[i]])) {
    droseu[[i]]$ReplicateVial <- as.factor(droseu[[i]]$ReplicateVial) }
  if ("ReplicateVialOld" %in% colnames(droseu[[i]])) {
    droseu[[i]]$ReplicateVialOld <- as.factor(droseu[[i]]$ReplicateVialOld) }
  if ("Censor" %in% colnames(droseu[[i]])) {
    droseu[[i]]$Censor <- as.factor(droseu[[i]]$Censor) }
  if ("ReplicateCage" %in% colnames(droseu[[i]])) {
    droseu[[i]]$ReplicateCage <- as.factor(droseu[[i]]$ReplicateCage) }
  if ("ReplicateCageOld" %in% colnames(droseu[[i]])) {
    droseu[[i]]$ReplicateCageOld <- as.factor(droseu[[i]]$ReplicateCageOld) }
  if ("ReplicateChamber" %in% colnames(droseu[[i]])) {
    droseu[[i]]$ReplicateChamber <- as.factor(droseu[[i]]$ReplicateChamber) }
  if ("ReplicateChamberOld" %in% colnames(droseu[[i]])) {
    droseu[[i]]$ReplicateChamberOld <- as.factor(droseu[[i]]$ReplicateChamberOld) }
}

##### fix read.csv behavior when Sex is F only, typically in pigmentation data

# female only traits: dia, pgm, pgm2, fec
# male only traits: la
# mixed sex traits: cet, dtp, via

droseu$dia$Sex <- sub("FALSE", "F", droseu$dia$Sex)
droseu$fec$Sex <- sub("FALSE", "F", droseu$fec$Sex)
droseu$pgm$Sex <- sub("FALSE", "F", droseu$pgm$Sex)
droseu$pgm2$Sex <- sub("FALSE", "F", droseu$pgm2$Sex)



##### calculate diapause proportions

droseu$dia <- droseu$dia %>%
  mutate(MostAdvancedStage = as.numeric(MostAdvancedStage),
         NumberOfEggs = as.numeric(NumberOfEggs),
         Max_Stage7 = ifelse(MostAdvancedStage <= 7 & NumberOfEggs == 0, 1, 0),
         Max_Stage8 = ifelse(MostAdvancedStage <= 8 & NumberOfEggs == 0, 1, 0),
         Max_Stage9 = ifelse(MostAdvancedStage <= 9 & NumberOfEggs == 0, 1, 0)) %>%
  group_by(across(-c(contains("Stage"), NumberOfEggs, Individual))) %>%
  summarise_at(vars(contains("Max_Stage")), mean) %>% 
  rename_at(vars(contains("Max_Stage")), ~ paste0("Prop_", .x)) %>%
  as.data.frame()

 

##### define all the trait variables


var_list <- c("CCRT_seconds", "ZT_hours_MESA", "ZT_hours_LSPR", "Period_MESA", "Period_LSPR", "Rhythmicity_LSPR_amplitude", "Rhythmicity_JTK_p_BH_corrected", "CSM_PropDead_ED", "Prop_Max_Stage7", "Prop_Max_Stage8", "Prop_Max_Stage9", "DT_EggAdult", "DT_EggPupa", "DW_micrograms", "NumberOfAdultsEclosed", "TimeDeath_min", "Period", "CircPhase", "AbsPhase", "ND", "Activity", "LSL_AgeAtDeath_days", "LSM_AgeAtDeath_days", "LSP_AgeAtDeath_days", "PercT4", "PercT5", "PercT6", "TotalPerc", "ScoreT4", "ScoreT5", "ScoreT6", "TotalScore", "AgeAtDeath_hours", "TL_micrometers", "ProportionEggtoAdultSurvival", "CentroidSizeLeft_micrometers", "CentroidSizeRight_micrometers")

##### turn all trait variables to numeric

for (i in 1:length(droseu)) {
  for (v in var_list) {
    if (v %in% colnames(droseu[[i]])) {
      droseu[[i]][,v] <- as.numeric(droseu[[i]][,v]) }
  }
}



##### arcsin transformation of proportion data

droseu$csm <- droseu$csm %>% 
  mutate(CSM_PropDead_ED_asin = asin(sqrt(CSM_PropDead_ED)))

droseu$dia <- droseu$dia %>% 
  mutate(Prop_Max_Stage7_asin = asin(sqrt(Prop_Max_Stage7)),
         Prop_Max_Stage8_asin = asin(sqrt(Prop_Max_Stage8)),
         Prop_Max_Stage9_asin = asin(sqrt(Prop_Max_Stage9)))

droseu$pgm <- droseu$pgm %>% 
  mutate(PercT4_asin = asin(sqrt(PercT4/100)),
         PercT5_asin = asin(sqrt(PercT5/100)),
         PercT6_asin = asin(sqrt(PercT6/100)),
         TotalPerc_asin = asin(sqrt(TotalPerc/100)))

droseu$via <- droseu$via %>% 
  mutate(ProportionEggtoAdultSurvival_asin = asin(sqrt(ProportionEggtoAdultSurvival)))


##### log2 transformation of ND data
# introduces -Inf values because two data points are 0

droseu$la <- droseu$la %>% mutate(ND_log2 = log2(ND))



##### update the trait variable list with transformed data

var_list_up <- c(var_list, "CSM_PropDead_ED_asin", "Prop_Max_Stage7_asin", "Prop_Max_Stage8_asin", "Prop_Max_Stage9_asin", "PercT4_asin", "PercT5_asin", "PercT6_asin", "TotalPerc_asin", "ProportionEggtoAdultSurvival_asin", "ND_log2")




##### run Shapiro and / or Anderson Darling tests

#shapiTest <- function(x, v) {
#  shapiro.test(select(x, all_of(v)) %>% unlist)
#}

#andersonDarlingTest <- function(x, v) {
#  ad.test(select(x, all_of(v)) %>% unlist)
#}

#traits <- vector("list", length = length(droseu))
#for (i in 1:length(droseu)) {
#  trait <- droseu[[i]]
#  trait_name <- names(droseu)[[i]]
#  trait_var <- colnames(trait)[colnames(trait) %in% var_list_up]
#  n_var <- length(trait_var)
#  if (n_var > 0) {
#    if ("Sex" %in% colnames(trait)) {
#      trait <- arrange(trait, Supervisor.PI, Sex)
#      lab_names <- unique(paste(trait$Supervisor.PI, trait$Sex, sep = "_"))
#      trait <- group_split(trait, Supervisor.PI, Sex) 
#    } else { 
#      trait <- arrange(trait, Supervisor.PI)
#      lab_names <- unique(trait$Supervisor.PI)
#      trait <- group_split(trait, Supervisor.PI) 
#    }
#    names(trait) <- lab_names
#    stats <- list()
#    for(v in 1:n_var){
#      l <- unlist(lapply(trait, nrow))
#      if (max(l) < 5000) {
#        test <- lapply(trait, shapiTest, v = trait_var[v])
#        s <- c("W_Statistic", "p")
#      } else {
#        test <- lapply(trait, andersonDarlingTest, v = trait_var[v]) 
#        s <- c("A_Statistic", "p")
#      }
#      stat <- as.data.frame(t(sapply(test, `[`, c("statistic", "p.value"))))
#      colnames(stat) <- s
#      stat$Lab <- lab_names
#      stat$Trait <- trait_name
#      stat$Variable <- trait_var[v]
#      stat$Sig <- ifelse(stat$p < 0.05, "TRUE", "FALSE")
#      stat <- relocate(stat, Trait, Lab, Variable)
#      rownames(stat) <- NULL
#      stats[[v]] <- stat
#    }
#    traits[[i]] <- bind_rows(stats)
#  }
#  names(traits)[[i]] <- trait_name
#}

#print(traits)




##### add some colors for plotting
# optional, but would ensure uniform colors for populations
# colors to be defined



col2hex <- function(cname) {
  colMat <- col2rgb(cname)
  rgb(red=colMat[1,]/255, green=colMat[2,]/255, blue=colMat[3,]/255)}


col_plot <- data.frame(Population = as.factor(c("AK", "GI", "KA", "MA", "MU", "RE", "UM", "VA", "YE")), Color = col2hex(met.brewer("Johnson", 9)))

droseu <- lapply(droseu, inner_join, col_plot)
#droseu <- lapply(droseu, arrange, Population_Lat)


##### save the data

saveRDS(droseu, file = "droseu_master_list_2022-04-05.rds")







##### quick histogramms



#droseu$via %>%
#  ggplot() +
#  geom_histogram(aes(x = ProportionEggtoAdultSurvival), binwidth = 0.05) +
#  facet_wrap(Supervisor.PI ~ .)



#droseu$dw %>%
#  ggplot() +
#  geom_histogram(aes(x = DW_micrograms), binwidth = 0.02) +
#  facet_grid(Sex ~ Supervisor.PI)


#droseu$wa %>%
#  ggplot() +
#  geom_histogram(aes(x = CentroidSizeLeft_micrometers), binwidth = 50) +
#  facet_grid(Sex ~ Supervisor.PI)











### add missing columns and fill them with NA for easier batch processing

#droseu2 <- droseu
#for (i in (1:length(droseu2))){
#  droseu2[[i]]$Trait <- toupper(names(droseu2)[i])
#  if (!"Line" %in% colnames(droseu2[[i]])) droseu2[[i]]$Line = NA
#  if (!"Batch" %in% colnames(droseu2[[i]])) droseu2[[i]]$Batch = 1
#  if (names(droseu2)[i] %in% c("dia", "fec")) droseu2[[i]]$Sex = "F"
#  if (names(droseu2)[i] %in% c("dtp", "via")) droseu2[[i]]$Sex = NA
#  if (names(droseu2)[i] %in% c("dia", "dw", "fec", "pgm", "pgm2", "tl")) droseu2[[i]]$ReplicateVial = paste(droseu2[[i]]$Batch, droseu2[[i]]$Line, 1, sep = "_")
#  if ("ReplicateCage" %in% colnames(droseu2[[i]])) droseu2[[i]]$ReplicateVial = droseu2[[i]]$ReplicateCage
#  if ("ReplicateChamber" %in% colnames(droseu2[[i]])) droseu2[[i]]$ReplicateVial = droseu2[[i]]$ReplicateChamber
#  if (!"ReplicateVial" %in% colnames(droseu2[[i]])) droseu2[[i]]$ReplicateVial = NA
#}

### keep columns of interest only and split by trait / PI combinations
#common_cols <- c("Trait", "Supervisor.PI", "Diet", "Batch", "Population", "Line", "Sex", "ReplicateVial")
#droseu2 <- lapply(droseu2, function(x) x[, colnames(x) %in% c(common_cols, var_list)])
#droseu2 <- do.call("rbind", droseu2)
#droseu2 <- split(droseu2, f = list(droseu2$Trait, droseu2$Supervisor.PI), drop = T)

#line <- vector("list", length = length(droseu2))
#for (i in 1:length(droseu2)){
#  v <- var_list[var_list %in% colnames(droseu2[[i]])]
#  line[[i]] <- droseu2[[i]] %>%
#    group_by(Supervisor.PI, Sex, Population, Line, ReplicateVial) %>%
#    summarise_at(v, mean) %>%
#    summarise_at(v, mean)
#}
#names(line) <- names(droseu2)



#traits <- vector("list", length = length(line))
#for (i in 1:length(line)) {
#  trait <- line[[i]]
#  trait_name <- names(line)[[i]]
#  trait_var <- colnames(trait)[colnames(trait) %in% var_list]
#  n_var <- length(trait_var)
#  if (n_var > 0) {
#    if ("Sex" %in% colnames(trait)) {
#      trait <- arrange(trait, Supervisor.PI, Sex)
#      lab_names <- unique(paste(trait$Supervisor.PI, trait$Sex, sep = "_"))
#      trait <- group_by(trait, Supervisor.PI, Sex) %>% group_split
#    } else { 
#      trait <- arrange(trait, Supervisor.PI)
#      lab_names <- unique(trait$Supervisor.PI)
#      trait <- group_by(trait, Supervisor.PI) %>% group_split
#    }
#    names(trait) <- lab_names
#    stats <- list()
#    for(v in 1:n_var){
#      l <- unlist(lapply(trait, nrow))
#     if (max(l) < 5000) {
#        test <- lapply(trait, shapiTest, v = trait_var[v])
#        s <- c("W_Statistic", "p")
#      } else {
#        test <- lapply(trait, andersonDarlingTest, v = trait_var[v]) 
#        s <- c("A_Statistic", "p")
#      }
#      stat <- as.data.frame(t(sapply(test, `[`, c("statistic", "p.value"))))
#      colnames(stat) <- s
#      stat$Lab <- lab_names
#      stat$Trait <- trait_name
#      stat$Variable <- trait_var[v]
#      stat$Sig <- ifelse(stat$p < 0.05, "TRUE", "FALSE")
#      stat <- relocate(stat, Trait, Lab, Variable)
#      rownames(stat) <- NULL
#      stats[[v]] <- stat
#    }
#    traits[[i]] <- bind_rows(stats)
#  }
#  names(traits)[[i]] <- trait_name
#}

#traits_line <- bind_rows(traits)




