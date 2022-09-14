
######################################################################################### 
######################  EXTRACT LINES ESTIMATES FROM LINEAR MODELS ######################
#########################################################################################

# We get Line estimates by extracting Line random coefficients from the mixed models using 
# https://github.com/m-clark/mixedup/blob/master/R/extract_random_coefs.R

# Some lab / traits cannot be used here because Line was not included as a random factor in the linear models either because there is no Line replication or Line was removed for a better model fit
# Via_Schmidt_lm_pop no Line replication
# Dia_Flatt_lm_pop no Line replication, used glm instead
# LA_AbsPhase_Tauber_lm_pop singular fit when Line is included
# LS_F_Flatt_lmer_pop phenotype was done on Pop
# LS_M_Flatt_lmer_pop phenotype was done on Pop

# Diapause has been analysed both with LMM and GLMM, the latter allowing to include Line in the model for the Flatt lab data. 



##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(mixedup)
library(foreach)


##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")

##### load data
lmers <- readRDS("LinearModelsPop/all_lmers_pop_list.rds")
length(lmers)
glmers <- readRDS("LinearModelsPop/all_glmers_pop_list.rds")
length(glmers)

##### select models were Line could be included as a random factor, see comments above
models <- c(lmers, glmers)
models <- models[!names(models) %in% c("Via_Schmidt_lm_pop", "Dia_Flatt_lm_pop", "LA_AbsPhase_Tauber_lm_pop", "LS_F_Flatt_lmer_pop", "LS_M_Flatt_lmer_pop")]
length(models)

##### all models at once

### first get random effects with mixedup::extract_random_effects. This is a wrapper for lme4::ranef, gives deviations from each population intercept (not from the global intercept)

line_re <- foreach(m = 1:length(models)) %do% {
  e <- extract_random_effects(models[[m]], re = "Line:Population") %>% 
    separate(group, c("Line", "Population"), ":") %>%
    rename(Value_re = value, SE_re = se) %>%
    dplyr::select(Population, Line, Value_re, SE_re) 
  info <- str_split(names(models)[m], "_") %>% unlist
  if (length(info) >= 5) info[1] <- paste(info[1], info[2], sep = "_")
  info[1] <- sub("_F", "", info[1])
  info[1] <- sub("_M", "", info[1])
  e <- mutate(e, Model = paste(info[length(info)-1], info[length(info)], sep = "_"), 
              Trait = info[1], Lab = info[length(info)-2], 
              Sex = info[length(info)-3]) %>%
    relocate(Model, Trait, Lab, Sex) 
  if (!unique(e$Sex) %in% c("F", "M")) e$Sex <- "NA" 
  if (length(unique(e$Trait)) == 1 & unique(e$Trait) %in% c("Dia", "Fec")) e$Sex <- "F"
  if (length(unique(e$Trait)) == 1 & grepl("Pgm", unique(e$Trait))) e$Sex <- "F"
  if (length(unique(e$Trait)) == 1 & grepl("LA_", unique(e$Trait))) e$Sex <- "B"
  e }
names(line_re) <- names(models)
saveRDS(line_re, "LinearModelsPop/all_models_line_random_effects_list.rds")
write.csv(bind_rows(line_re), file = "LinearModelsPop/all_models_line_random_effects.csv", row.names = F)


### use population estimates (coefs) to get the lines coefs

pops_estimates <- read.csv("LinearModelsPop/all_models_pop_estimates.csv") %>% mutate(Model = paste(Model, Predictor, sep = "_")) %>% dplyr::select(-Predictor)

lines_re <- read.csv("LinearModelsPop/all_models_line_random_effects.csv")

lines_coefs <- inner_join(pops_estimates, lines_re, by = c("Model", "Trait", "Lab", "Sex", "Population")) %>%
  mutate(SE = sqrt(SE^2 + SE_re^2)) %>% 
  mutate(Coef = Estimate + Value_re) %>%
  dplyr::select(Model, Trait, Lab, Sex, Population, Line, Coef, SE)

write.csv(lines_coefs, "LinearModelsPop/all_models_line_random_coefs.csv", row.names = F)






##### by trait

### first get random effects with mixedup::extract_random_effects. This is a wrapper for lme4::ranef, gives deviations from each population intercept (not from the global intercept)

lmers <- list.files(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "_lmers_pop.rds")
glmers <- list.files(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "_glmers_pop.rds")

paths <- c(lmers, glmers)

by_trait <- foreach(p = paths) %do% {
  models <- readRDS(p)
  models <- models[grep("lmer", names(models))]
  models <- models[!names(models) %in% c("LS_F_Flatt_lmer_pop", "LS_M_Flatt_lmer_pop")]
  out_rds <- sub("pop.rds", "line_random_effects.rds", p)
  out_txt <- sub("pop.rds", "line_random_effects.txt", p)
  out_txt_w <- sub("pop.rds", "line_random_effects_wide.txt", p)
  line_re <- foreach(m = 1:length(models), .combine = "rbind") %do% {
    e <- extract_random_effects(models[[m]], re = "Line:Population") %>% 
      separate(group, c("Line", "Population"), ":") %>%
      rename(Value_re = value, SE_re = se) %>%
      dplyr::select(Population, Line, Value_re, SE_re) 
    info <- str_split(names(models)[m], "_") %>% unlist
    if (length(info) >= 5) info[1] <- paste(info[1], info[2], sep = "_")
    info[1] <- sub("_F", "", info[1])
    info[1] <- sub("_M", "", info[1])
    e <- mutate(e, Model = paste(info[length(info)-1], info[length(info)], sep = "_"), 
                Trait = info[1], Lab = info[length(info)-2], 
                Sex = info[length(info)-3]) %>%
      relocate(Model, Trait, Lab, Sex) 
    if (!unique(e$Sex) %in% c("F", "M")) e$Sex <- "NA" 
    if (length(unique(e$Trait)) == 1 & unique(e$Trait) %in% c("Dia", "Fec")) e$Sex <- "F"
    if (length(unique(e$Trait)) == 1 & grepl("Pgm", unique(e$Trait))) e$Sex <- "F"
    if (length(unique(e$Trait)) == 1 & grepl("LA_", unique(e$Trait))) e$Sex <- "B"
    e }
  saveRDS(line_re, file = out_rds)
  write.table(line_re, file = out_txt, row.names = F, quote = F, sep = "\t")
  #line_re_wide <- pivot_wider(line_re, names_from = Lab, values_from = c(Value_re, SE_re))
  #write.table(line_re_wide, file = out_txt_w, row.names = F, quote = F, sep = "\t")
}


### use population estimates (coefs) to get the lines coefs

re <- list.files(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "_line_random_effects.rds")

for (i in re){
  
  trait_line_re <- readRDS(i)
  
  trait_pop_es <- readRDS(sub("_line_random_effects.rds", "_pop_model_estimates.rds", i)) %>% mutate(Model = paste(Model, Predictor, sep = "_")) %>% dplyr::select(-Predictor)
  
  trait_line_coefs <- inner_join(trait_pop_es, trait_line_re, by = c("Model", "Trait", "Lab", "Sex", "Population")) %>%
    mutate(SE = sqrt(SE^2 + SE_re^2)) %>%
    mutate(Coef = Estimate + Value_re) %>%
    dplyr::select(Model, Trait, Lab, Sex, Population, Line, Coef, SE)
  
  saveRDS(trait_line_coefs, file = sub("_effects.rds", "_coefs.rds", i))
  write.table(trait_line_coefs, file = sub("_effects.rds", "_coefs.txt", i), row.names = F, quote = F, sep = "\t")
  #trait_line_coefs_wide <- pivot_wider(trait_line_coefs, names_from = Lab, values_from = c(Coef, SE))
  #write.table(trait_line_coefs_wide, file = sub("_effects.rds", "_coefs_wide.txt", i), row.names = F, quote = F, sep = "\t")
}


   


  