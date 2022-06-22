
############################################################################### 
######################  EXTRACT LINES ESTIMATES FROM LMM ######################
###############################################################################




##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
#library(lme4)
#library(lsmeans)
#library(afex)
#library(multcomp)
library(mixedup)
library(foreach)


##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")

##### load data
lmers <- readRDS("LinearModelsPop/all_lmers_list_pop.rds")
glmers <- readRDS("LinearModelsPop/all_glmers_list_pop.rds")


##### select models were Line could be included as a random factor

models <- c(lmers, glmers)
models <- models[grep("lmer", names(models))] # removes 3 lm
models <- models[!names(models) %in% c("LS_F_Flatt_lmer_pop", "LS_M_Flatt_lmer_pop")]


##### all models at once

line_re <- foreach(m = 1:length(models)) %do% {
  e <- extract_random_effects(models[[m]], re = "Line:Population") %>% 
    separate(group, c("Line", "Population"), ":") %>%
    rename(Estimate = value, SE = se) %>%
    dplyr::select(Population, Line, Estimate, SE) 
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
saveRDS(line_re, "LinearModelsPop/all_models_list_line_random_effects.rds")
write.csv(bind_rows(line_re), file = "LinearModelsPop/all_models_line_random_effects.csv", row.names = F)


##### by trait

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
      rename(Estimate = value, SE = se) %>%
      dplyr::select(Population, Line, Estimate, SE) 
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
  line_re_wide <- pivot_wider(line_re, names_from = Lab, values_from = c(Estimate, SE))
  write.table(line_re_wide, file = out_txt_w, row.names = F, quote = F, sep = "\t")
}

   
  
  
  
  
  
  
  


#lines_re <- list()
#for (m  in 1:length(models)) {
#  e <- extract_random_effects(models[[m]], re = "Population:Line") %>% 
#    mutate(Model = sub("_lmer_pop", "", names(models)[m])) %>%
#    separate(group, c("Line", "Population"), ":") %>%
#    rename(Estimate = value, SE = se) %>%
#    dplyr::select(Model, Population, Line, Estimate, SE) 
#  info <- str_split(unique(e$Model), "_") %>% unlist
#  if (length(info) >= 3) info[1] <- paste(info[1], info[2], sep = "_")
#  info[1] <- sub("_F", "", info[1])
#  info[1] <- sub("_M", "", info[1])
#  e <- mutate(e, Trait = info[1], Lab = info[length(info)], Sex = info[length(info)-1]) %>%
#    relocate(Model, Trait, Lab, Sex) 
#  if (!unique(e$Sex) %in% c("F", "M")) e$Sex <- "NA" 
#  if (length(unique(e$Trait)) == 1 & unique(e$Trait) %in% c("Dia", "Fec")) e$Sex <- "F"
#  if (length(unique(e$Trait)) == 1 & grepl("Pgm", unique(e$Trait))) e$Sex <- "F"
#  if (length(unique(e$Trait)) == 1 & unique(e$Trait) %in% c("LA")) e$Sex <- "B"
#  lines_re[[m]] <- e
#  }
  

wl_lme4 <- lme4::lmer(CentroidSizeLeft_micrometers ~ Population + (1|Line), data = filter(droseu$wa, Supervisor.PI == "Onder"))

wl_re <- extract_random_effects(wl_lme4, re = "Line")

head(wl_re)

wl_coef <- extract_random_coefs(wl_lme4, re = "Line")

head(wl_coef)

cor.test(wl_re$value, wl_coef$value)





wa_afex <- afex::lmer(CentroidSizeLeft_micrometers ~ Population + (1|Line), data = filter(droseu$wa, Supervisor.PI == "Onder"))

extract_random_effects(w, re = "Line")
extract_random_coefs(w, re = "Line")


g <- readRDS("LinearModelsPop/Diapause/Dia_glmers_pop.rds")

extract_random_effects(g[[1]], re = "Line:Population")
extract_random_coefs(g[[1]], re = "Line:Population")





