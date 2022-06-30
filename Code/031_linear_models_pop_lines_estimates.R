
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
lmers <- readRDS("LinearModelsPop/all_lmers_pop_list.rds")
glmers <- readRDS("LinearModelsPop/all_glmers_pop_list.rds")


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
saveRDS(line_re, "LinearModelsPop/all_models_line_random_effects_list.rds")
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

   
  
  
####### compound values for line estimates ######

##### source functions
source("Code/functions.R")

library(meta)

models <- list.files(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "lmers_line_random_effects.rds")


##### meta loop

for (f in 1:length(models)) {
  
  m <- readRDS(models[f])
  m_effects <- makeEffects(m)
  
  out_rds <- sub("line_random", "line_compound_random", models[f])
  out_txt <- sub(".rds", ".txt", out_rds)
  
  traits <- unique(m_effects$Trait)
  meta_trait_res <- list()
  meta_trait_out <- list()
  
  for(tr in 1:length(traits)) {
    trait <- filter(m_effects, Trait == traits[tr])
    sexes <- unique(trait$Sex)
    meta_sex_res <- list()
    meta_sex_out <- list()
  
    for (s in 1:length(sexes)) {
      msex <- metagen(data = filter(trait, Sex == sexes[s]), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE, method.tau = "DL") ### change method
      mres <- update.meta(msex, subgroup = Line, tau.common = FALSE)
      mout <- data.frame(Trait = traits[tr], 
                         Population = str_sub(mres$bylevs, 1, 2), 
                         Line = mres$bylevs, Sex = sexes[s], 
                         Value = mres$TE.random.w, 
                         SE = mres$seTE.random.w, 
                         LLM = mres$lower.random.w, 
                         ULM = mres$upper.random.w, 
                         N_lab = mres$k.w)
      meta_sex_res[[s]] <- mres
      mnane <- str_match(out_rds, '([^/]+)(?:/[^/]+){0}$')[,1]
      names(meta_sex_res)[s] <- sub(".rds", paste0("_", traits[tr], "_", sexes[s]), mnane)
      meta_sex_out[[s]] <- mout
    }
    meta_trait_out[[tr]] <- bind_rows(meta_sex_out)
    meta_trait_res[[tr]] <- meta_sex_res
  }
  meta_trait_out <- bind_rows(meta_trait_out)
  meta_trait_res <- unlist(meta_trait_res, recursive=FALSE)
  
  saveRDS(meta_trait_out, out_rds)
  write.table(meta_trait_out, out_txt, row.names = F, quote = F, sep = "\t")
  capture.output(meta_trait_res, file = sub(".txt", "_meta.txt", out_txt)) 
  
}



metas_list <- rdsBatchReaderToList(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "line_compound_random_effects.rds")

metas_list <- metas_list[-5]


lines_wide <- bind_rows(metas_list) %>%
  dplyr::select(Trait, Population, Line, Sex, Value) %>%
  pivot_wider(values_from = Value, names_from = c(Trait, Sex)) 



ggpairs(select(lines_wide, contains("_F")))




lines_wide <- bind_rows(metas_list) %>%
  dplyr::select(Trait, Population, Line, Sex, Value) %>%
  dplyr::group_by(Population, Line, Trait, Sex) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L) 




