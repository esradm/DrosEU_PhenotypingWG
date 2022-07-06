
#######################################################################################
######################  RUN ALL META ANALYSES FOR LINE ESTIMATES ######################
#######################################################################################


##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(meta)
library(ggpubr)
library(MetBrewer)


##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")

##### source functions
source("Code/functions.R")


##### get all the line estimates, from all the traits
lines_estimates <- list.files(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "lmers_line_random_effects.rds")


##### meta loop to get compound estimates per trait
# at the moment it is quick and dirty because it also runs a meta on traits for which only one lab is involved, returning the original line estomates as compound estimates. This is correct but maybe we should do that in a better way.

for (f in 1:length(lines_estimates)) {
  
  m <- readRDS(lines_estimates[f])
  m_effects <- makeEffects(m)
  
  out_rds <- sub("line_random", "line_compound_random", lines_estimates[f])
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


##### combine traits

compound_effects <- rdsBatchReaderToList(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "line_compound_random_effects.rds")

compound_effects <- compound_effects[-5] # remove diapause lmer

compound_effects_wide <- bind_rows(compound_effects) %>%
  dplyr::select(Trait, Population, Line, Sex, Value) %>%
  pivot_wider(values_from = Value, names_from = c(Trait, Sex)) 

saveRDS(compound_effects, "LinearModelsPop/all_models_line_compound_random_effects_list.rds")
write.csv(compound_effects_wide, file = "LinearModelsPop/all_models_line_compound_random_effects.csv", row.names = F)


##### quick plot

library(GGally)
ggpairs(select(compound_effects_wide, contains("_F")))






