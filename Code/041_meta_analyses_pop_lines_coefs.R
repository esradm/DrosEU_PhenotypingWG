
#################################################################################
################# META ANALYSES TO GET COMPOUND LINE ESTIMATES ##################
#################################################################################

 




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

##### set and create main output directory
meta_dir <- "MetaAnalyses"
dir.create(meta_dir, showWarnings = F) 

##### set input directory containing the linear model outputs
lmer_dir <- "LinearModelsPop"
lmer_pat <- "lmers_line_random_coefs.rds"



############# PARTIAL DATA ############# 

##### define data that is considered partial and that should be removed from meta analyses

# Posnien in WA and TL 3 lines per pop
# Ritchie in WA and TL 5 lines per pop, we lose TL Males if removed
# Grath in Via and DT_A 10 lines for 3 pops

# would remove Posnien for sure, Ritchie and Grath to be discussed

#partial_labs <- c("Posnien", "Grath", "Ritchie")
partial_labs <- "Posnien"



############# SUBGROUP META ANALYSES ############# 

##### get all the line estimates, from all the traits
lines_coefs <- list.files(path = lmer_dir, recursive = T, full.names = T, pattern = lmer_pat)
print(lines_coefs)

##### remove Dia lmers, keep glmers instead
lines_coefs <- lines_coefs[-grep("Dia_lmers", lines_coefs)]


##### meta loop to get compound estimates per trait

for (f in 1:length(lines_coefs)) {
  # get linear models path
  fpath <- lines_coefs[f]
  # create output directory
  trait_dir <- str_split(fpath, "/", simplify = T)[2]
  trait_dir <- paste(meta_dir, trait_dir, sep = "/")
  dir.create(trait_dir, showWarnings = F) 
  # read linear models output
  lc <- readRDS(fpath) %>% rename(Estimate = Coef)
  # remove previoulsy defined partial data
  lc_complete <- filter(lc, !Lab %in% partial_labs)
  # make effects
  lc_effects <- makeEffects(lc_complete)
  # get the subtraits
  traits <- unique(lc_effects$Trait)
  # loop over subtraits
  for(tr in 1:length(traits)) {
    trait <- filter(lc_effects, Trait == traits[tr])
    # get the sexes
    sexes <- unique(trait$Sex)
    # loop over sexes and perform meta
    for (s in 1:length(sexes)) {
      # define output files names
      mod <- ifelse(grepl( "glmers", fpath), "glmers", "lmers") 
      out_sex_rds <- file.path(trait_dir, paste(traits[tr], sexes[s], mod, "line_meta.rds", sep ="_")) 
      out_sex_txt <- sub(".rds", "_summary.txt", out_sex_rds)
      out_sex_comp_rds <- sub(".rds", "_compound_random_coefs.rds", out_sex_rds)
      out_sex_comp_txt <- sub(".rds", ".txt", out_sex_comp_rds)
      # meta analysis
      msex <- metagen(data = filter(trait, Sex == sexes[s]), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE, method.tau = "REML")
      # subgroup analysis
      mres <- update.meta(msex, subgroup = Line, tau.common = FALSE)
      # save meta object and get meta summary for each sex
      saveRDS(mres, file = out_sex_rds)
      capture.output(summary(mres), file = out_sex_txt) 
      # extract compound estimates from the meta object to a data.frame
      comp_estimates <- data.frame(Models = mod, 
                                   Trait = traits[tr], 
                                   Sex = sexes[s],
                                   Population = str_sub(mres$bylevs, 1, 2),
                                   Line = mres$bylevs,
                                   Estimate = mres$TE.random.w, 
                                   SE = mres$seTE.random.w, 
                                   LLEst = mres$lower.random.w, 
                                   ULEst = mres$upper.random.w,
                                   Q = mres$Q.b.random, 
                                   P = mres$pval.Q.b.random,
                                   N_lab = mres$k.w,
                                   N_lab_av = mean(mres$k.w))
      # save compound estimates object as rds and txt
      saveRDS(comp_estimates, file = out_sex_comp_rds)
      write.table(comp_estimates, out_sex_comp_txt, row.names = F, quote = F, sep = "\t")
    }
  }
}




############# COMBINE META ANALYSES OUTPUTS ############# 

##### as a global list
all_metas_line_list <- rdsBatchReaderToList(path = meta_dir, recursive = T, full.names = T, pattern = "lmers_line_meta_compound_random_coefs.rds")
saveRDS(all_metas_line_list, file.path(file = meta_dir, "all_models_line_meta_compound_random_coefs_list.rds"))

##### as collapsed list
all_metas_line <- bind_rows(all_metas_line_list)
saveRDS(all_metas_line, file = file.path(meta_dir, "all_models_line_meta_compound_random_coefs.rds"))
write.csv(all_metas_line, file = file.path(meta_dir, "all_models_line_meta_compound_random_coefs.csv"), row.names = F)

##### in wide format
all_metas_line_wide <- all_metas_line %>%
  dplyr::select(Trait, Population, Line, Sex, Estimate) %>%
  pivot_wider(values_from = Estimate, names_from = c(Trait, Sex)) 
write.csv(all_metas_line_wide, file = file.path(meta_dir, "all_models_line_meta_compound_random_coefs_wide.csv"), row.names = F)





##### quick plots for trait correlations

#library(GGally)
#corF <- ggpairs(select(compound_effects_wide, contains(c("_F", "NA"))))
#ggsave(corF, file = "~/Desktop/trait_cor_F.pdf", width = 14, height = 14)

#corM <- ggpairs(select(compound_effects_wide, contains(c("_M", "NA"))))
#ggsave(corM, file = "~/Desktop/trait_cor_M.pdf", width = 14, height = 14)










