
######################################################################
###################### POPULATION META ANALYSES ######################
######################################################################



# CCRT_F_meta.rds contains the output of subgroup meta
# CCRT_F_meta_summary.txt summary of subgroup meta as txt file 
# CCRT_F_meta_pop_compound_estimates.rds compound estimates
# CCRT_F_meta_pop_compound_estimates.txt compound estimates as txt file



##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(meta)
library(ggpubr)
library(MetBrewer)
library(ggrepel)

##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")

##### source functions
source("Code/functions.R")

##### set and create main output directory
meta_dir <- "MetaAnalyses"
dir.create(meta_dir, showWarnings = F) 

##### set input directory containing the linear model outputs
lmer_dir <- "LinearModelsPop"
lmer_pat <- "lmers_pop_model_estimates.rds"



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
pops_coefs <- list.files(path = lmer_dir, recursive = T, full.names = T, pattern = lmer_pat)
print(pops_coefs)

##### remove Dia lmers, keep glmers instead
pops_coefs <- pops_coefs[-grep("Dia_lmers", pops_coefs)]



##### meta loop to get compound estimates per trait

for (f in 1:length(pops_coefs)) {
  # get linear models path
  fpath <- pops_coefs[f]
  # create output directory
  trait_dir <- str_split(fpath, "/", simplify = T)[2]
  trait_dir <- paste(meta_dir, trait_dir, sep = "/")
  dir.create(trait_dir, showWarnings = F) 
  # read linear models output
  m <- readRDS(fpath)
  # remove previoulsy defined partial data
  m_complete <- filter(m, !Lab %in% partial_labs)
  # make effects
  m_effects <- makeEffects(m_complete)
  # get the subtraits
  traits <- unique(m_effects$Trait)
  # loop over subtraits
  for(tr in 1:length(traits)) {
    trait <- filter(m_effects, Trait == traits[tr])
    # get the sexes
    sexes <- unique(trait$Sex)
    # loop over sexes and perform meta
    for (s in 1:length(sexes)) {
      # define output files names
      mod <- ifelse(grepl( "glmers", fpath), "glmers", "lmers") 
      out_sex_rds <- file.path(trait_dir, paste(traits[tr], sexes[s], mod, "pop_meta.rds", sep ="_")) 
      out_sex_txt <- sub(".rds", "_summary.txt", out_sex_rds)
      out_sex_comp_rds <- sub(".rds", "_compound_estimates.rds", out_sex_rds)
      out_sex_comp_txt <- sub(".rds", ".txt", out_sex_comp_rds)
      # meta analysis
      msex <- metagen(data = filter(trait, Sex == sexes[s]), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE, method.tau = "REML")
      # subgroup analysis
      mres <- update.meta(msex, subgroup = Population, tau.common = FALSE)
      # save meta object and get meta summary for each sex
      saveRDS(mres, file = out_sex_rds)
      capture.output(summary(mres), file = out_sex_txt) 
      # extract compound estimates from the meta object to a data.frame
      comp_estimates <- data.frame(Models = mod, 
                                   Trait = traits[tr], 
                                   Sex = sexes[s],
                                   Population = str_sub(mres$bylevs, 1, 2), 
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
all_metas_pop_list <- rdsBatchReaderToList(path = meta_dir, recursive = T, full.names = T, pattern = "lmers_pop_meta_compound_estimates.rds")
saveRDS(all_metas_pop_list, file.path(file = meta_dir, "all_models_pop_meta_compound_estimates_list.rds"))

##### as collapsed list
all_metas_pop <- bind_rows(all_metas_pop_list)
saveRDS(all_metas_pop, file = file.path(meta_dir, "all_models_pop_meta_compound_estimates.rds"))
write.csv(all_metas_pop, file = file.path(meta_dir, "all_models_pop_meta_compound_estimates.csv"), row.names = F)




############# PVALUES AND MULTIPLE TESTING CORRECTION ############# 

##### read in meta results
all_metas_pop_list <- rdsBatchReaderToList(path = meta_dir, recursive = T, full.names = T, pattern = "lmers_pop_meta_compound_estimates.rds")

##### define function to extract p values
getPvalue <- function(x) {
    x %>% mutate(Min_lab = min(N_lab), Max_lab = max(N_lab)) %>% 
    select(Models, Trait, Sex, Q, P, Min_lab, Max_lab) %>%
    distinct()
}

##### apply the function to all meta results
all_metas_pvalues <- lapply(all_metas_pop_list, getPvalue) %>%
  bind_rows()

##### define number of traits and apply Bonferroni correction
#n_traits <- 13 # unique traits we run meta on
#n_traits <- 22 # unique traits we run meta on + sex
#n_traits <- 17 - 1 # unique traits + subtraits we run meta on
n_traits <- 27 -1 # unique traits + subtraits we run meta on + sex
# DW_M is also not meta analysed, hence - 1

all_metas_pvalues_multi <- all_metas_pvalues %>%
  mutate(P_bonf = P * n_traits,
         P_bonf = ifelse(P_bonf > 1, 1, P_bonf))

##### output pvalues
saveRDS(all_metas_pvalues_multi, file.path(meta_dir, "all_models_pop_meta_pvalues.rds"))
write.csv(all_metas_pvalues_multi, file.path(meta_dir, "all_models_pop_meta_pvalues.csv"), row.names = F)


############# PVALUE Q PLOT ############# 

all_metas_pvalues_multi_filt <- all_metas_pvalues_multi %>% 
  filter(Min_lab > 1) 

pvalue_plot <- all_metas_pvalues_multi_filt %>%
  ggplot(aes(x = Q, y = -log10(P))) +
  geom_point(size = 4, alpha = 0.3, pch = 19) +
  theme_bw(16) +
  geom_hline(yintercept = -log10(0.05/nrow(all_metas_pvalues_multi_filt)), linetype = 2, size = 0.5, col = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, size = 0.5) +
  labs(title = "Meta analyses Q and P values", x = "Q value", y = "-log10(Pvalue)") +
  geom_label_repel(data = mutate(all_metas_pvalues_multi_filt, Label = ifelse(P >= 0.05, NA, paste(Trait, Sex))), aes(x = Q, y = -log10(P), label = Label), size = 3, box.padding = 0.45, label.padding = 0.45, point.padding = 0, segment.color = 'grey50', max.overlaps = 10, min.segment.length = 0, seed = 1, force = 20) +
  annotate("text", x = max(all_metas_pvalues_multi_filt$Q), y = -log10(0.05/nrow(all_metas_pvalues_multi_filt)) + 0.2, label = "Bonferroni correction", color = "red", hjust = 1, vjust = 0)

ggsave(pvalue_plot, filename = file.path(meta_dir, "all_models_pop_meta_pvalues.pdf"), height = 5, width = 5)





############# FACET PLOT ############# 


##### list all metas
metas <- readRDS("MetaAnalyses/all_models_pop_meta_compound_estimates.rds")

##### read in population info
pops <- readRDS("InfoTables/DrosEU_Populations.rds")

##### prepare data for plotting
metas_prep <- metas %>% 
  mutate(Group = paste(Trait, Sex, Models, sep = " "),
         Q_plot = paste0("italic(Q) == ", round(Q, 2)),
         P_plot = ifelse(P * n_traits < 0.001, "italic(p) < 0.001", paste0("italic(p) == ", round(P * n_traits, 3))),
         P_plot = ifelse(P * n_traits > 1, "italic(p) == 1", P_plot)) %>%
  group_by(Group) %>% 
  mutate(xpos = min(LLEst), ypos = 11.2, ypos2 = 10.2, N_plot = paste("N ==", round(N_lab_av, 1))) %>%
  filter(N_lab_av >= 2) %>%
  mutate(Population = factor(Population, levels = pops$by_lat$Population),
         y = as.numeric(Population))

##### prepare additional text
stats_text <- select(metas_prep, Group, Q_plot, P_plot, N_plot, xpos, ypos, ypos2) %>%
  distinct()

##### facet plot
p_meta_CI_facet <- ggplot(data = metas_prep, aes(x = Estimate, y = y, color = Population)) +
  facet_wrap(Group ~ ., scales = "free_x", ncol = 5) +
  theme_bw() +
  geom_point(shape = 15, size = 2) +
  geom_errorbarh(aes(xmin = LLEst, xmax = ULEst), height = 0) +
  droseu_color_scale +
  scale_y_continuous(name = "Population", breaks = 1:length(unique(metas_prep$Population)), labels = unique(metas_prep$Population)) +
  labs(x = "Population summary effect", y = "Population", title = "Subgroup meta analyses results", subtitle = "Populations summary effects with 95% CI, Q and P values and average number of labs (N)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text(data = stats_text, aes(x = xpos, y = ypos, label = paste(Q_plot, P_plot, sep = "~~")), color = "black", parse = T, hjust = 0, size = 3, vjust = 1) +
  geom_text(data = stats_text, aes(x = xpos, y = ypos2, label = N_plot), color = "black", parse = T, hjust = 0, size = 3, vjust = 1) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  theme(legend.position = "none") +
  expand_limits(y = 11)

##### save facet plot
ggsave(p_meta_CI_facet, filename = "MetaAnalyses/all_models_pop_meta_summary_effect.pdf", width = 8.27, height = 11.69)
ggsave(p_meta_CI_facet, filename = "MetaAnalyses/all_models_pop_meta_summary_effect.png", width = 8.27, height = 11.69)



############# INDVIDUAL PLOTS ############# 

# re usage of facet plot prep data

##### list all metas
metas <- list.files(path = meta_dir, recursive = T, full.names = T, pattern = "lmers_pop_meta_compound_estimates.rds")

for (i in 1:length(metas)){
  fpath <- metas[i]
  # set output names
  p_out_pdf <- sub("compound_estimates.rds", "summary_effect.pdf", fpath)
  p_out_png <- sub(".pdf", ".png", p_out_pdf)
  # read meta in and format Q and P, order Population by latitude
  m <- readRDS(fpath)
  m <- inner_join(select(m, Models, Trait, Sex, Population), metas_prep)
  if (nrow(m) > 0) {
    title_text <- paste(unique(m$Trait), unique(m$Sex), "summary effect with 95% CI")
    qvalue <- unique(m$Q_plot)
    pvalue <- unique(m$P_plot)
    nvalue <- unique(m$N_plot)

    p_meta_CI <- ggplot(data = m, aes(x = Estimate, y = y, color = Population)) +
      theme_bw(16) +
      geom_point(shape = 15, size = 6) +
      geom_errorbarh(aes(xmin = LLEst, xmax = ULEst), height = 0, lwd = 1.5) +
      droseu_color_scale +
      scale_y_continuous(name = "Population", breaks = 1:length(unique(m$Population)), labels = unique(m$Population)) +
      labs(x = "Population summary effect", y = "Population", title = title_text) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      geom_text(data = m, aes(x = unique(xpos), y = unique(ypos2), label = paste(qvalue, pvalue, nvalue, sep = "~~")), color = "black", parse = T, hjust = 0, size = 5, vjust = 1) +
      theme(legend.position = "none") +
      theme(plot.title = element_text(size = 16)) +
      expand_limits(y = 10)
  
    ggsave(p_meta_CI, filename = p_out_pdf, width = 5, height = 5)
    ggsave(p_meta_CI, filename = p_out_png, width = 5, height = 5)
  }
}


































################### DEPRECATED ####################




##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(meta)
library(ggpubr)
library(MetBrewer)
library(foreach)


##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")


##### load data
droseu <- readRDS("Data/droseu_master_list_2022-05-02.rds")
estimates <- readRDS("LinearModelsPop/all_models_pop_estimates_list.rds")
pops <- readRDS("InfoTables/DrosEU_Populations.rds")


##### create output directory
meta_dir <- "MetaAnalyses"
dir.create(meta_dir, showWarnings = F) 


##### source functions
source("Code/functions.R")



############# VIABILITY #############

# create output directory
out_dir <- "Viability"
dir.create(file.path(meta_dir, out_dir), showWarnings = F) 

# get effects
Via_effects <- makeEffects(estimates$via_lmer)

# run meta
Via_meta <- metagen(data = filter(Via_effects, Sex == "NA"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)
# sub groups
Via_meta <- update.meta(Via_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(Via_meta, file = file.path(meta_dir, out_dir, "Via_meta.rds"))



############# DEVELOPMENT TIME #############

# create output directory
out_dir <- "DevelopmentTime"
dir.create(file.path(meta_dir, out_dir), showWarnings = F) 

# get effects
DT_A_effects <- makeEffects(filter(estimates$dt_lmer, Trait == "DT_A"))

# run meta
DT_A_F_meta <- metagen(data = filter(DT_A_effects, Sex == "F"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)
# sub groups
DT_A_F_meta <- update.meta(DT_A_F_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(DT_A_F_meta, file = file.path(meta_dir, out_dir, "DT_A_F_meta.rds"))


# run meta
DT_A_M_meta <- metagen(data = filter(DT_A_effects, Sex == "M"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)
# sub groups
DT_A_M_meta <- update.meta(DT_A_M_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(DT_A_M_meta, file = file.path(meta_dir, out_dir, "DT_A_M_meta.rds"))





############# DRY WEIGHT #############

# create output directory
out_dir <- "DryWeight"
dir.create(file.path(meta_dir, out_dir), showWarnings = F) 

# get effects
DW_effects <- makeEffects(estimates$dw_lmer)

# run meta
DW_F_meta <- metagen(data = filter(DW_effects, Sex == "F"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)
# sub groups
DW_F_meta <- update.meta(DW_F_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(DW_F_meta, file = file.path(meta_dir, out_dir, "DW_F_meta.rds"))


# run meta
DW_M_meta <- metagen(data = filter(DW_effects, Sex == "M"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)
# sub groups
DW_M_meta <- update.meta(DW_M_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(DW_M_meta, file = file.path(meta_dir, out_dir, "DW_M_meta.rds"))






############# THORAX LENGTH #############

# create output directory
out_dir <- "ThoraxLength"
dir.create(file.path(meta_dir, out_dir), showWarnings = F) 

# get effects
TL_effects <- makeEffects(estimates$tl_lmer)

# run meta
TL_F_meta <- metagen(data = filter(TL_effects, Sex == "F"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)
# sub groups
TL_F_meta <- update.meta(TL_F_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(TL_F_meta, file = file.path(meta_dir, out_dir, "TL_F_meta.rds"))


# run meta
TL_M_meta <- metagen(data = filter(TL_effects, Sex == "M"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)
# sub groups
TL_M_meta <- update.meta(TL_M_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(TL_M_meta, file = file.path(meta_dir, out_dir, "TL_M_meta.rds"))





############# WING AREA #############

# create output directory
out_dir <- "WingArea"
dir.create(file.path(meta_dir, out_dir), showWarnings = F) 

# get effects
WA_effects <- makeEffects(estimates$wa_lmer)

# run meta for both L and R for females
#WA_F_meta <- metagen(data = filter(WA_effects, Sex == "F"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)
# sub groups
#WA_F_meta <- update.meta(WA_F_meta, subgroup = Population, tau.common = FALSE)

# save output
#saveRDS(WA_F_meta, file = file.path(meta_dir, out_dir, "WA_F_meta.rds"))

# run meta for both L and R for males
#WA_M_meta <- metagen(data = filter(WA_effects, Sex == "M"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)
# sub groups
#WA_M_meta <- update.meta(WA_M_meta, subgroup = Population, tau.common = FALSE)

# save output
#saveRDS(WA_M_meta, file = file.path(meta_dir, out_dir, "WA_M_meta.rds"))




# run meta for L for females
WA_L_F_meta <- metagen(data = filter(WA_effects, Sex == "F", Trait == "WA_L"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)
# sub groups
WA_L_F_meta <- update.meta(WA_L_F_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(WA_L_F_meta, file = file.path(meta_dir, out_dir, "WA_L_F_meta.rds"))

# run meta for L for males
WA_L_M_meta <- metagen(data = filter(WA_effects, Sex == "M", Trait == "WA_L"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)
# sub groups
WA_L_M_meta <- update.meta(WA_L_M_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(WA_L_M_meta, file = file.path(meta_dir, out_dir, "WA_L_M_meta.rds"))


# run meta for R for females
WA_R_F_meta <- metagen(data = filter(WA_effects, Sex == "F", Trait == "WA_R"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)
# sub groups
WA_R_F_meta <- update.meta(WA_R_F_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(WA_R_F_meta, file = file.path(meta_dir, out_dir, "WA_R_F_meta.rds"))


# run meta for R for males
WA_R_M_meta <- metagen(data = filter(WA_effects, Sex == "M", Trait == "WA_R"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)
# sub groups
WA_R_M_meta <- update.meta(WA_R_M_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(WA_R_M_meta, file = file.path(meta_dir, out_dir, "WA_R_M_meta.rds"))











############# FECUNDITY #############

# create output directory
out_dir <- "Fecundity"
dir.create(file.path(meta_dir, out_dir), showWarnings = F) 

# get effects
Fec_effects <- makeEffects(estimates$fec_lmer)

# run meta
Fec_meta <- metagen(data = filter(Fec_effects, Sex == "F"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)

# sub groups
Fec_meta <- update.meta(Fec_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(Fec_meta, file = file.path(meta_dir, out_dir, "Fec_meta.rds"))





############# LIFESPAN #############

# create output directory
out_dir <- "Lifespan"
dir.create(file.path(meta_dir, out_dir), showWarnings = F) 

# get effects
LS_effects <- makeEffects(estimates$ls_lmer)

# run meta for females
LS_F_meta <- metagen(data = filter(LS_effects, Sex == "F"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)

# sub groups
LS_F_meta <- update.meta(LS_F_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(LS_F_meta, file = file.path(meta_dir, out_dir, "LS_F_meta.rds"))

# run meta for males
LS_M_meta <- metagen(data = filter(LS_effects, Sex == "M"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)

# sub groups
LS_M_meta <- update.meta(LS_M_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(LS_M_meta, file = file.path(meta_dir, out_dir, "LS_M_meta.rds"))





############# COLD-SHOCK MORTALITY #############

# create output directory
out_dir <- "ColdShock"
dir.create(file.path(meta_dir, out_dir), showWarnings = F) 

# get effects
CSM_effects <- makeEffects(estimates$csm_lmer)

# run meta for females
CSM_F_meta <- metagen(data = filter(CSM_effects, Sex == "F"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)

# sub groups
CSM_F_meta <- update.meta(CSM_F_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(CSM_F_meta, file = file.path(meta_dir, out_dir, "CSM_F_meta.rds"))

# run meta for males
CSM_M_meta <- metagen(data = filter(CSM_effects, Sex == "M"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)

# sub groups
CSM_M_meta <- update.meta(CSM_M_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(CSM_M_meta, file = file.path(meta_dir, out_dir, "CSM_M_meta.rds"))






############# CHILL-COMA MORTALITY #############

# create output directory
out_dir <- "ChillComa"
dir.create(file.path(meta_dir, out_dir), showWarnings = F) 

# get effects
CCRT_effects <- makeEffects(estimates$ccrt_lmer)

# run meta for females
CCRT_F_meta <- metagen(data = filter(CCRT_effects, Sex == "F"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)

# sub groups
CCRT_F_meta <- update.meta(CCRT_F_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(CCRT_F_meta, file = file.path(meta_dir, out_dir, "CCRT_F_meta.rds"))

# run meta for males
CCRT_M_meta <- metagen(data = filter(CCRT_effects, Sex == "M"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)

# sub groups
CCRT_M_meta <- update.meta(CCRT_M_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(CCRT_M_meta, file = file.path(meta_dir, out_dir, "CCRT_M_meta.rds"))








############# HEAT-SHOCK MORTALITY #############

# create output directory
out_dir <- "HeatShock"
dir.create(file.path(meta_dir, out_dir), showWarnings = F) 

# get effects
HSM_effects <- makeEffects(estimates$hsm_lmer)

# run meta for females
HSM_F_meta <- metagen(data = filter(HSM_effects, Sex == "F"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)

# sub groups
HSM_F_meta <- update.meta(HSM_F_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(HSM_F_meta, file = file.path(meta_dir, out_dir, "HSM_F_meta.rds"))

# run meta for males
HSM_M_meta <- metagen(data = filter(HSM_effects, Sex == "M"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)

# sub groups
HSM_M_meta <- update.meta(HSM_M_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(HSM_M_meta, file = file.path(meta_dir, out_dir, "HSM_M_meta.rds"))









############# DIAPAUSE #############

# create output directory
out_dir <- "Diapause"
dir.create(file.path(meta_dir, out_dir), showWarnings = F) 

# get effects
Dia_lmer_effects <- makeEffects(estimates$dia_lmer)

# run meta for females
Dia_lmer_meta <- metagen(data = filter(Dia_lmer_effects, Sex == "F"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)

# sub groups
Dia_lmer_meta <- update.meta(Dia_lmer_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(Dia_lmer_meta, file = file.path(meta_dir, out_dir, "Dia_lmer_meta.rds"))


# get effects
Dia_glmer_effects <- makeEffects(estimates$dia_glmer)

# run meta for females
Dia_glmer_meta <- metagen(data = filter(Dia_glmer_effects, Sex == "F"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)

# sub groups
Dia_glmer_meta <- update.meta(Dia_glmer_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(Dia_glmer_meta, file = file.path(meta_dir, out_dir, "Dia_glmer_meta.rds"))










############# STARVATION RESISTANCE #############

# create output directory
out_dir <- "Starvation"
dir.create(file.path(meta_dir, out_dir), showWarnings = F) 

# get effects
SR_effects <- makeEffects(estimates$sr_lmer)

# run meta for females
SR_F_meta <- metagen(data = filter(SR_effects, Sex == "F"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)

# sub groups
SR_F_meta <- update.meta(SR_F_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(SR_F_meta, file = file.path(meta_dir, out_dir, "SR_F_meta.rds"))

# run meta for males
SR_M_meta <- metagen(data = filter(SR_effects, Sex == "M"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)

# sub groups
SR_M_meta <- update.meta(SR_M_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(SR_M_meta, file = file.path(meta_dir, out_dir, "SR_M_meta.rds"))




############# PIGMENTATION #############

# create output directory
out_dir <- "Pigmentation"
dir.create(file.path(meta_dir, out_dir), showWarnings = F) 

# get effects
Pgm_effects <- makeEffects(estimates$pgm_lmer)



# run meta for females - T4
Pgm_T4_meta <- metagen(data = filter(Pgm_effects, Sex == "F", Trait == "Pgm_T4"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)

# sub groups
Pgm_T4_meta <- update.meta(Pgm_T4_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(Pgm_T4_meta, file = file.path(meta_dir, out_dir, "Pgm_T4_meta.rds"))





# run meta for females - T5
Pgm_T5_meta <- metagen(data = filter(Pgm_effects, Sex == "F", Trait == "Pgm_T5"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)

# sub groups
Pgm_T5_meta <- update.meta(Pgm_T5_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(Pgm_T5_meta, file = file.path(meta_dir, out_dir, "Pgm_T5_meta.rds"))




# run meta for females - T6
Pgm_T6_meta <- metagen(data = filter(Pgm_effects, Sex == "F", Trait == "Pgm_T6"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)

# sub groups
Pgm_T6_meta <- update.meta(Pgm_T6_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(Pgm_T6_meta, file = file.path(meta_dir, out_dir, "Pgm_T6_meta.rds"))




# run meta for females - Total
Pgm_Total_meta <- metagen(data = filter(Pgm_effects, Sex == "F", Trait == "Pgm_Total"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)

# sub groups
Pgm_Total_meta <- update.meta(Pgm_Total_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(Pgm_Total_meta, file = file.path(meta_dir, out_dir, "Pgm_Total_meta.rds"))


# run meta for females, w/o Schmidt's data
#Pgm_Total_meta <- metagen(data = filter(Pgm_effects, Sex == "F", Trait == "Pgm_Total",  Lab != "Schmidt"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)

# sub groups
#Pgm_Total_meta <- update.meta(Pgm_Total_meta, subgroup = Population, tau.common = FALSE)


# run meta for females, all three labs
#Pgm_meta <- metagen(data = filter(Pgm_effects, Sex == "F", Trait == "Pgm_Total"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)

# sub groups
#Pgm_meta <- update.meta(Pgm_meta, subgroup = Population, tau.common = FALSE)

# save output
#saveRDS(Pgm_meta, file = file.path(meta_dir, out_dir, "SR_F_meta.rds"))




############# output all metas ############# 

all_metas <- ls()[grep("_meta", ls())]
save(all_metas, file = file.path(meta_dir, "all_metas_list.Rdata"))




######### summaries

metas <- list.files(path = "MetaAnalyses", recursive = T, full.names = T, pattern = "meta.rds")


for (i in 1:length(metas)){
  f <- metas[i]
  s_out_txt <- sub(".rds", "_summary.txt", f)
  m <- readRDS(f)
  s <- summary(m)
  capture.output(s, file = s_out_txt)
}


######### multiple testing correction

rdsBatchReaderToList <- function(...) {
  temp = list.files(...)
  tnames <- str_split(temp, "/", simplify = T)[,3]
  tnames <- str_replace(tnames, ".rds", "")
  tlist <- lapply(temp, readRDS)
  names(tlist) <- tnames
  return(tlist)
}

metas_list <- rdsBatchReaderToList(path = "MetaAnalyses", recursive = T, full.names = T, pattern = "meta.rds")



combineMetaPValues <- function(x) {
  pvals <- list()
  for (meta in 1:length(x)) {
    info <- str_split(names(x)[meta], "_") %>% unlist
    if (length(info) >= 3) info[1] <- paste(info[1], info[2], sep = "_")
    info[1] <- sub("_F", "", info[1])
    info[1] <- sub("_M", "", info[1])
    info[1] <- sub("_lmer", "", info[1])
    info[1] <- sub("_glmer", "", info[1])
    info <- data.frame(Trait = info[1], 
                       Sex = info[length(info)-1],
                       Model = ifelse(grepl("glmer", names(x)[meta]), "glmer", "lmer"),
                       Q = x[[meta]]$Q.b.random,
                       P = x[[meta]]$pval.Q.b.random)
    if (!info$Sex %in% c("F", "M")) info$Sex <- "NA" 
    if (info$Trait %in% c("Dia", "Fec")) info$Sex <- "F"
    if (grepl("Pgm_", info$Trait)) info$Sex <- "F"
    if (grepl("LA_", info$Trait)) info$Sex <- "B" 
    pvals[[meta]] <- info }
  bind_rows(pvals) }

n_traits <- 20

metas_pvalues <- combineMetaPValues(metas_list) %>%
  mutate(P_bonf = P * n_traits,
         P_bonf = ifelse(P_bonf > 1, 1, P_bonf))

write.csv(metas_pvalues, "MetaAnalyses/all_metas_pvalues.csv", row.names = F)
saveRDS(metas_pvalues, "MetaAnalyses/all_metas_pvalues.rds")




#metas_pvalues <- bind_rows(Meta = names(metas_list), Q = lapply(metas_list, function(x) x$Q.b.random) %>% unlist(), P = lapply(metas_list, function(x) x$pval.Q.b.random) %>% unlist()) %>% mutate(P_bonf = P * n_traits) %>% mutate(P_bonf = ifelse(P_bonf > 1, 1, P_bonf))



######### plots

myColors <- met.brewer("Johnson", 9)
names(myColors) <- as.factor(c("AK", "GI", "KA", "MA", "MU", "RE", "UM", "VA", "YE"))
colScale <- scale_colour_manual(name = "Population", values = myColors)

metas <- list.files(path = "MetaAnalyses", recursive = T, full.names = T, pattern = "meta.rds")

for (i in 1:length(metas)){
  f <- metas[i]
  p_out_pdf <- sub(".rds", "_summary_effect.pdf", f)
  p_out_png <- sub(".rds", "_summary_effect.png", f)
  m <- readRDS(f)
  m <- data.frame(Population = m$bylevs, Mstar = m$TE.random.w, SEMstar = m$seTE.random.w, Q = m$Q.b.random, p = m$pval.Q.b.random, LLMstar = m$lower.random.w, ULMstar = m$upper.random.w ) %>% mutate(Q_plot = paste0("italic(Q) == ", round(Q, 2)), P_plot = ifelse(metas_pvalues$P_bonf[i] < 0.001, "italic(p) < 0.001", paste0("italic(p) == ", round(metas_pvalues$P_bonf[i], 3))), Population = factor(Population, levels = pops$by_lat$Population))
  trait_name <- str_match(f, '([^/]+)(?:/[^/]+){0}$')[,1]
  trait_name <- sub("_meta.rds", "", trait_name)
  
  #ann <- data.frame(x = c(-Inf, -Inf), y = c(8, 9), l = c(unique(m$P_plot), unique(m$Q_plot)))
  p_meta_SE <- ggplot(data = m, aes(x = Mstar, y = 1:length(Population), color = Population)) +
    theme_bw() +
    geom_point(size = 8, shape = 15) +
    geom_errorbarh(aes(xmax = Mstar + SEMstar, xmin = Mstar - SEMstar), height = 0) +
    colScale +
    scale_y_continuous(name = "Population", breaks = 1:length(m$Population), labels = m$Population) +
    labs(x = "Summary effect", title = paste(trait_name, "meta with SE")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    annotate("text", x = -Inf, y = Inf, label = unique(m$Q_plot), hjust=-0.2, vjust=3.2, parse = T) +
    annotate("text", x = -Inf, y = Inf, label = unique(m$P_plot), hjust=-0.4, vjust=4.2, parse = T) +
    theme(legend.position = "none")
  
  p_meta_CI <- ggplot(data = m, aes(x = Mstar, y = 1:length(Population), color = Population)) +
    theme_bw() +
    geom_point(size = 8, shape = 15) +
    geom_errorbarh(aes(xmax = ULMstar, xmin = LLMstar), height = 0) +
    colScale +
    scale_y_continuous(name = "Population", breaks = 1:length(m$Population), labels = m$Population) +
    labs(x = "Summary effect", title = paste(trait_name, "meta with 95% CI")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    annotate("text", x = -Inf, y = Inf, label = unique(m$Q_plot), hjust=-0.2, vjust=3.2, parse = T) +
    annotate("text", x = -Inf, y = Inf, label = unique(m$P_plot), hjust=-0.4, vjust=4.2, parse = T) +
    theme(legend.position = "none")
  
  p_meta <- ggarrange(p_meta_SE, p_meta_CI)
  
  pdf(p_out_pdf, width = 8, height = 5)
  print(p_meta)
  dev.off()
  
  png(filename = p_out_png, height = 1500, width = 2400, res = 300)
  print(p_meta)
  dev.off()
}
  

#####



metas <- list.files(path = "MetaAnalyses", recursive = T, full.names = T, pattern = "meta.rds")

for (i in 1:length(metas)){
  f <- metas[i]
  p_out_rds <- sub(".rds", "_pop_compound_estimates.rds", f)
  p_out_txt <- sub(".rds", "_pop_compound_estimates.txt", f)
  m <- readRDS(f)
  m <- data.frame(Population = m$bylevs, Mstar = m$TE.random.w, SEMstar = m$seTE.random.w, Q = m$Q.b.random, p = m$pval.Q.b.random, LLMstar = m$lower.random.w, ULMstar = m$upper.random.w, N_lab = m$k.w)
  saveRDS(m, p_out_rds)
  write.table(m, p_out_txt, row.names = F, quote = F, sep = "\t")
}  



all_metas_pop <- rdsBatchReaderToList(path = "MetaAnalyses", recursive = T, full.names = T, pattern = "_pop_compound_estimates.rds")

# do not consider lmer for diapause
#all_metas_pop <- all_metas_pop[-grep("Dia_lmer", names(all_metas_pop))]


all_metas_pop_comp <- foreach(m = 1:length(all_metas_pop)) %do% {
  info <- str_split(names(all_metas_pop)[m], "_") %>% unlist
  if (length(info) >= 7 | info[1] %in% c("Dia", "Pgm")) info[1] <- paste(info[1], info[2], sep = "_")
  info[1] <- sub("_F", "", info[1])
  info[1] <- sub("_M", "", info[1])
  e <- mutate(all_metas_pop[[m]], Trait = info[1], Sex = info[length(info)-4]) %>%
    relocate(Trait, Sex) 
  if (!unique(e$Sex) %in% c("F", "M")) e$Sex <- "NA" 
  if (length(unique(e$Trait)) == 1 & unique(e$Trait) %in% c("Dia", "Fec")) e$Sex <- "F"
  if (length(unique(e$Trait)) == 1 & grepl("Pgm", unique(e$Trait))) e$Sex <- "F"
  if (length(unique(e$Trait)) == 1 & grepl("LA_", unique(e$Trait))) e$Sex <- "B"
  e }
names(all_metas_pop_comp) <- names(all_metas_pop)
saveRDS(all_metas_pop_comp, "MetaAnalyses/all_metas_pop_compound_estimates_list.rds")
write.csv(bind_rows(all_metas_pop_comp), file = "MetaAnalyses/all_metas_pop_compound_estimates.csv", row.names = F)



#str(m , list.len = length(m))



####### meta regression
# just a trial to see how we could use meta regression to take into account the effect of diet

library(metafor)

diets <- read.csv("InfoTables/DrosEU_Diets_Feb22.csv") %>%
  select(PI.Lab.head, Diet, P.C) %>%
  rename(Lab = PI.Lab.head) %>%
  mutate(Lab = ifelse(Lab == "Stamenkovic-Radak", "StamenkovicRadak", Lab)) %>%
  distinct()


### viability

# get effects
Via_effects <- makeEffects(estimates$via_lmer)

# add diet info
Via_effects <- inner_join(Via_effects, diets)


# run meta
Via_meta <- metagen(data = filter(Via_effects, Sex == "NA"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE, method.tau = "ML")

# this gives the same as the subgroup analysis
Via_meta_reg <- metareg(Via_meta, ~ Population, method.tau = "REML")


Via_meta_reg <- rma(yi = Y, sei = SE, data = Via_effects %>% filter(Lab != "Test"), method = "ML", mods = ~ Population)

Via_meta_reg_full <- rma(yi = Y, sei = SE, data = Via_effects, method = "ML", mods = ~ Population + Lab)

anova(Via_meta_reg_full, Via_meta_reg)

Via_meta_reg_int <- rma(yi = Y, sei = SE, data = Via_effects, method = "ML", mods = ~ Population * P.C)





############# WING AREA #############


# get effects
WA_effects <- makeEffects(estimates$wa_lmer)

# add diet info
WA_effects <- inner_join(WA_effects, diets)


WA_L_F_meta <- metagen(data = filter(WA_effects, Lab != "Posnien", Sex == "F", Trait == "WA_L"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE, method.tau = "ML")
# sub groups
WA_L_F_meta <- update.meta(WA_L_F_meta, subgroup = Population, tau.common = FALSE)



WA_L_F_meta_reg <- rma(yi = Y, sei = SE, data = WA_effects %>% filter(Lab != "Posnien", Sex == "F", Trait == "WA_L"), method = "ML", mods = ~ Population)

WA_L_F_meta_reg_full <- rma(yi = Y, sei = SE, data = WA_effects %>% filter(Lab != "Posnien", Sex == "F", Trait == "WA_L"), method = "ML", mods = ~ Population + P.C)

anova(WA_L_F_meta_reg, WA_L_F_meta_reg_full)




############# THORAX LENGTH #############


# get effects
TL_effects <- makeEffects(estimates$tl_lmer)

# add diet info
TL_effects <- inner_join(TL_effects, diets)


TL_F_meta <- metagen(data = filter(TL_effects, Lab != "Posnien", Sex == "F"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE, 
                     method.tau = "ML")
# sub groups
TL_F_meta <- update.meta(TL_F_meta, subgroup = Population, tau.common = FALSE)



TF_F_meta_reg <- rma(yi = Y, sei = SE, data = TL_effects %>% filter(Lab != "Posnien", Sex == "F"), method = "ML", mods = ~ Population)

TF_F_meta_reg_full <- rma(yi = Y, sei = SE, data = TL_effects %>% filter(Lab != "Posnien", Sex == "F"), method = "ML", mods = ~ Population +  (1|P.C))

anova(TF_F_meta_reg, TF_F_meta_reg_full)





m <- data.frame(Population = WA_L_F_meta$bylevs, Mstar = WA_L_F_meta$TE.random.w, SEMstar = WA_L_F_meta$seTE.random.w, Q = WA_L_F_meta$Q.b.random, p = WA_L_F_meta$pval.Q.b.random, LLMstar = WA_L_F_meta$lower.random.w, ULMstar = WA_L_F_meta$upper.random.w) %>% mutate(Q_plot = paste0("italic(Q) == ", round(Q, 2)), P_plot = ifelse(p < 0.001, "italic(p) < 0.001", paste0("italic(p) == ", round(p, 3))), Population = factor(Population, levels = pops$by_lat$Population))

p_meta_SE <- ggplot(data = m, aes(x = Mstar, y = 1:length(Population), color = Population)) +
  theme_bw() +
  geom_point(size = 8, shape = 15) +
  geom_errorbarh(aes(xmax = Mstar + SEMstar, xmin = Mstar - SEMstar), height = 0) +
  colScale +
  scale_y_continuous(name = "Population", breaks = 1:length(m$Population), labels = m$Population) +
  labs(x = "Summary effect") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  annotate("text", x = -Inf, y = Inf, label = unique(m$Q_plot), hjust=-0.2, vjust=3.2, parse = T) +
  annotate("text", x = -Inf, y = Inf, label = unique(m$P_plot), hjust=-0.4, vjust=4.2, parse = T) +
  theme(legend.position = "none")

p_meta_CI <- ggplot(data = m, aes(x = Mstar, y = 1:length(Population), color = Population)) +
  theme_bw() +
  geom_point(size = 8, shape = 15) +
  geom_errorbarh(aes(xmax = ULMstar, xmin = LLMstar), height = 0) +
  colScale +
  scale_y_continuous(name = "Population", breaks = 1:length(m$Population), labels = m$Population) +
  labs(x = "Summary effect") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  annotate("text", x = -Inf, y = Inf, label = unique(m$Q_plot), hjust=-0.2, vjust=3.2, parse = T) +
  annotate("text", x = -Inf, y = Inf, label = unique(m$P_plot), hjust=-0.4, vjust=4.2, parse = T) +
  theme(legend.position = "none")

p_meta <- ggarrange(p_meta_SE, p_meta_CI)

pdf("MetaAnalyses/WingArea/WA_L_F_meta_summary_effect_wo_partial_data.pdf", width = 8, height = 5)
print(p_meta)
dev.off()


ggplot(data = Via_effects, aes(x = P.C, y = Y)) + geom_point(aes(color = Population))
ggplot(data = Via_effects, aes(x = Lab, y = Y)) + geom_point(aes(color = Population))
ggplot(data = Via_effects, aes(x = Lab, y = Y)) + geom_boxplot()


ggplot(data = TL_effects, aes(x = P.C, y = Y)) + geom_point(aes(color = Population))
ggplot(data = TL_effects, aes(x = Lab, y = Y)) + geom_point(aes(color = Population))
ggplot(data = TL_effects, aes(x = Lab, y = Y)) + geom_boxplot()










