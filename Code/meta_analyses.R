
####################################################################
######################  RUN ALL META ANALYSES ######################
####################################################################




##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(meta)
library(ggpubr)


##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")


##### source functions
source("Functions/lab_correlations_functions.R")
source("Functions/meta_analysis_functions.R")

##### load data
droseu <- readRDS("Data/droseu_master_list_2022-03-25.rds")
estimates <- readRDS("LinearModelsPop/all_model_estimates.rds")
pops <- readRDS("InfoTables/DrosEU_Populations.rds")


##### create output directory
meta_dir <- "MetaAnalyses"
dir.create(meta_dir, showWarnings = F) 


##### define functions

makeEffects <- function(x) {
  x <- mutate(x, Population = factor(Population, levels = c("YE","RE","GI","MU","MA","UM","KA","VA","AK")), Lab = as.factor(Lab), V = SE^2, Study = paste(Population, Lab, sep = "_"))
  x <- relocate(x, Trait, Population, Sex, Lab, Study) %>% arrange(Population) %>% dplyr::rename(Y = Estimate) }





############# VIABILITY #############

# create output directory
out_dir <- "Viability"
dir.create(file.path(meta_dir, out_dir), showWarnings = F) 

# get effects
Via_effects <- makeEffects(estimates$via)

# run meta
Via_meta <- metagen(data = filter(Via_effects, Sex == "F"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)
# sub groups
Via_meta <- update.meta(Via_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(Via_meta, file = file.path(meta_dir, out_dir, "Via_meta.rds"))



############# DEVELOPMENT TIME #############

# create output directory
out_dir <- "DevelopmentTime"
dir.create(file.path(meta_dir, out_dir), showWarnings = F) 

# get effects
DT_A_effects <- makeEffects(filter(estimates$dta, Trait == "DT_A"))

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
DW_effects <- makeEffects(estimates$dw)

# run meta
DT_F_meta <- metagen(data = filter(DW_effects, Sex == "F"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)
# sub groups
DT_F_meta <- update.meta(DT_F_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(DT_F_meta, file = file.path(meta_dir, out_dir, "DT_F_meta.rds"))


# run meta
DT_M_meta <- metagen(data = filter(DW_effects, Sex == "M"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)
# sub groups
DT_M_meta <- update.meta(DT_M_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(DT_M_meta, file = file.path(meta_dir, out_dir, "DT_M_meta.rds"))






############# DRY WEIGHT #############

# create output directory
out_dir <- "ThoraxLength"
dir.create(file.path(meta_dir, out_dir), showWarnings = F) 

# get effects
TL_effects <- makeEffects(estimates$tl)

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
WA_effects <- makeEffects(estimates$wa)

# run meta for both L and R for females
WA_F_meta <- metagen(data = filter(WA_effects, Sex == "F"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)
# sub groups
WA_F_meta <- update.meta(WA_F_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(WA_F_meta, file = file.path(meta_dir, out_dir, "WA_F_meta.rds"))

# run meta for both L and R for males
WA_M_meta <- metagen(data = filter(WA_effects, Sex == "M"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)
# sub groups
WA_M_meta <- update.meta(WA_M_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(WA_M_meta, file = file.path(meta_dir, out_dir, "WA_M_meta.rds"))




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
Fec_effects <- makeEffects(estimates$fec)

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
LS_effects <- makeEffects(estimates$ls)

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
CSM_effects <- makeEffects(estimates$csm)

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










############# HEAT-SHOCK MORTALITY #############

# create output directory
out_dir <- "HeatShock"
dir.create(file.path(meta_dir, out_dir), showWarnings = F) 

# get effects
HSM_effects <- makeEffects(estimates$hsm)

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
Dia_effects <- makeEffects(estimates$dia)

# run meta for females
Dia_meta <- metagen(data = filter(Dia_effects, Sex == "F"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)

# sub groups
Dia_meta <- update.meta(Dia_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(Dia_meta, file = file.path(meta_dir, out_dir, "Dia_meta.rds"))











############# STARVATION RESISTANCE #############

# create output directory
out_dir <- "Starvation"
dir.create(file.path(meta_dir, out_dir), showWarnings = F) 

# get effects
SR_effects <- makeEffects(estimates$sr)

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
Pgm_effects <- makeEffects(estimates$pgm)

# run meta for females, w/o Schmidt's data
Pgm_meta <- metagen(data = filter(Pgm_effects, Sex == "F", Trait == "Pgm_Total",  Lab != "Schmidt"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)

# sub groups
Pgm_meta <- update.meta(Pgm_meta, subgroup = Population, tau.common = FALSE)

# save output
saveRDS(Pgm_meta, file = file.path(meta_dir, out_dir, "SR_F_meta.rds"))



# run meta for females, all three labs
#Pgm_meta <- metagen(data = filter(Pgm_effects, Sex == "F", Trait == "Pgm_Total"), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE)

# sub groups
#Pgm_meta <- update.meta(Pgm_meta, subgroup = Population, tau.common = FALSE)

# save output
#saveRDS(Pgm_meta, file = file.path(meta_dir, out_dir, "SR_F_meta.rds"))




############# ouput all metas ############# 

all_metas <- ls()[grep("_meta", ls())]
saveRDS(all_metas, file = file.path(meta_dir, "all_metas.rds"))






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

metas <- rdsBatchReaderToList(path = "MetaAnalyses", recursive = T, full.names = T, pattern = "meta.rds")

metas_pvalues <- bind_rows(Meta = names(metas), Pvalue = lapply(metas, function(x) x$pval.Q.b.random) %>% unlist())
metas_pvalues <- mutate(metas_pvalues, Pvalue_bonf = Pvalue * 12) %>%
  mutate(Pvalue_bonf = ifelse(Pvalue_bonf > 1, 1, Pvalue_bonf))



######### plots


for (i in 1:length(metas)){
  f <- metas[i]
  p_out_pdf <- sub(".rds", "_summary_effect.pdf", f)
  p_out_png <- sub(".rds", "_summary_effect.png", f)
  m <- readRDS(f)
  m <- data.frame(Population = m$bylevs, Mstar = m$TE.random.w, SEMstar = m$seTE.random.w, Q = m$Q.b.random, p = m$pval.Q.b.random, LLMstar = m$lower.random.w, ULMstar = m$upper.random.w ) %>% mutate(Q_plot = paste0("italic(Q) == ", round(Q, 2)), p_plot = ifelse(p < 0.001, "italic(p) < 0.001", paste0("italic(p) == ", round(p, 4))), Population = factor(Population, levels = pops$by_lat$Population))

  
  p_meta_SE <- ggplot(data = m, aes(x = Mstar, y = 1:length(Population), color = Population)) +
    theme_bw() +
    geom_point(size = 4, shape = 15) +
    geom_errorbarh(aes(xmax = Mstar + SEMstar, xmin = Mstar - SEMstar), height = 0) +
    scale_color_manual(values = pops$by_lat$Color) +
    scale_y_continuous(name = "Population", breaks = 1:length(m$Population), labels = m$Population) +
    labs(x = "Summary effect", title = "Pop. summary effect with SE") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    annotate("text", x = -Inf, y = Inf, label = unique(m$Q_plot), hjust=-0.2, vjust=3.7, parse = T) +
    #annotate("text", x = -Inf, y = Inf, label = unique(m$p_plot), hjust=-0.2, vjust=4.2, parse = T) +
    theme(legend.position = "none")
  
  p_meta_CI <- ggplot(data = m, aes(x = Mstar, y = 1:length(Population), color = Population)) +
    theme_bw() +
    geom_point(size = 4, shape = 15) +
    geom_errorbarh(aes(xmax = ULMstar, xmin = LLMstar), height = 0) +
    scale_color_manual(values = pops$by_lat$Color) +
    scale_y_continuous(name = "Population", breaks = 1:length(m$Population), labels = m$Population) +
    labs(x = "Summary effect", title = "Pop. summary effect with 95% CI") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    annotate("text", x = -Inf, y = Inf, label = unique(m$Q_plot), hjust=-0.2, vjust=3.7, parse = T) +
    #annotate("text", x = -Inf, y = Inf, label = unique(m$p_plot), hjust=-0.2, vjust=4.2, parse = T) +
    theme(legend.position = "none")
  
  p_meta <- ggarrange(p_meta_SE, p_meta_CI)
  
  pdf(p_out_pdf, width = 8, height = 5)
  print(p_meta)
  dev.off()
  
  png(filename = p_out_png, height = 1500, width = 2400, res = 300)
  print(p_meta)
  dev.off()
}
  







#str(m , list.len = length(m))


