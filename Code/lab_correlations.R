
###############################################################
######################  LAB COORELATIONS ######################
###############################################################


##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)


##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")


##### source functions
source("Functions/lab_correlations_functions.R")


##### load data
droseu <- readRDS("Data/droseu_master_list_2022-03-25.rds")
estimates <- readRDS("LinearModelsPop/all_model_estimates.rds")


##### create output directory
cor_dir <- "LabCorrelations"
dir.create(cor_dir, showWarnings = F) 




##### plots

estimates <- list.files(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "pop_model_estimates.rds")

for (i in 1:length(estimates)){
  f <- estimates[i]
  e <- readRDS(f)
  if (length(unique(e$Lab)) > 1) {
    d_out <- sub("LinearModelsPop", cor_dir, f)
    d_out <- str_match(d_out, '(.*[^/]+)(?:/[^/]+){1}$')[,2]
    dir.create(d_out, showWarnings = F) 
    ts <- group_split(e, Trait, Sex)
    for (j in 1:length(ts)){
      if (length(unique(ts[[j]]$Lab)) > 1) {
        ts_out_pdf <- file.path(d_out, paste0(unique(ts[[j]]$Trait), "_", unique(ts[[j]]$Sex), "_lab_correlations.pdf"))
        pdf(ts_out_pdf)
        scatterPlotMatrix(filter(ts[[j]], Trait == unique(ts[[j]]$Trait), Sex == unique(ts[[j]]$Sex)))
        dev.off()
        ts_out_png <- sub("pdf", "png", ts_out_pdf)
        png(ts_out_png, 2100, 2100, res = 300)
        scatterPlotMatrix(filter(ts[[j]], Trait == unique(ts[[j]]$Trait), Sex == unique(ts[[j]]$Sex)))
        dev.off()
      }
    }
  }
}
  
   







