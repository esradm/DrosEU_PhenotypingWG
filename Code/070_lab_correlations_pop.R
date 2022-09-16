
###############################################################
######################  LAB CORRELATIONS ######################
###############################################################


##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(rstatix)

##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")


##### source functions
source("Functions/lab_correlations_functions.R")



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
        ts_out_pdf <- file.path(d_out, paste0(unique(ts[[j]]$Trait), "_", unique(ts[[j]]$Sex), "_", unique(ts[[j]]$Model), "_lab_correlations.pdf"))
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
  


######


estimates <- list.files(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "pop_model_estimates.rds")

estimates <- estimates[-grep("LA_lmers", estimates)] 
estimates <- estimates[-grep("Dia_lmers", estimates)] 




cor_pearson <- list()
for (i in 1:length(estimates)) {
  # get data in and reformat
  f <- estimates[i]
  e <- readRDS(f) %>%
    dplyr::select(Trait, Sex, Population, Lab, Estimate) %>%
    pivot_wider(names_from = c("Lab"), values_from = Estimate) %>%
    group_split(Trait, Sex)
  # get correlations
  co <- list()
  for (j in 1:length(e)) {
    cols_ok <- !apply(e[[j]], 2, function(x) all(is.na(x)))
    e_sub <- dplyr::select(e[[j]], names(cols_ok)[cols_ok])
    if (ncol(e_sub) > 4) {
      if (ncol(e_sub) == 5) {
        co[[j]] <- cor_test(e_sub, -c(Sex, Population, Trait), method = "pearson")
      } else {
        co[[j]] <- cor_mat(e_sub, -c(Sex, Population, Trait), method = "pearson") %>% cor_gather()
      }
      co[[j]] <- filter(co[[j]], var1 != var2) %>% 
        rename(Lab1 = var1, Lab2 = var2, R = cor, P = p) %>%
        mutate(Sex = unique(e_sub$Sex),
               Trait = unique(e_sub$Trait)) %>%
        dplyr::select(Trait, Sex, Lab1, Lab2, R, P) 
      
      co[[j]]$Labs <- "NA"
      for (r in 1:nrow(co[[j]])) { 
        co[[j]]$Labs[r] <- co[[j]][r, c("Lab1", "Lab2")] %>% 
          unlist %>% sort %>% paste(collapse = "_") 
      }
      
      co[[j]] <- group_by(co[[j]], Trait, Sex, Labs) %>%
        distinct(Labs, .keep_all = TRUE)
      
    }
  }
  cor_pearson[[i]] <- bind_rows(co)
}
cor_pearson <- bind_rows(cor_pearson)

cor_pearson <- cor_pearson %>% 
  mutate(Method = "Pearson",
         Group = paste(Trait, Sex),
         Cutoff = as.factor(ifelse(P < 0.05, 0, 1)))







pop_cor <- ggplot(cor_pearson, aes(x = R, y = -log10(P))) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, size = 0.5, col = "red") +
  geom_vline(xintercept = 0, size = 0.5, col = "grey80") +
  geom_point(aes(color = Cutoff)) +
  #scale_color_manual(values = c("black", "red")) +
  facet_wrap(Group ~ ., ncol = 7) +
  theme_bw(14) +
  labs(x = "Pearson's correlation coefficient", title = "Lab correlations - Population level - Pearson's coefficients") +
  theme(legend.position = "none")


##### save facet plot
ggsave(pop_cor, filename = "LabCorrelations/lab_correlation_pop_pearson.pdf", width = 11.69, height = 8.27)











cor_spearman <- list()
for (i in 1:length(estimates)) {
  # get data in and reformat
  f <- estimates[i]
  e <- readRDS(f) %>%
    dplyr::select(Trait, Sex, Population, Lab, Estimate) %>%
    pivot_wider(names_from = c("Lab"), values_from = Estimate) %>%
    group_split(Trait, Sex)
  # get correlations
  co <- list()
  for (j in 1:length(e)) {
    cols_ok <- !apply(e[[j]], 2, function(x) all(is.na(x)))
    e_sub <- dplyr::select(e[[j]], names(cols_ok)[cols_ok])
    if (ncol(e_sub) > 4) {
      if (ncol(e_sub) == 5) {
        co[[j]] <- cor_test(e_sub, -c(Sex, Population, Trait), method = "spearman")
      } else {
        co[[j]] <- cor_mat(e_sub, -c(Sex, Population, Trait), method = "spearman") %>% cor_gather()
      }
      co[[j]] <- filter(co[[j]], var1 != var2) %>% 
        rename(Lab1 = var1, Lab2 = var2, R = cor, P = p) %>%
        mutate(Sex = unique(e_sub$Sex),
               Trait = unique(e_sub$Trait)) %>%
        dplyr::select(Trait, Sex, Lab1, Lab2, R, P) 
      
      co[[j]]$Labs <- "NA"
      for (r in 1:nrow(co[[j]])) { 
        co[[j]]$Labs[r] <- co[[j]][r, c("Lab1", "Lab2")] %>% 
          unlist %>% sort %>% paste(collapse = "_") 
      }
      
      co[[j]] <- group_by(co[[j]], Trait, Sex, Labs) %>%
        distinct(Labs, .keep_all = TRUE)
      
    }
  }
  cor_spearman[[i]] <- bind_rows(co)
}
cor_spearman <- bind_rows(cor_spearman)

cor_spearman <- cor_spearman %>% 
  mutate(Method = "Spearman",
         Group = paste(Trait, Sex),
         Cutoff = as.factor(ifelse(P < 0.05, 0, 1)))







pop_cor <- ggplot(cor_spearman, aes(x = R, y = -log10(P))) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, size = 0.5, col = "red") +
  geom_vline(xintercept = 0, size = 0.5, col = "grey80") +
  geom_point(aes(color = Cutoff)) +
  #scale_color_manual(values = c("black", "red")) +
  facet_wrap(Group ~ ., ncol = 7) +
  theme_bw(14) +
  labs(x = "Spearman's rank correlation coefficient", title = "Lab correlations - Population level - Spearman's rank coefficients") +
  theme(legend.position = "none")
        

##### save facet plot
ggsave(pop_cor, filename = "LabCorrelations/lab_correlation_pop_spearman.pdf", width = 11.69, height = 8.27)




cor_spearman <- list()
for (i in 1:length(estimates)) {
  # get data in and reformat
  f <- estimates[i]
  e <- readRDS(f) %>%
    dplyr::select(Trait, Sex, Population, Lab, Estimate) %>%
    pivot_wider(names_from = c("Lab"), values_from = Estimate) %>%
    group_split(Trait, Sex)
  # get correlations
  co <- list()
  for (j in 1:length(e)) {
    cols_ok <- !apply(e[[j]], 2, function(x) all(is.na(x)))
    e_sub <- dplyr::select(e[[j]], names(cols_ok)[cols_ok])
    if (ncol(e_sub) > 4) {
      if (ncol(e_sub) == 5) {
        co[[j]] <- cor_test(e_sub, -c(Sex, Population, Trait), method = "spearman")
      } else {
        co[[j]] <- cor_mat(e_sub, -c(Sex, Population, Trait), method = "spearman") %>% cor_gather()
      }
      co[[j]] <- filter(co[[j]], var1 != var2) %>% 
        rename(Lab1 = var1, Lab2 = var2, R = cor, P = p) %>%
        mutate(Sex = unique(e_sub$Sex),
               Trait = unique(e_sub$Trait)) %>%
        dplyr::select(Trait, Sex, Lab1, Lab2, R, P) 
      
      co[[j]]$Labs <- "NA"
      for (r in 1:nrow(co[[j]])) { 
        co[[j]]$Labs[r] <- co[[j]][r, c("Lab1", "Lab2")] %>% 
          unlist %>% sort %>% paste(collapse = "_") 
      }
      
      co[[j]] <- group_by(co[[j]], Trait, Sex, Labs) %>%
        distinct(Labs, .keep_all = TRUE)
      
    }
  }
  cor_spearman[[i]] <- bind_rows(co)
}
cor_ <- bind_rows(cor_spearman)

cor_spearman <- cor_spearman %>% 
  mutate(Method = "Spearman",
         Group = paste(Trait, Sex),
         Cutoff = as.factor(ifelse(P < 0.05, 0, 1)))









