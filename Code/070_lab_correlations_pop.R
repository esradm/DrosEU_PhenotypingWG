
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
dir.create(cor_dir, showWarnings = FALSE)




########### PER TRAIT PAIRWISE PLOTS ##########


##### POP LEVEL

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
        png(ts_out_png, 800, 800, res = 120)
        scatterPlotMatrix(filter(ts[[j]], Trait == unique(ts[[j]]$Trait), Sex == unique(ts[[j]]$Sex)))
        dev.off()
      }
    }
  }
}



##### LINE LEVEL

line_estimates <- list.files(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "lmers_line_random_coefs.rds")

line_estimates <- line_estimates[-grep("Dia_lmers", line_estimates)] 



for (i in 1:length(line_estimates)){
  f <- line_estimates[i]
  e <- readRDS(f) %>% mutate(Predictor = "NA") %>% rename(Estimate = Coef)
  if (length(unique(e$Lab)) > 1) {
    d_out <- sub("LinearModelsPop", cor_dir, f)
    d_out <- str_match(d_out, '(.*[^/]+)(?:/[^/]+){1}$')[,2]
    dir.create(d_out, showWarnings = F) 
    ts <- group_split(e, Trait, Sex)
    for (j in 1:length(ts)){
      if (length(unique(ts[[j]]$Lab)) > 1) {
        ts_out_pdf <- file.path(d_out, paste0(unique(ts[[j]]$Trait), "_", unique(ts[[j]]$Sex), "_", unique(ts[[j]]$Model), "_lab_line_pearson_correlations.pdf"))
        pdf(ts_out_pdf)
        scatterPlotMatrix(filter(ts[[j]], Trait == unique(ts[[j]]$Trait), Sex == unique(ts[[j]]$Sex)))
        dev.off()
        ts_out_png <- sub("pdf", "png", ts_out_pdf)
        png(ts_out_png, 840, 840, res = 120)
        scatterPlotMatrix(filter(ts[[j]], Trait == unique(ts[[j]]$Trait), Sex == unique(ts[[j]]$Sex)))
        dev.off()
      }
    }
  }
}





########### POP LEVEL PEARSON FOR ALL TRAITS ##########

estimates <- list.files(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "pop_model_estimates.rds")

estimates <- estimates[-grep("LA_lmers", estimates)] 
estimates <- estimates[-grep("Dia_lmers", estimates)] 



##### get coefficients and p values
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

##### combine to data.frame
cor_pearson <- bind_rows(cor_pearson) %>% 
  mutate(Method = "Pearson",
         Group = paste(Trait, Sex),
         Cutoff = as.factor(ifelse(P < 0.05, 0, 1)))

##### get stats for plotting
stats_text <- cor_pearson %>% 
  group_by(Trait, Sex) %>%
  mutate(nl = paste0("n = ", length(unique(Lab1)) + 1), 
         nc = paste0("nc = ", n()), 
         ncs = paste0("ncs = ", sum(Cutoff == 0))) %>%
  ungroup() %>%
  mutate(ymax = max(-log10(P)),
         xmin = min(R)) %>%
  select(Group, nl, nc, ncs, xmin, ymax) %>%
  distinct()

##### plot
pop_cor <- ggplot(cor_pearson, aes(x = R, y = -log10(P))) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, size = 0.5, col = "red") +
  geom_vline(xintercept = 0, size = 0.5, col = "grey80") +
  geom_point(aes(fill = Cutoff), size = 2, pch = 21, color = "black") +
  facet_wrap(Group ~ ., ncol = 7) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1)) +
  theme_bw(14) +
  labs(x = "Pearson's correlation coefficient", title = "Pairwise lab correlations - Population level - Pearson's coefficients", subtitle = "n, number of labs; nc, number of pairwise lab combinations; ncs, number of significant correlations") +
  theme(legend.position = "none",
        panel.grid = element_blank()) +
  geom_text(data = stats_text, aes(label = paste(nl, nc, ncs, sep = "\n"), x = xmin, y = ymax), hjust = 0, vjust = 1, size = 3.5)

##### save facet plot
ggsave(pop_cor, filename = "LabCorrelations/lab_correlation_pop_pearson.png", width = 11.69, height = 8.27, dpi = 120)









########### POP LEVEL SPEARMAN FOR ALL TRAITS ##########

estimates <- list.files(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "pop_model_estimates.rds")

estimates <- estimates[-grep("LA_lmers", estimates)] 
estimates <- estimates[-grep("Dia_lmers", estimates)] 


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



stats_text <- cor_spearman %>% 
  group_by(Trait, Sex) %>%
  mutate(nl = paste0("n = ", length(unique(Lab1)) + 1), 
         nc = paste0("nc = ", n()), 
         ncs = paste0("ncs = ", sum(Cutoff == 0))) %>%
  ungroup() %>%
  mutate(ymax = max(-log10(P)),
         xmin = min(R)) %>%
  select(Group, nl, nc, ncs, xmin, ymax) %>%
  distinct()


pop_cor <- ggplot(cor_spearman, aes(x = R, y = -log10(P))) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, size = 0.5, col = "red") +
  geom_vline(xintercept = 0, size = 0.5, col = "grey80") +
  geom_point(aes(fill = Cutoff), size = 2, pch = 21, color = "black") +
  facet_wrap(Group ~ ., ncol = 7) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1)) +
  theme_bw(14) +
  labs(x = "Spearman's correlation coefficient", title = "Pairwise lab correlations - Population level - Spearman's coefficients", subtitle = "n, number of labs; nc, number of pairwise lab combinations; ncs, number of significant correlations") +
  theme(legend.position = "none",
        panel.grid = element_blank()) +
  geom_text(data = stats_text, aes(label = paste(nl, nc, ncs, sep = "\n"), x = xmin, y = ymax), hjust = 0, vjust = 1, size = 3.5)


##### save facet plot
ggsave(pop_cor, filename = "LabCorrelations/lab_correlation_pop_spearman.png", width = 11.69, height = 8.27, dpi = 120)





########### LINE LEVEL PEARSON FOR ALL TRAITS ##########



line_estimates <- list.files(path = "LinearModelsPop", recursive = T, full.names = T, pattern = "lmers_line_random_coefs.rds")

line_estimates <- line_estimates[-grep("LA_lmers", line_estimates)] 
line_estimates <- line_estimates[-grep("Dia_lmers", line_estimates)] 




cor_pearson <- list()
for (i in 1:length(line_estimates)) {
  # get data in and reformat
  f <- line_estimates[i]
  e <- readRDS(f) %>%
    rename(Estimate = Coef) %>%
    dplyr::select(Trait, Sex, Line, Lab, Estimate) %>%
    pivot_wider(names_from = c("Lab"), values_from = Estimate) %>%
    group_split(Trait, Sex)
  # get correlations
  co <- list()
  for (j in 1:length(e)) {
    cols_ok <- !apply(e[[j]], 2, function(x) all(is.na(x)))
    e_sub <- dplyr::select(e[[j]], names(cols_ok)[cols_ok])
    if (ncol(e_sub) > 4) {
      if (ncol(e_sub) == 5) {
        co[[j]] <- cor_test(e_sub, -c(Sex, Line, Trait), method = "pearson")
      } else {
        co[[j]] <- cor_mat(e_sub, -c(Sex, Line, Trait), method = "pearson") %>% cor_gather()
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



stats_text <- cor_pearson %>% 
  group_by(Trait, Sex) %>%
  mutate(nl = paste0("n = ", length(unique(Lab1)) + 1), 
         nc = paste0("nc = ", n()), 
         ncs = paste0("ncs = ", sum(Cutoff == 0))) %>%
  ungroup() %>%
  mutate(ymax = max(-log10(P)),
         xmin = min(R)) %>%
  select(Group, nl, nc, ncs, xmin, ymax) %>%
  distinct()



line_cor <- ggplot(cor_pearson, aes(x = R, y = -log10(P))) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, size = 0.5, col = "red") +
  geom_vline(xintercept = 0, size = 0.5, col = "grey80") +
  geom_point(aes(fill = Cutoff), size = 2, pch = 21, color = "black") +
  facet_wrap(Group ~ ., ncol = 7) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5)) +
  theme_bw(14) +
  labs(x = "Pearson's correlation coefficient", title = "Pairwise lab correlations - Line level - Pearson's coefficients", subtitle = "n, number of labs; nc, number of pairwise lab combinations; ncs, number of significant correlations") +
  theme(legend.position = "none",
        panel.grid = element_blank()) +
  geom_text(data = stats_text, aes(label = paste(nl, nc, ncs, sep = "\n"), x = xmin, y = ymax), hjust = 0, vjust = 1, size = 3.5)



ggsave(line_cor, filename = "LabCorrelations/lab_correlation_line_pearson.png", width = 11.69, height = 8.27, dpi = 120)
















