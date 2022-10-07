
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

estimates <- list.files(
  path = "LinearModelsPop", recursive = TRUE, full.names = TRUE,
  pattern = "pop_model_estimates.rds"
)

for (i in seq_len(length(estimates))) {
  f <- estimates[i]
  e <- readRDS(f)
  if (length(unique(e$Lab)) > 1) {
    d_out <- sub("LinearModelsPop", cor_dir, f)
    d_out <- str_match(d_out, "(.*[^/]+)(?:/[^/]+){1}$")[, 2]
    dir.create(d_out, showWarnings = FALSE)
    ts <- group_split(e, Trait, Sex)
    for (j in seq_len(length(ts))) {
      if (length(unique(ts[[j]]$Lab)) > 1) {
        ts_out_pdf <- file.path(d_out, paste0(
          unique(ts[[j]]$Trait), "_", unique(ts[[j]]$Sex), "_",
          unique(ts[[j]]$Model), "_lab_correlations.pdf"
        ))
        pdf(ts_out_pdf)
        scatterPlotMatrix(filter(
          ts[[j]], Trait == unique(ts[[j]]$Trait),
          Sex == unique(ts[[j]]$Sex)
        ))
        dev.off()
        ts_out_png <- sub("pdf", "png", ts_out_pdf)
        png(ts_out_png, 800, 800, res = 120)
        scatterPlotMatrix(filter(
          ts[[j]], Trait == unique(ts[[j]]$Trait),
          Sex == unique(ts[[j]]$Sex)
        ))
        dev.off()
      }
    }
  }
}



##### LINE LEVEL

line_estimates <- list.files(
  path = "LinearModelsPop", recursive = TRUE, full.names = TRUE,
  pattern = "lmers_line_random_coefs.rds"
)

line_estimates <- line_estimates[-grep("Dia_lmers", line_estimates)]

for (i in seq_len(length(line_estimates))) {
  f <- line_estimates[i]
  e <- readRDS(f) %>%
    mutate(Predictor = "NA") %>%
    rename(Estimate = Coef)
  if (length(unique(e$Lab)) > 1) {
    d_out <- sub("LinearModelsPop", cor_dir, f)
    d_out <- str_match(d_out, "(.*[^/]+)(?:/[^/]+){1}$")[, 2]
    dir.create(d_out, showWarnings = FALSE)
    ts <- group_split(e, Trait, Sex)
    for (j in seq_len(length(ts))) {
      if (length(unique(ts[[j]]$Lab)) > 1) {
        ts_out_pdf <- file.path(d_out, paste0(
          unique(ts[[j]]$Trait), "_", unique(ts[[j]]$Sex), "_",
          unique(ts[[j]]$Model), "_lab_line_pearson_correlations.pdf"
        ))
        pdf(ts_out_pdf)
        scatterPlotMatrix(filter(
          ts[[j]], Trait == unique(ts[[j]]$Trait),
          Sex == unique(ts[[j]]$Sex)
        ))
        dev.off()
        ts_out_png <- sub("pdf", "png", ts_out_pdf)
        png(ts_out_png, 840, 840, res = 120)
        scatterPlotMatrix(filter(
          ts[[j]], Trait == unique(ts[[j]]$Trait),
          Sex == unique(ts[[j]]$Sex)
        ))
        dev.off()
      }
    }
  }
}





########### POP LEVEL PEARSON FOR ALL TRAITS ##########

estimates <- list.files(
  path = "LinearModelsPop", recursive = TRUE, full.names = TRUE,
  pattern = "pop_model_estimates.rds"
)

estimates <- estimates[-grep("LA_lmers", estimates)]
estimates <- estimates[-grep("Dia_lmers", estimates)]

##### get coefficients and p values

lab_correlation <- function(estimates_path = estimates,
  method = "pearson") {
    correlation_list <- list(length = length(estimates_path))
    for (i in seq_len(length(estimates_path))) {
      f <- estimates_path[i]
      e <- readRDS(f)
      if ("Coef" %in% colnames(e)) {
        e <- dplyr::rename(e, Estimate = Coef)
      }
      cols_to_rm <- c("Model", "Predictor", "SE")
      if ("Line" %in% colnames(e)) {
        cols_to_rm <- c("Model", "Predictor", "SE", "Population")
      }
      e <- e[, !colnames(e) %in% cols_to_rm] %>%
        tidyr::pivot_wider(names_from = Lab, values_from = Estimate) %>%
        dplyr::group_split(Trait, Sex)
      co <- list()
      for (j in seq_len(length(e))) {
        cols_ok <- !apply(e[[j]], 2, function(x) all(is.na(x)))
        e_sub <- dplyr::select(e[[j]], names(cols_ok)[cols_ok])
        if (ncol(e_sub) > 4) {
          e_sub2 <- e_sub[, !colnames(e_sub) %in% c("Sex", "Population", "Line", "Trait")]
          if (ncol(e_sub2) == 2) {
            co[[j]] <- cor_test(e_sub2, method = method[1])
          } else {
            co[[j]] <- cor_mat(e_sub2, method = method[1]) %>%
              cor_gather()
          }
        co[[j]] <- filter(co[[j]], var1 != var2) %>%
          dplyr::rename(Lab1 = var1, Lab2 = var2, R = cor, P = p) %>%
          mutate(
          Sex = unique(e_sub$Sex),
          Trait = unique(e_sub$Trait)
          ) %>%
          dplyr::select(Trait, Sex, Lab1, Lab2, R, P)

        co[[j]]$Labs <- "NA"
        for (r in seq_len(nrow(co[[j]]))) {
          co[[j]]$Labs[r] <- co[[j]][r, c("Lab1", "Lab2")] %>%
            unlist() %>%
            sort() %>%
            paste(collapse = "_")
        }

        co[[j]] <- dplyr::group_by(co[[j]], Trait, Sex, Labs) %>%
          dplyr::distinct("Labs", .keep_all = TRUE)
      }
    }
    correlation_list[[i]] <- dplyr::bind_rows(co)
    print(co)
  }
  correlation_list
}

pop_pearson_list <- lab_correlation(estimates, "pearson")

##### combine to data.frame

pop_pearson <- bind_rows(pop_pearson_list) %>%
  mutate(
    Method = "Pearson",
    Group = paste(Trait, Sex),
    Cutoff = factor(ifelse(P < 0.05, 0, 1), levels = c(0, 1))
  )

##### get stats for plotting

pop_pearson_stats_text <- pop_pearson %>%
  group_by(Trait, Sex) %>%
  mutate(
    nl = paste0("n = ", length(unique(Lab1)) + 1),
    nc = paste0("nc = ", n()),
    ncs = paste0("ncs = ", sum(Cutoff == 0))
  ) %>%
  ungroup() %>%
  mutate(
    ymax = max(-log10(P)),
    xmin = min(R)
  ) %>%
  select(Group, nl, nc, ncs, xmin, ymax) %>%
  distinct()


##### plot

pop_pearson_facet <- ggplot(pop_pearson, aes(x = R, y = -log10(P))) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, size = 0.5, col = "red") +
  geom_vline(xintercept = 0, size = 0.5, col = "grey80") +
  geom_point(aes(fill = Cutoff), size = 2, pch = 21, color = "black") +
  facet_wrap(Group ~ ., ncol = 7) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1)) +
  theme_bw(14) +
  labs(
    x = "Pearson's correlation coefficient",
    title = "Pairwise lab correlations - Population level - Pearson's coefficients",
    subtitle = "n, number of labs; nc, number of pairwise lab combinations; ncs, number of significant correlations"
  ) +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  ) +
  geom_text(
    data = pop_pearson_stats_text,
    aes(label = paste(nl, nc, ncs, sep = "\n"), x = xmin, y = ymax),
    hjust = 0, vjust = 1, size = 3.5
  )

##### save facet plot

ggsave(
  pop_pearson_facet,
  filename = "LabCorrelations/lab_correlation_pop_pearson.png",
  width = 11.69, height = 8.27, dpi = 120
)

ggsave(
  pop_pearson_facet,
  filename = "LabCorrelations/lab_correlation_pop_pearson.pdf",
  width = 11.69, height = 8.27
)










########### POP LEVEL SPEARMAN FOR ALL TRAITS ##########

estimates <- list.files(
  path = "LinearModelsPop", recursive = TRUE, full.names = TRUE,
  pattern = "pop_model_estimates.rds"
)

estimates <- estimates[-grep("LA_lmers", estimates)]
estimates <- estimates[-grep("Dia_lmers", estimates)]

pop_spearman_list <- lab_correlation(estimates, "spearman")

pop_spearman <- bind_rows(pop_spearman_list) %>%
  mutate(
    Method = "Spearman",
    Group = paste(Trait, Sex),
    Cutoff = factor(ifelse(P < 0.05, 0, 1), levels = c(0, 1))
  )

pop_spearman_stats_text <- pop_spearman %>%
  group_by(Trait, Sex) %>%
  mutate(
    nl = paste0("n = ", length(unique(Lab1)) + 1),
    nc = paste0("nc = ", n()),
    ncs = paste0("ncs = ", sum(Cutoff == 0))
  ) %>%
  ungroup() %>%
  mutate(
    ymax = max(-log10(P)),
    xmin = min(R)
  ) %>%
  select(Group, nl, nc, ncs, xmin, ymax) %>%
  distinct()


pop_spearman_facet <- ggplot(pop_spearman, aes(x = R, y = -log10(P))) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, size = 0.5, col = "red") +
  geom_vline(xintercept = 0, size = 0.5, col = "grey80") +
  geom_point(aes(fill = Cutoff), size = 2, pch = 21, color = "black") +
  facet_wrap(Group ~ ., ncol = 7) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1)) +
  theme_bw(14) +
  labs(x = "Spearman's correlation coefficient",
    title = "Pairwise lab correlations - Population level - Spearman's coefficients",
    subtitle = "n, number of labs; nc, number of pairwise lab combinations; ncs, number of significant correlations") +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  ) +
  geom_text(
    data = pop_spearman_stats_text,
    aes(label = paste(nl, nc, ncs, sep = "\n"), x = xmin, y = ymax),
    hjust = 0, vjust = 1, size = 3.5
  )


##### save facet plot

ggsave(pop_spearman_facet,
  filename = "LabCorrelations/lab_correlation_pop_spearman.png",
  width = 11.69, height = 8.27, dpi = 120
)

ggsave(pop_spearman_facet,
  filename = "LabCorrelations/lab_correlation_pop_spearman.pdf",
  width = 11.69, height = 8.27
)





########### LINE LEVEL PEARSON FOR ALL TRAITS ##########


line_estimates <- list.files(
  path = "LinearModelsPop", recursive = TRUE, full.names = TRUE,
  pattern = "lmers_line_random_coefs.rds"
)

line_estimates <- line_estimates[-grep("LA_lmers", line_estimates)]
line_estimates <- line_estimates[-grep("Dia_lmers", line_estimates)]

line_pearson_list <- lab_correlation(line_estimates, "pearson")

line_pearson <- bind_rows(line_pearson_list) %>%
  mutate(
    Method = "Pearson",
    Group = paste(Trait, Sex),
    Cutoff = factor(ifelse(P < 0.05, 0, 1), levels = c(0, 1))
  )



line_pearson_stats_text <- line_pearson %>%
  group_by(Trait, Sex) %>%
  mutate(
    nl = paste0("n = ", length(unique(Lab1)) + 1),
    nc = paste0("nc = ", n()),
    ncs = paste0("ncs = ", sum(Cutoff == 0))
  ) %>%
  ungroup() %>%
  mutate(
    ymax = max(-log10(P)),
    xmin = min(R)
  ) %>%
  select(Group, nl, nc, ncs, xmin, ymax) %>%
  distinct()



line_pearson_facet <- ggplot(line_pearson, aes(x = R, y = -log10(P))) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, size = 0.5, col = "red") +
  geom_vline(xintercept = 0, size = 0.5, col = "grey80") +
  geom_point(aes(fill = Cutoff), size = 2, pch = 21, color = "black") +
  facet_wrap(Group ~ ., ncol = 7) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5)) +
  theme_bw(14) +
  labs(x = "Pearson's correlation coefficient",
    title = "Pairwise lab correlations - Line level - Pearson's coefficients",
    subtitle = "n, number of labs; nc, number of pairwise lab combinations; ncs, number of significant correlations") +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  ) +
  geom_text(
    data = line_pearson_stats_text,
    aes(label = paste(nl, nc, ncs, sep = "\n"), x = xmin, y = ymax),
    hjust = 0, vjust = 1, size = 3.5
  )

ggsave(line_pearson_facet,
  filename = "LabCorrelations/lab_correlation_line_pearson.png",
  width = 11.69, height = 8.27, dpi = 120
)

ggsave(line_pearson_facet,
  filename = "LabCorrelations/lab_correlation_line_pearson.pdf",
  width = 11.69, height = 8.27
)


########### LINE LEVEL SPEARMAN FOR ALL TRAITS ##########

line_spearman_list <- lab_correlation(line_estimates, "spearman")

line_spearman <- bind_rows(line_spearman_list) %>%
  mutate(
    Method = "Spearman",
    Group = paste(Trait, Sex),
    Cutoff = factor(ifelse(P < 0.05, 0, 1), levels = c(0, 1)),
    P_plot = ifelse(P == 0, 10^-25, P)
  )



line_spearman_stats_text <- line_spearman %>%
  group_by(Trait, Sex) %>%
  mutate(
    nl = paste0("n = ", length(unique(Lab1)) + 1),
    nc = paste0("nc = ", n()),
    ncs = paste0("ncs = ", sum(Cutoff == 0))
  ) %>%
  ungroup() %>%
  mutate(
    ymax = max(-log10(P_plot)),
    xmin = min(R)
  ) %>%
  select(Group, nl, nc, ncs, xmin, ymax) %>%
  distinct()



line_spearman_facet <- ggplot(line_spearman, aes(x = R, y = -log10(P_plot))) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, size = 0.5, col = "red") +
  geom_vline(xintercept = 0, size = 0.5, col = "grey80") +
  geom_point(aes(fill = Cutoff), size = 2, pch = 21, color = "black") +
  facet_wrap(Group ~ ., ncol = 7) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5)) +
  theme_bw(14) +
  labs(x = "Spearman's correlation coefficient",
    title = "Pairwise lab correlations - Line level - Spearman's coefficients",
    subtitle = "n, number of labs; nc, number of pairwise lab combinations; ncs, number of significant correlations",
    y = "-log10(P)") +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  ) +
  geom_text(
    data = line_spearman_stats_text,
    aes(label = paste(nl, nc, ncs, sep = "\n"), x = xmin, y = ymax),
    hjust = 0, vjust = 1, size = 3.5
  )



ggsave(line_spearman_facet,
  filename = "LabCorrelations/lab_correlation_line_spearman.png",
  width = 11.69, height = 8.27, dpi = 120
)

ggsave(line_spearman_facet,
  filename = "LabCorrelations/lab_correlation_line_spearman.pdf",
  width = 11.69, height = 8.27
)
