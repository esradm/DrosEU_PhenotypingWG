
##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(cowplot)
library(ggpubr)
library(ggrepel)

##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")




traits <- read.csv("InfoTables/trait_names.csv")

traits <- dplyr::select(traits, Trait, Title, Plot_order)
traits$Trait[traits$Trait == "LA_ND"] <- "LA_NDlog2"
traits$Trait[traits$Trait == "LSM"] <- "LS"

r2s <- readRDS("LinearModelsPop/all_models_pop_r2.rds")
pvals <- readRDS("LinearModelsPop/all_models_pop_pvalues.rds")

pr2 <- inner_join(r2s, pvals) %>%
  inner_join(traits) %>%
  mutate(
    Group = paste(Trait, Sex, sep = " "),
    Sig = as.factor(ifelse(P < 0.05, 1, 0))
  ) %>%
  arrange(desc(Plot_order)) %>%
  dplyr::select(Trait, Sex, Group, Marg_r2, Sig, Title, Plot_order)


pr2$Title_sex <- paste(pr2$Title, pr2$Sex, sep = " - ")
pr2$Title_sex <- gsub(" - NA", "", pr2$Title_sex)
pr2$Title_sex <- factor(pr2$Title_sex, levels = unique(pr2$Title_sex))



marg_r2 <- ggplot(data = pr2, aes(x = Marg_r2, y = Title_sex, fill = Sig)) +
  geom_point(size = 2, shape = 21, alpha = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c("grey50", "red")) +
  theme(
    panel.grid.major.y = element_line(size = 0.5),
    legend.position = "none",
    axis.text = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.title = element_text(size = 8)
  ) +
  #scale_y_discrete(labels = unique(pr2$Title_sex)) +
  labs(
    x = "Marginal R2"
  )




line_pearson_list <- readRDS("LabCorrelations/lab_correlation_line_pearson_list.rds")



line_pearson <- bind_rows(line_pearson_list) %>%
  mutate(
    Method = "Pearson",
    Group = paste(Trait, Sex),
    Cutoff = factor(ifelse(P < 0.05, 0, 1), levels = c(0, 1))
  ) %>%
  ungroup() %>%
  dplyr::select(Group, R, Cutoff)


no_pearson <- data.frame(
  Group = pr2$Group[!pr2$Group %in% line_pearson$Group],
  R = 2,
  Cutoff = as.factor(0)
)

line_pearson <- bind_rows(line_pearson, no_pearson) %>%
  inner_join(dplyr::select(pr2, Group, Title_sex)) %>%
  distinct()






line_pearson_alt <- ggplot(line_pearson, aes(x = R, y = Title_sex)) +
  geom_vline(xintercept = 0, size = 0.5, col = "grey80") +
  geom_point(aes(fill = Cutoff), size = 2, pch = 21, alpha = 0.5) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1)) +
  theme_classic() +
  scale_fill_manual(values = c("red", "grey50")) +
  coord_cartesian(x = c(-0.5, 0.75)) +
  labs(
    x = "Pearson's R",
    y = "Trait"
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    panel.grid.major.y = element_line(size = 0.5)
  )








##### meta pvalues
metas_pvals <- readRDS("MetaAnalyses/all_models_pop_meta_pvalues.rds")





metas_pvals <- metas_pvals %>%
  mutate(
    Group = paste(Trait, Sex),
    Cutoff = factor(ifelse(P_bonf < 0.05, 0, 1), levels = c(0, 1))
  ) %>%
  dplyr::select(Group, Q, Cutoff)




no_metas <- data.frame(
  Group = pr2$Group[!pr2$Group %in% metas_pvals$Group],
  Q = -100,
  Cutoff = as.factor(0)
)

metas_pvals <- bind_rows(metas_pvals, no_metas) %>%
  inner_join(dplyr::select(pr2, Group, Title_sex)) %>%
  distinct()




metas_alt <- ggplot(metas_pvals, aes(x = Q, y = Title_sex)) +
  geom_point(aes(fill = Cutoff), size = 2.5, pch = 21, alpha = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c("red", "grey50")) +
  coord_cartesian(x = c(-10, 110)) +
  scale_x_continuous(breaks = c(0, 50, 100)) +
  labs(x = "Q") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    panel.grid.major.y = element_line(size = 0.5)
  )



r2 <- marg_r2 +
  labs(x = expression(Marginal~italic(R^2))) +
  scale_fill_manual(values = c("white","red"))

r <- line_pearson_alt +
  labs(x = expression(Correlation~coefficient~(italic(r)))) +
  scale_fill_manual(values = c("red", "white")) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.y = element_blank()
  )

q <- metas_alt +
  labs(x = expression(Heterogeneity~(italic(Q)))) +
  scale_fill_manual(values = c("red", "white")) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.y = element_blank()
  )

aligned_plots <- align_plots(r2, r, q, align = "h")

ps <- ggdraw() +
  draw_plot(aligned_plots[[1]], x = 0, y = 0, width = 0.58, height = 1) +
  draw_plot(aligned_plots[[2]], x = 0.59, y = 0, width = 0.25, height = 1) +
  draw_plot(aligned_plots[[3]], x = 0.85, y = 0, width = 0.15, height = 1) +
  draw_plot_label(c("A", "B", "C"),
    x = c(0, 0.57, 0.83),
    y = c(1, 1, 1), size = 12
  )

ggsave(ps,
  filename = "Figures/figure2_v5.png",
  width = 6.3, height = 4.4
)









############## below is deprecated

############# Q AND P VALUES PLOT #############

bh_thresh <- sum(sort(metas_pvals$P) < 0.05) / nrow(metas_pvals) * 0.05
bonf_thresh <- 0.05 / nrow(metas_pvals)


pvalue_plot <- metas_pvals %>%
  ggplot(aes(x = Q, y = -log10(P))) +
  geom_point(size = 1.5, alpha = 0.3, pch = 19) +
  theme_classic() +
  geom_hline(yintercept = -log10(bonf_thresh), linetype = 2, size = 0.5, col = "red") +
  geom_hline(yintercept = -log10(bh_thresh), linetype = 2, size = 0.5) +
  labs(x = "Q value", y = "-log10(Pvalue)") +
  geom_label_repel(
    data = mutate(metas_pvals, Label = ifelse(P_bh >= 0.05, NA, paste(Trait, Sex))),
    aes(x = Q, y = -log10(P), label = Label),
    size = 2, box.padding = 0.2, label.padding = 0.2,
    point.padding = 0, segment.color = "grey50",
    max.overlaps = Inf, min.segment.length = 0, seed = 1, force = 100
  ) +
  annotate(
    "text", x = max(metas_pvals$Q), y = -log10(bonf_thresh) + 0.2, size = 2,
    label = "Bonferroni threshold", color = "red", hjust = 1, vjust = 0
  ) +
  annotate(
    "text", x = max(metas_pvals$Q), y = -log10(bh_thresh) + 0.2, size = 2,
    label = "BH threshold", hjust = 1, vjust = 0
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
  )





ps <- ggdraw() +
  draw_plot(marg_r2, x = 0, y = 0.4, width = 0.63, height = 0.6) +
  draw_plot(line_pearson_alt +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    ), x = 0.68, y = 0.4, width = 0.35, height = 0.6) +
  draw_plot(pvalue_plot, x = 0, y = 0, width = 0.4, height = 0.4) +
  draw_plot_label(c("A", "B", "C"),
    x = c(0, 0.65, 0),
    y = c(1, 1, 0.4), size = 12
  )

#ggsave(ps,
#  filename = "Figures/figure2_v1.png",
#  width = 6.3, height = 6.3
#)






