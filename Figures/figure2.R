
##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(cowplot)
library(ggpubr)

##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")




traits <- read.csv("InfoTables/trait_names.csv")

traits <- dplyr::select(traits, Trait, Title)
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
  dplyr::select(Trait, Sex, Group, Marg_r2, Sig, Title)

pr2$Title_sex <- paste(pr2$Title, pr2$Sex, sep = " - ")
pr2$Title_sex <- gsub(" - NA", "", pr2$Title_sex)



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

#ggsave(marg_r2, filename = file.path(lmer_dir, "marginal_r2.pdf"), width = 10, height = 12)
#ggsave(marg_r2, filename = file.path(lmer_dir, "marginal_r2.png"), width = 10, height = 12)




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





ps <- ggdraw() +
  draw_plot(marg_r2, x = 0, y = 0, width = 0.63, height = 1) +
  draw_plot(line_pearson_alt +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  ), x = 0.68, y = 0, width = 0.35, height = 1) +
  draw_plot_label(c("A", "B"),
    x = c(0, 0.65),
    y = c(1, 1), size = 12
  )

ggsave(ps,
  filename = "Figures/figure2_v1.png",
  width = 6.3, height = 4
)


ggsave(ps,
  filename = "graphics/chrom_viab_figure_v4.0.pdf",
  width = 7, height = 9, dpi = 300, units = "in"
)


ggarrange(marg_r2, line_pearson_alt,  common.legend = TRUE)
plot_grid(marg_r2, line_pearson_alt, labels = c('A', 'B'))

cowplot::plot_grid(marg_r2, line_pearson_alt +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  ),
nrow = 1,
labels = "auto",
rel_widths = c(2, 1))







ggsave(line_pearson_alt,
  filename = file.path(cor_dir, "lab_correlation_line_pearson_v2.png"),
  width = 10, height = 12
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
