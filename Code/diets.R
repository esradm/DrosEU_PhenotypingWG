


##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(ggbeeswarm)
library(ggdist)
library(ggpubr)
#library(lme4)
#library(lsmeans)
#library(afex)
#library(multcomp)



##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")


##### source functions
#source("Functions/lab_correlations_functions.R")

##### load data
droseu <- readRDS("Data/droseu_master_list_2022-05-02.rds")
diets <- read.csv("InfoTables/DrosEU_Diets_Feb22.csv")



p1 <- diets %>% select(PI.Lab.head,  P.C, Diet) %>% distinct() %>% 
  ggplot(aes(x = Diet, y = P.C, color = Diet)) +
  geom_boxplot(outlier.colour = NA) +
  geom_quasirandom(size = 4, alpha = 0.5, width = .3) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none") +
  ggtitle("All ratios") +
  theme(plot.title = element_text(size=22))

p2 <- diets %>% filter(P.C < 1) %>% select(PI.Lab.head,  P.C, Diet) %>% distinct() %>%
  ggplot(aes(x = Diet, y = P.C, color = Diet)) +
  geom_boxplot(outlier.colour = NA) +
  geom_quasirandom(size = 4, alpha = 0.5, width = .3) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none") +
  ggtitle("Without extreme ratios") +
  theme(plot.title = element_text(size=22))


ps <- ggarrange(p1, p2)

ggsave(ps, filename = "InfoTables/DrosEU_Diets_PC_ratios.pdf", width = 10, height = 8)





labs <- lapply(droseu, function(x) as.data.frame(unique(x$Supervisor.PI)))
traits <- rep(names(droseu), lapply(labs, nrow))
d <- data.frame(trait = traits, lab = as.character(unlist(bind_rows(labs))))
rownames(d) <- NULL
d$lab[d$lab == "StamenkovicRadak"] <- "Stamenkovic-Radak"


d <- left_join(d, diets %>% rename(lab = PI.Lab.head) %>% select(lab, Diet, P.C) %>% distinct)


p3 <- d %>% filter(!is.na(P.C)) %>%
  ggplot(aes(x = P.C, y = trait, color = Diet)) +
  geom_point(size = 3, alpha = 0.5) +
  theme_classic(base_size = 24) +
  labs(y = "Traits", title = "All ratios") +
  theme(panel.grid.major.y = element_line(size = 0.5)) +
  theme(plot.title = element_text(size=22))

p4 <- d %>% filter(P.C < 1) %>%
  ggplot(aes(x = P.C, y = trait, color = Diet)) +
  geom_point(size = 3, alpha = 0.5) +
  theme_classic(base_size = 24) +
  labs(y = "Traits", title = "Without extreme ratios") +
  theme(panel.grid.major.y = element_line(size = 0.5)) +
  theme(plot.title = element_text(size=22))


ggsave(ggarrange(p3, p4, common.legend = TRUE), filename = "InfoTables/DrosEU_Diets_PC_ratios_traits.pdf", width = 10, height = 8)


