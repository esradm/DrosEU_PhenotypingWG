
##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(MetBrewer)


##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")

##### create directory
dir.create("Inversions")

##### load data
invs <- read.csv("./Data/DrosEU_lines_inversions_typing.csv")
pops <- readRDS("InfoTables/DrosEU_Populations.rds")
droseu <- readRDS("Data/droseu_master_list_2022-05-02.rds")
pop_comp <- read.csv("MetaAnalyses/all_models_pop_meta_compound_estimates.csv")
line_comp <- read.csv("MetaAnalyses/all_models_line_meta_compound_random_coefs.csv")

##### define colors for plotting

myColors <- met.brewer("Johnson", 9)
names(myColors) <- as.factor(c("AK", "GI", "KA", "MA", "MU", "RE", "UM", "VA", "YE"))
colScale <- scale_colour_manual(name = "Population", values = myColors)




##### inv average frequency

invs_pop <- invs %>%
  mutate(
    population = str_extract(line_id, "^.{2}"),
    In2Lt = ifelse(In2Lt == 0.5 | In2Lt == 1, 1, 0)
  ) %>%
  group_by(population) %>%
  summarise_at(vars(In3RP:In3Rmo), mean) %>%
  rename(Population = population)



##### pop level

pop_list <- filter(pop_comp, Trait != "Dia_lmer") %>%
  dplyr::select(Trait, Sex, Population, Estimate) %>%
  inner_join(invs_pop) %>%
  inner_join(pops$by_lat) %>%
  group_split(Trait, Sex)




pop_2lt_pearson <- list()
for (i in 1:length(pop_list)) {
  cortest <- cor.test(pop_list[[i]]$Estimate, pop_list[[i]]$In2Lt)
  pop_2lt_pearson[[i]] <- data.frame(
    Trait = unique(pop_list[[i]]$Trait), 
    Sex = unique(pop_list[[i]]$Sex),
    R = cortest$estimate,
    P = cortest$p.value,
    Method = "pearson")
}

pop_2lt_pearson <- bind_rows(pop_2lt_pearson) %>%
  mutate(
    Sex = factor(Sex, levels = c("F", "M", "B", "NA")),
    Label = paste(Trait, Sex, sep = "_")
  )

write.csv(pop_2lt_pearson, "Inversions/pop_2lt_pearson_correlations.csv", row.names = F)


pop_lat_pearson_facet <- bind_rows(pop_list) %>%
  mutate(Label = paste(Trait, Sex, sep = "_")) %>%
  ggplot(aes(x = In2Lt, y = Estimate)) +
  geom_point(aes(color = Population), size = 2) +
  colScale +
  facet_wrap(Label ~., scales = "free", ncol = 7) +
  geom_smooth(method = "lm", se = F, color = "black", size = 0.5) +
  stat_cor(method = "pearson", label.x.npc = 0, label.y.npc = 0.05, size = 3) +
  labs(
    title = "Pearson correlations with In(2L)t frequencies - Population level",
    x = "In(2L)t frequency",
    y = "Population estimates"
  ) +
  theme_bw(14)

ggsave(
  pop_lat_pearson_facet,
  filename = "Inversions/pop_2lt_pearson_correlations_facets.png",
  height = 10, width = 14, dpi = 120
)