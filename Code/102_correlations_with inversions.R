
##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(MetBrewer)
library(lme4)
library(afex)


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


write.csv(invs_pop,
  file = "Inversions/inversion_frequencies_pops.csv",
  row.names = FALSE
)



##### pop level

pop_list <- filter(pop_comp, Trait != "Dia_lmer") %>%
  dplyr::select(Trait, Sex, Population, Estimate) %>%
  inner_join(invs_pop) %>%
  inner_join(pops$by_lat) %>%
  group_split(Trait, Sex)



### 2Lt


pop_2lt_pearson <- list()
for (i in 1:length(pop_list)) {
  cortest <- cor.test(pop_list[[i]]$Estimate, pop_list[[i]]$In2Lt)
  pop_2lt_pearson[[i]] <- data.frame(
    Trait = unique(pop_list[[i]]$Trait),
    Sex = unique(pop_list[[i]]$Sex),
    R = cortest$estimate,
    P = cortest$p.value,
    P_bonf = cortest$p.value * length(pop_list)
  )
}

pop_2lt_pearson <- bind_rows(pop_2lt_pearson) %>%
  mutate(
    P_bonf = ifelse(P_bonf > 1, 1, P_bonf),
    Method = "pearson",
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






### 3RP



pop_3rp_pearson <- list()
for (i in 1:length(pop_list)) {
  cortest <- cor.test(pop_list[[i]]$Estimate, pop_list[[i]]$In3RP)
  pop_3rp_pearson[[i]] <- data.frame(
    Trait = unique(pop_list[[i]]$Trait),
    Sex = unique(pop_list[[i]]$Sex),
    R = cortest$estimate,
    P = cortest$p.value,
    P_bonf = cortest$p.value * length(pop_list)
  )
}

pop_3rp_pearson <- bind_rows(pop_3rp_pearson) %>%
  mutate(
    P_bonf = ifelse(P_bonf > 1, 1, P_bonf),
    Method = "pearson",
    Sex = factor(Sex, levels = c("F", "M", "B", "NA")),
    Label = paste(Trait, Sex, sep = "_")
  )


write.csv(pop_3rp_pearson, "Inversions/pop_3RP_pearson_correlations.csv", row.names = F)


pop_lat_pearson_facet <- bind_rows(pop_list) %>%
  mutate(Label = paste(Trait, Sex, sep = "_")) %>%
  ggplot(aes(x = In3RP, y = Estimate)) +
  geom_point(aes(color = Population), size = 2) +
  colScale +
  facet_wrap(Label ~., scales = "free", ncol = 7) +
  geom_smooth(method = "lm", se = F, color = "black", size = 0.5) +
  stat_cor(method = "pearson", label.x.npc = 0, label.y.npc = 0.05, size = 3) +
  labs(
    title = "Pearson correlations with In(3R)P frequencies - Population level",
    x = "In(3R)P frequency",
    y = "Population estimates"
  ) +
  theme_bw(14)

ggsave(
  pop_lat_pearson_facet,
  filename = "Inversions/pop_3RP_pearson_correlations_facets.png",
  height = 10, width = 14, dpi = 120
)






### 3LP



pop_3lp_pearson <- list()
for (i in 1:length(pop_list)) {
  cortest <- cor.test(pop_list[[i]]$Estimate, pop_list[[i]]$In3LP)
  pop_3lp_pearson[[i]] <- data.frame(
    Trait = unique(pop_list[[i]]$Trait),
    Sex = unique(pop_list[[i]]$Sex),
    R = cortest$estimate,
    P = cortest$p.value,
    P_bonf = cortest$p.value * length(pop_list)
  )
}

pop_3lp_pearson <- bind_rows(pop_3lp_pearson) %>%
  mutate(
    P_bonf = ifelse(P_bonf > 1, 1, P_bonf),
    Method = "pearson",
    Sex = factor(Sex, levels = c("F", "M", "B", "NA")),
    Label = paste(Trait, Sex, sep = "_")
  )


write.csv(pop_3lp_pearson, "Inversions/pop_3LP_pearson_correlations.csv", row.names = F)


pop_lat_pearson_facet <- bind_rows(pop_list) %>%
  mutate(Label = paste(Trait, Sex, sep = "_")) %>%
  ggplot(aes(x = In3LP, y = Estimate)) +
  geom_point(aes(color = Population), size = 2) +
  colScale +
  facet_wrap(Label ~., scales = "free", ncol = 7) +
  geom_smooth(method = "lm", se = F, color = "black", size = 0.5) +
  stat_cor(method = "pearson", label.x.npc = 0, label.y.npc = 0.05, size = 3) +
  labs(
    title = "Pearson correlations with In(3L)P frequencies - Population level",
    x = "In(3L)P frequency",
    y = "Population estimates"
  ) +
  theme_bw(14)

ggsave(
  pop_lat_pearson_facet,
  filename = "Inversions/pop_3LP_pearson_correlations_facets.png",
  height = 10, width = 14, dpi = 120
)





### 2RNs



pop_2rns_pearson <- list()
for (i in 1:length(pop_list)) {
  cortest <- cor.test(pop_list[[i]]$Estimate, pop_list[[i]]$In2RNs)
  pop_2rns_pearson[[i]] <- data.frame(
    Trait = unique(pop_list[[i]]$Trait),
    Sex = unique(pop_list[[i]]$Sex),
    R = cortest$estimate,
    P = cortest$p.value,
    P_bonf = cortest$p.value * length(pop_list)
  )
}

pop_2rns_pearson <- bind_rows(pop_2rns_pearson) %>%
  mutate(
    P_bonf = ifelse(P_bonf > 1, 1, P_bonf),
    Method = "pearson",
    Sex = factor(Sex, levels = c("F", "M", "B", "NA")),
    Label = paste(Trait, Sex, sep = "_")
  )


write.csv(pop_2rns_pearson, "Inversions/pop_2RNs_pearson_correlations.csv", row.names = F)


pop_lat_pearson_facet <- bind_rows(pop_list) %>%
  mutate(Label = paste(Trait, Sex, sep = "_")) %>%
  ggplot(aes(x = In2RNs, y = Estimate)) +
  geom_point(aes(color = Population), size = 2) +
  colScale +
  facet_wrap(Label ~., scales = "free", ncol = 7) +
  geom_smooth(method = "lm", se = F, color = "black", size = 0.5) +
  stat_cor(method = "pearson", label.x.npc = 0, label.y.npc = 0.05, size = 3) +
  labs(
    title = "Pearson correlations with In(2R)Ns frequencies - Population level",
    x = "In(2R)Ns frequency",
    y = "Population estimates"
  ) +
  theme_bw(14)

ggsave(
  pop_lat_pearson_facet,
  filename = "Inversions/pop_2RNs_pearson_correlations_facets.png",
  height = 10, width = 14, dpi = 120
)




### 3Rmo



pop_3rmo_pearson <- list()
for (i in 1:length(pop_list)) {
  cortest <- cor.test(pop_list[[i]]$Estimate, pop_list[[i]]$In3Rmo)
  pop_3rmo_pearson[[i]] <- data.frame(
    Trait = unique(pop_list[[i]]$Trait),
    Sex = unique(pop_list[[i]]$Sex),
    R = cortest$estimate,
    P = cortest$p.value,
    P_bonf = cortest$p.value * length(pop_list)
  )
}

pop_3rmo_pearson <- bind_rows(pop_3rmo_pearson) %>%
  mutate(
    P_bonf = ifelse(P_bonf > 1, 1, P_bonf),
    Method = "pearson",
    Sex = factor(Sex, levels = c("F", "M", "B", "NA")),
    Label = paste(Trait, Sex, sep = "_")
  )


write.csv(pop_3rmo_pearson, "Inversions/pop_3Rmo_pearson_correlations.csv", row.names = F)


pop_lat_pearson_facet <- bind_rows(pop_list) %>%
  mutate(Label = paste(Trait, Sex, sep = "_")) %>%
  ggplot(aes(x = In3Rmo, y = Estimate)) +
  geom_point(aes(color = Population), size = 2) +
  colScale +
  facet_wrap(Label ~., scales = "free", ncol = 7) +
  geom_smooth(method = "lm", se = F, color = "black", size = 0.5) +
  stat_cor(method = "pearson", label.x.npc = 0, label.y.npc = 0.05, size = 3) +
  labs(
    title = "Pearson correlations with In(3R)mo frequencies - Population level",
    x = "In(3R)mo frequency",
    y = "Population estimates"
  ) +
  theme_bw(14)

ggsave(
  pop_lat_pearson_facet,
  filename = "Inversions/pop_3Rmo_pearson_correlations_facets.png",
  height = 10, width = 14, dpi = 120
)





# line level GLMs


## four lines "RE9"  "UM22" "UM24" "AK19" have not been typed for inversions

line_invs <- inner_join(
  line_comp,
  invs %>% select(contains("In")) %>% rename(Line = line_id)
) %>%
  mutate(In2Lt = ifelse(In2Lt == 0.5 | In2Lt == 1, 1, 0)) %>%
  group_split(Trait, Sex)


# In3RP

line_3rp_lm <- list()
for (i in 1:length(line_invs)) {
  inv_lm <- lm(Estimate ~ as.factor(In3RP), data = line_invs[[i]])
  line_3rp_lm[[i]] <- data.frame(
    Trait = unique(line_invs[[i]]$Trait),
    Sex = unique(line_invs[[i]]$Sex),
    R = summary(inv_lm)$r.squared,
    P = summary(inv_lm)$coefficients[2,4],
    P_bonf = summary(inv_lm)$coefficients[2,4] * length(line_invs)
  )
}

line_3rp_lm <- bind_rows(line_3rp_lm) %>%
  mutate(
    P_bonf = ifelse(P_bonf > 1, 1, P_bonf),
    Sex = factor(Sex, levels = c("F", "M", "B", "NA")),
    Label = paste(Trait, Sex, sep = "_")
  )


write.csv(line_3rp_lm, "Inversions/line_3RP_lms.csv", row.names = F)



line_3rp_lm_facet <- bind_rows(line_invs) %>%
  mutate(Label = paste(Trait, Sex, sep = "_")) %>%
  ggplot(aes(x = as.factor(In3RP), y = Estimate, color = Population)) +
  geom_quasirandom(size = 1, alpha = 0.5, width = .2) +
  colScale +
  facet_wrap(Label ~., scales = "free", ncol = 7) +
  labs(
    title = "Effect of In(3R)P on measured traits - Line level",
    x = "In(3R)P",
    y = "Line estimates"
  ) +
  theme_bw(14)

ggsave(
  line_3rp_lm_facet,
  filename = "Inversions/line_3RP_lm_facets.png",
  height = 10, width = 14, dpi = 120
)






# In3Rmo

line_3rmo_lm <- list()
for (i in 1:length(line_invs)) {
  inv_lm <- lm(Estimate ~ as.factor(In3Rmo), data = line_invs[[i]])
  line_3rmo_lm[[i]] <- data.frame(
    Trait = unique(line_invs[[i]]$Trait),
    Sex = unique(line_invs[[i]]$Sex),
    R = summary(inv_lm)$r.squared,
    P = summary(inv_lm)$coefficients[2,4],
    P_bonf = summary(inv_lm)$coefficients[2,4] * length(line_invs)
  )
}

line_3rmo_lm <- bind_rows(line_3rmo_lm) %>%
  mutate(
    P_bonf = ifelse(P_bonf > 1, 1, P_bonf),
    Sex = factor(Sex, levels = c("F", "M", "B", "NA")),
    Label = paste(Trait, Sex, sep = "_")
  )

write.csv(line_3rmo_lm, "Inversions/line_3Rmo_lms.csv", row.names = F)


line_3rmo_lm_facet <- bind_rows(line_invs) %>%
  mutate(Label = paste(Trait, Sex, sep = "_")) %>%
  ggplot(aes(x = as.factor(In3Rmo), y = Estimate, color = Population)) +
  geom_quasirandom(size = 1, alpha = 0.5, width = .2) +
  colScale +
  facet_wrap(Label ~., scales = "free", ncol = 7) +
  labs(
    title = "Effect of In(3R)mo on measured traits - Line level",
    x = "In(3R)mo",
    y = "Line estimates"
  ) +
  theme_bw(14)

ggsave(
  line_3rmo_lm_facet,
  filename = "Inversions/line_3Rmo_lm_facets.png",
  height = 10, width = 14, dpi = 120
)



# In2Lt

line_2Lt_lm <- list()
for (i in 1:length(line_invs)) {
  inv_lm <- lm(Estimate ~ as.factor(In2Lt), data = line_invs[[i]])
  line_2Lt_lm[[i]] <- data.frame(
    Trait = unique(line_invs[[i]]$Trait),
    Sex = unique(line_invs[[i]]$Sex),
    R = summary(inv_lm)$r.squared,
    P = summary(inv_lm)$coefficients[2,4],
    P_bonf = summary(inv_lm)$coefficients[2,4] * length(line_invs)
  )
}

line_2Lt_lm <- bind_rows(line_2Lt_lm) %>%
  mutate(
    P_bonf = ifelse(P_bonf > 1, 1, P_bonf),
    Sex = factor(Sex, levels = c("F", "M", "B", "NA")),
    Label = paste(Trait, Sex, sep = "_")
  )


write.csv(line_2Lt_lm, "Inversions/line_2Lt_lms.csv", row.names = F)



line_2Lt_lm_facet <- bind_rows(line_invs) %>%
  mutate(Label = paste(Trait, Sex, sep = "_")) %>%
  ggplot(aes(x = as.factor(In2Lt), y = Estimate, color = Population)) +
  geom_quasirandom(size = 1, alpha = 0.5, width = .2) +
  colScale +
  facet_wrap(Label ~., scales = "free", ncol = 7) +
  labs(
    title = "Effect of In(2L)t on measured traits - Line level",
    x = "In(2L)t",
    y = "Line estimates"
  ) +
  theme_bw(14)

ggsave(
  line_2Lt_lm_facet,
  filename = "Inversions/line_2Lt_lm_facets.png",
  height = 10, width = 14, dpi = 120
)





# In3LP

line_3lp_lm <- list()
for (i in 1:length(line_invs)) {
  inv_lm <- lm(Estimate ~ as.factor(In3LP), data = line_invs[[i]])
  line_3lp_lm[[i]] <- data.frame(
    Trait = unique(line_invs[[i]]$Trait),
    Sex = unique(line_invs[[i]]$Sex),
    R = summary(inv_lm)$r.squared,
    P = summary(inv_lm)$coefficients[2,4],
    P_bonf = summary(inv_lm)$coefficients[2,4] * length(line_invs)
  )
}

line_3lp_lm <- bind_rows(line_3lp_lm) %>%
  mutate(
    P_bonf = ifelse(P_bonf > 1, 1, P_bonf),
    Sex = factor(Sex, levels = c("F", "M", "B", "NA")),
    Label = paste(Trait, Sex, sep = "_")
  )


write.csv(line_3lp_lm, "Inversions/line_3LP_lms.csv", row.names = F)



line_3lp_lm_facet <- bind_rows(line_invs) %>%
  mutate(Label = paste(Trait, Sex, sep = "_")) %>%
  ggplot(aes(x = as.factor(In3LP), y = Estimate, color = Population)) +
  geom_quasirandom(size = 1, alpha = 0.5, width = .2) +
  colScale +
  facet_wrap(Label ~., scales = "free", ncol = 7) +
  labs(
    title = "Effect of In(3L)P on measured traits - Line level",
    x = "In(3L)P",
    y = "Line estimates"
  ) +
  theme_bw(14)

ggsave(
  line_3lp_lm_facet,
  filename = "Inversions/line_3LP_lm_facets.png",
  height = 10, width = 14, dpi = 120
)




# In2Rns

line_2rns_lm <- list()
for (i in 1:length(line_invs)) {
  inv_lm <- lm(Estimate ~ as.factor(In2RNs), data = line_invs[[i]])
  line_2rns_lm[[i]] <- data.frame(
    Trait = unique(line_invs[[i]]$Trait),
    Sex = unique(line_invs[[i]]$Sex),
    R = summary(inv_lm)$r.squared,
    P = summary(inv_lm)$coefficients[2,4],
    P_bonf = summary(inv_lm)$coefficients[2,4] * length(line_invs)
  )
}

line_2rns_lm <- bind_rows(line_2rns_lm) %>%
  mutate(
    P_bonf = ifelse(P_bonf > 1, 1, P_bonf),
    Sex = factor(Sex, levels = c("F", "M", "B", "NA")),
    Label = paste(Trait, Sex, sep = "_")
  )


write.csv(line_2rns_lm, "Inversions/line_2RNs_lms.csv", row.names = F)



line_2rns_lm_facet <- bind_rows(line_invs) %>%
  mutate(Label = paste(Trait, Sex, sep = "_")) %>%
  ggplot(aes(x = as.factor(In2RNs), y = Estimate, color = Population)) +
  geom_quasirandom(size = 1, alpha = 0.5, width = .2) +
  colScale +
  facet_wrap(Label ~., scales = "free", ncol = 7) +
  labs(
    title = "Effect of In(2R)Ns on measured traits - Line level",
    x = "In(2R)Ns",
    y = "Line estimates"
  ) +
  theme_bw(14)

ggsave(
  line_2rns_lm_facet,
  filename = "Inversions/line_2RNs_lm_facets.png",
  height = 10, width = 14, dpi = 120
)


