
###############################################################
###################### SOME CORRELATIONS ######################
###############################################################

### IN PROGRESS


# note to self: size patterns are clearly influenced by RE and AK, respectively the least and the most viable populations - could viability cause difference in larval density that would translate into size differences?

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
dir.create("GeoCorrelations")

##### load data
pops <- readRDS("InfoTables/DrosEU_Populations.rds")
droseu <- readRDS("Data/droseu_master_list_2022-05-02.rds")
pop_comp <- read.csv("MetaAnalyses/all_models_pop_meta_compound_estimates.csv")
line_comp <- read.csv("MetaAnalyses/all_models_line_meta_compound_random_coefs.csv")

##### define colors for plotting

myColors <- met.brewer("Johnson", 9)
names(myColors) <- as.factor(c("AK", "GI", "KA", "MA", "MU", "RE", "UM", "VA", "YE"))
colScale <- scale_colour_manual(name = "Population", values = myColors)



##### pop level

pop_geo_list <- filter(pop_comp, Trait != "Dia_lmer") %>%
  dplyr::select(Trait, Sex, Population, Estimate) %>%
  inner_join(pops$by_lat) %>%
  group_split(Trait, Sex)

# add traits measured only in 1 lab (traits that have not been through meta)
# also removing LA_AbsPhase for the time being since we don't have Line random coefs for this trait - for Pop and Line plots to be comparable

#single_lab <- filter(pop_estimates, Lab == "Tauber" | Trait == "DT_P") %>%
#  filter(Trait != "LA_AbsPhase") %>%
#  dplyr::select(Trait, Sex, Population, Estimate) %>% 
#  inner_join(pops$by_lat)

#pop_geo_list <- bind_rows(pop_geo, single_lab) %>%
#  group_split(Trait, Sex)
  

### latitude

pop_lat_pearson <- list()
for (i in 1:length(pop_geo_list)) {
  cortest <- cor.test(pop_geo_list[[i]]$Estimate, pop_geo_list[[i]]$Latitude)
  pop_lat_pearson[[i]] <- data.frame(Trait = unique(pop_geo_list[[i]]$Trait), 
                                 Sex = unique(pop_geo_list[[i]]$Sex),
                                 R = cortest$estimate,
                                 P = cortest$p.value,
                                 Method = "pearson")
}
pop_lat_pearson <- bind_rows(pop_lat_pearson) %>%
  mutate(Sex = factor(Sex, levels = c("F", "M", "B", "NA")),
         Label = paste(Trait, Sex, sep = "_"))

write.csv(pop_lat_pearson, "GeoCorrelations/pop_lat_pearson_correlations.csv", row.names = F)


pop_lat_pearson_plot <- ggplot(data = pop_lat_pearson, aes(x = R, y = -log10(P))) +
  geom_point(size = 4, alpha = 0.3, pch = 19) +
  theme_bw(16) +
  geom_hline(yintercept = -log10(0.05/nrow(pop_lat_pearson)), linetype = 2, size = 0.5, col = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, size = 0.5) +
  labs(title = "Pearson correlations with Latitude\nPopulation level", colour = "Sex", x = "Pearson's R", y = "-log10(Pvalue)") +
  geom_label_repel(data = mutate(pop_lat_pearson, Label = ifelse(P >= 0.05, NA, Label)), 
                          aes(x = R, y = -log10(P), label = Label),
                   size = 3, box.padding = 0.45, label.padding = 0.45, point.padding = 0,
                   segment.color = 'grey50', max.overlaps = 60, min.segment.length = 0,
                   seed = 1, force = 5)
  
ggsave(pop_lat_pearson_plot, filename = "GeoCorrelations/pop_lat_pearson_correlations.png", height = 7, width = 7, dpi = 120)


pop_lat_pearson_facet <- bind_rows(pop_geo_list) %>%
  mutate(Label = paste(Trait, Sex, sep = "_")) %>%
  ggplot(aes(x = Latitude, y = Estimate)) +
  geom_point(aes(color = Population), size = 2) +
  colScale +
  facet_wrap(Label ~., scales = "free", ncol = 7) +
  geom_smooth(method = "lm", se = F, color = "black", size = 0.5) +
  stat_cor(method = "pearson", label.x.npc = 0, label.y.npc = 0.05, size = 3) +
  labs(title = "Pearson correlations with Latitude - Population level", y = "Population estimates") +
  theme_bw(14)

ggsave(pop_lat_pearson_facet, filename = "GeoCorrelations/pop_lat_pearson_correlations_facets.png", height = 10, width = 14, dpi = 120)
  




### longitude

pop_lon_pearson <- list()
for (i in 1:length(pop_geo_list)) {
  cortest <- cor.test(pop_geo_list[[i]]$Estimate, pop_geo_list[[i]]$Longitude)
  pop_lon_pearson[[i]] <- data.frame(Trait = unique(pop_geo_list[[i]]$Trait), 
                                     Sex = unique(pop_geo_list[[i]]$Sex),
                                     R = cortest$estimate,
                                     P = cortest$p.value,
                                     Method = "pearson")
}
pop_lon_pearson <- bind_rows(pop_lon_pearson) %>%
  mutate(Sex = factor(Sex, levels = c("F", "M", "B", "NA")),
         Label = paste(Trait, Sex, sep = "_"))



write.csv(pop_lon_pearson, "GeoCorrelations/pop_lon_pearson_correlations.csv", row.names = F)


pop_lon_pearson_plot <- ggplot(data = pop_lon_pearson, aes(x = R, y = -log10(P))) +
  geom_point(size = 4, alpha = 0.3, pch = 19) +
  theme_bw(16) +
  geom_hline(yintercept = -log10(0.05/nrow(pop_lon_pearson)), linetype = 2, size = 0.5, col = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, size = 0.5) +
  labs(title = "Pearson correlations with Longitude\nPopulation level", colour = "Sex", x = "Pearson's R", y = "-log10(Pvalue)") +
  geom_label_repel(data = mutate(pop_lon_pearson, Label = ifelse(P >= 0.05, NA, Label)), 
                   aes(x = R, y = -log10(P), label = Label),
                   size = 3, box.padding = 0.45, label.padding = 0.45, point.padding = 0,
                   segment.color = 'grey50', max.overlaps = 60, min.segment.length = 0,
                   seed = 1, force = 5)

ggsave(pop_lon_pearson_plot, filename = "GeoCorrelations/pop_lon_pearson_correlations.png", height = 7, width = 7, dpi = 120)


pop_lon_pearson_facet <- bind_rows(pop_geo_list) %>%
  mutate(Label = paste(Trait, Sex, sep = "_")) %>%
  ggplot(aes(x = Longitude, y = Estimate)) +
  geom_point(aes(color = Population), size = 2) +
  colScale +
  facet_wrap(Label ~., scales = "free", ncol = 7) +
  geom_smooth(method = "lm", se = F, color = "black", size = 0.5) +
  stat_cor(method = "pearson", label.x.npc = 0, label.y.npc = 0.05, size = 3) +
  labs(title = "Pearson correlations with Longitude - Population level", y = "Population estimates") +
  theme_bw(14)

ggsave(pop_lon_pearson_facet, filename = "GeoCorrelations/pop_lon_pearson_correlations_facets.png", height = 10, width = 14, dpi = 120)













##### line level

line_geo_list <- inner_join(select(line_comp, Trait, Sex, Population, Estimate), pops$by_lat) %>%
  group_split(Trait, Sex)


### latitude

line_lat_pearson <- list()
for (i in 1:length(line_geo_list)) {
  cortest <- cor.test(line_geo_list[[i]]$Estimate, line_geo_list[[i]]$Latitude)
  line_lat_pearson[[i]] <- data.frame(Trait = unique(line_geo_list[[i]]$Trait), 
                                 Sex = unique(line_geo_list[[i]]$Sex),
                                 R = cortest$estimate,
                                 P = cortest$p.value,
                                 Method = "pearson")
}
line_lat_pearson <- bind_rows(line_lat_pearson) %>%
  mutate(Sex = factor(Sex, levels = c("F", "M", "B", "NA")),
         Label = paste(Trait, Sex, sep = "_"))


write.csv(line_lat_pearson, "GeoCorrelations/line_lat_pearson_correlations.csv", row.names = F)


line_lat_pearson_plot <- ggplot(data = line_lat_pearson, aes(x = R, y = -log10(P))) +
  geom_point(size = 4, alpha = 0.3, pch = 19) +
  theme_bw(16) +
  geom_hline(yintercept = -log10(0.05/nrow(line_lat_pearson)), linetype = 2, size = 0.5, col = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, size = 0.5) +
  labs(title = "Pearson correlations with Latitude\nLine level", colour = "Sex", x = "Pearson's R", y = "-log10(Pvalue)") +
  geom_label_repel(data = mutate(line_lat_pearson, Label = ifelse(P >= 0.05, NA, Label)), 
                   aes(x = R, y = -log10(P), label = Label),
                   size = 3, box.padding = 0.45, label.padding = 0.45, point.padding = 0,
                   segment.color = 'grey50', max.overlaps = 60, min.segment.length = 0,
                   seed = 1, force = 5)

ggsave(line_lat_pearson_plot, filename = "GeoCorrelations/line_lat_pearson_correlations.png", height = 7, width = 7, dpi = 120)




line_lat_pearson_facet <- bind_rows(line_geo_list) %>%
  mutate(Label = paste(Trait, Sex, sep = "_")) %>%
  ggplot(aes(x = Latitude, y = Estimate)) +
  geom_point(aes(color = Population), size = 2) +
  colScale +
  facet_wrap(Label ~., scales = "free", ncol = 7) +
  geom_smooth(method = "lm", se = F, color = "black", size = 0.5) +
  stat_cor(method = "pearson", label.x.npc = 0, label.y.npc = 0.05, size = 3) +
  labs(title = "Pearson correlations with Latitude - Line level", y = "Line random coefficients") +
  theme_bw(14)

ggsave(line_lat_pearson_facet, filename = "GeoCorrelations/line_lat_pearson_correlations_facets.png", height = 10, width = 14, dpi = 120)



### longitude

line_lon_pearson <- list()
for (i in 1:length(line_geo_list)) {
  cortest <- cor.test(line_geo_list[[i]]$Estimate, line_geo_list[[i]]$Longitude)
  line_lon_pearson[[i]] <- data.frame(Trait = unique(line_geo_list[[i]]$Trait), 
                                      Sex = unique(line_geo_list[[i]]$Sex),
                                      R = cortest$estimate,
                                      P = cortest$p.value,
                                      Method = "pearson")
}
line_lon_pearson <- bind_rows(line_lon_pearson) %>%
  mutate(Sex = factor(Sex, levels = c("F", "M", "B", "NA")),
         Label = paste(Trait, Sex, sep = "_"))


write.csv(line_lon_pearson, "GeoCorrelations/line_lon_pearson_correlations.csv", row.names = F)


line_lon_pearson_plot <- ggplot(data = line_lon_pearson, aes(x = R, y = -log10(P))) +
  geom_point(size = 4, alpha = 0.3, pch = 19) +
  theme_bw(16) +
  geom_hline(yintercept = -log10(0.05/nrow(line_lon_pearson)), linetype = 2, size = 0.5, col = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, size = 0.5) +
  labs(title = "Pearson correlations with Longitude\nLine level", colour = "Sex", x = "Pearson's R", y = "-log10(Pvalue)") +
  geom_label_repel(data = mutate(line_lon_pearson, Label = ifelse(P >= 0.05, NA, Label)), 
                   aes(x = R, y = -log10(P), label = Label),
                   size = 3, box.padding = 0.45, label.padding = 0.45, point.padding = 0,
                   segment.color = 'grey50', max.overlaps = 60, min.segment.length = 0,
                   seed = 1, force = 10)

ggsave(line_lon_pearson_plot, filename = "GeoCorrelations/line_lon_pearson_correlations.png", height = 7, width = 7, dpi = 120)




line_lon_pearson_facet <- bind_rows(line_geo_list) %>%
  mutate(Label = paste(Trait, Sex, sep = "_")) %>%
  ggplot(aes(x = Longitude, y = Estimate)) +
  geom_point(aes(color = Population), size = 2) +
  colScale +
  facet_wrap(Label ~., scales = "free", ncol = 7) +
  geom_smooth(method = "lm", se = F, color = "black", size = 0.5) +
  stat_cor(method = "pearson", label.x.npc = 0, label.y.npc = 0.05, size = 3) +
  labs(title = "Pearson correlations with Longitude - Line level", y = "Line random coefficients") +
  theme_bw(14)

ggsave(line_lon_pearson_facet, filename = "GeoCorrelations/line_lon_pearson_correlations_facets.png", height = 10, width = 14, dpi = 120)







#pl <- align_plots(pop_lat_pearson_plot, line_lat_pearson_plot, align="v")
#ggsave(pl[[1]], filename = "GeoCorrelations/pop_lat_pearson_correlations.pdf", height = 7, width = 7)
#ggsave(pl[[2]], filename = "GeoCorrelations/line_lat_pearson_correlations.pdf", height = 7, width = 7)







