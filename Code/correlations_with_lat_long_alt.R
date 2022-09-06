
###############################################################
###################### SOME CORRELATIONS ######################
###############################################################

### IN PROGRESS

##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(ggrepel)
library(ggpubr)

##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")



#### fix data
pops_estimates <- read.csv("LinearModelsPop/all_models_pop_estimates.csv") %>% mutate(Model = paste(Model, Predictor, sep = "_")) %>% select(-Predictor)
 
lines_re <- read.csv("LinearModelsPop/all_models_line_random_effects.csv") %>%
  rename(Re = Estimate, SE_re = SE)

lines_coefs <- inner_join(pops_estimates, lines_re) %>%
  mutate(SE = sqrt(SE^2 + SE_re^2)) %>% # check this one
  mutate(Coef = Estimate + Re) %>%
  select(Model, Trait, Lab, Sex, Population, Line, Coef, SE)
  
write.csv(lines_coefs, "LinearModelsPop/all_models_line_random_coefs.csv", row.names = F)




# where is LA absphase in lines re. It is a lm, Lines not included in there.





 
pops_estimates$Trait[pops_estimates$Trait == "Dia"] <- paste(pops_estimates$Trait[pops_estimates$Trait == "Dia"], pops_estimates$Model[pops_estimates$Trait == "Dia"], sep = "_")

pops_estimates$Trait <- gsub("Dia_lm", "Dia_lmer", pops_estimates$Trait)
pops_estimates$Trait <- gsub("Dia_lmerer", "Dia_lmer", pops_estimates$Trait)


##### load data
pops <- readRDS("InfoTables/DrosEU_Populations.rds")
droseu <- readRDS("Data/droseu_master_list_2022-05-02.rds")
pop_comp <- read.csv("MetaAnalyses/all_metas_pop_compound_estimates.csv")
lines_comp <- readRDS("LinearModelsPop/all_models_line_compound_random_effects_list.rds")



##### pop level

pop_geo <- inner_join(select(pop_comp, Trait, Sex, Population, Mstar), pops$by_lat) %>%
  group_split(Trait, Sex)

pearson_cor <- list()
for (i in 1:length(pop_geo)) {
  cortest <- cor.test(pop_geo[[i]]$Mstar, pop_geo[[i]]$Latitude)
  pearson_cor[[i]] <- data.frame(Trait = unique(pop_geo[[i]]$Trait), 
                                 Sex = unique(pop_geo[[i]]$Sex),
                                 R = cortest$estimate,
                                 P = cortest$p.value,
                                 Method = "pearson")
}
pearson_cor <- bind_rows(pearson_cor)


ggplot(data = pearson_cor, aes(x = R, y = -log10(P))) +
  geom_point(aes(color = Sex), size = 3) +
  theme_classic() +
  geom_hline(yintercept = -log10(0.05/25), linetype = 2, size = 0.15) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, size = 0.15) +
  labs(title = "Pearson correlations between Population and Latitude", colour = "Sex", x = "Pearson's R", y = "-log10(Pvalue)") +
  geom_label_repel(aes(x = R, y = -log10(P), label = Trait),
                   size = 2.5, box.padding = 0.25, label.padding = 0.15, point.padding = 0,
                   segment.color = 'grey50', max.overlaps = 60, min.segment.length = 0,
                   seed = 1, force = 5)
  

ggsave(pval_slope_lat_lon_eu, 
       file = "", 
       width = 6, height = 6, dpi = 300)







spearman_cor <- list()
for (i in 1:length(pop_geo)) {
  cortest <- cor.test(rank(pop_geo[[i]]$Mstar), rank(pop_geo[[i]]$Latitude))
  spearman_cor[[i]] <- data.frame(Trait = unique(pop_geo[[i]]$Trait), 
                                 Sex = unique(pop_geo[[i]]$Sex),
                                 R = cortest$estimate,
                                 P = cortest$p.value,
                                 Method = "spearman")
}
spearman_cor <- bind_rows(spearman_cor)




##### line level

lines_geo <- bind_rows(lines_comp) %>%
  select(Trait, Sex, Population, Line, Value) 

%>%
  inner_join(pops$by_lat) %>%
  group_split(Trait, Sex)




ggplot(data = filter(lines_geo, Trait == "TL" & Sex == "F"), aes(x = Population, y = Value)) +
  geom_point(aes(color = Population)) +
  theme_classic() +
  geom_hline(yintercept = -log10(0.05/25), linetype = 2, size = 0.15) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, size = 0.15) +
  labs(title = "Pearson correlations between Population and Latitude", colour = "Sex", x = "Pearson's R", y = "-log10(Pvalue)") +
  geom_label_repel(aes(x = R, y = -log10(P), label = Trait),
                   size = 2.5, box.padding = 0.25, label.padding = 0.15, point.padding = 0,
                   segment.color = 'grey50', max.overlaps = 60, min.segment.length = 0,
                   seed = 1, force = 5)


a <- filter(lines_geo, Trait == "TL" & Sex == "F") %>%
  group_by(Population) %>% summarise(m = mean(Value)) %>%
  inner_join(pops$by_lat) 

%>%
  ggplot(aes(x = Latitude, y = m)) +
  geom_point(aes(color = Population))
  



pearson_cor_lines <- list()
for (i in 1:length(lines_geo)) {
  cortest <- cor.test(lines_geo[[i]]$Value, lines_geo[[i]]$Latitude)
  pearson_cor_lines[[i]] <- data.frame(Trait = unique(lines_geo[[i]]$Trait), 
                                 Sex = unique(lines_geo[[i]]$Sex),
                                 R = cortest$estimate,
                                 P = cortest$p.value,
                                 Method = "pearson")
}
pearson_cor_lines <- bind_rows(pearson_cor_lines)



spearman_cor <- list()
for (i in 1:length(pop_geo)) {
  cortest <- cor.test(rank(pop_geo[[i]]$Mstar), rank(pop_geo[[i]]$Latitude))
  spearman_cor[[i]] <- data.frame(Trait = unique(pop_geo[[i]]$Trait), 
                                  Sex = unique(pop_geo[[i]]$Sex),
                                  R = cortest$estimate,
                                  P = cortest$p.value,
                                  Method = "spearman")
}
spearman_cor <- bind_rows(spearman_cor)





a <- lmer(TL_micrometers ~ Population + (1|Line), data = filter(droseu$tl, Supervisor.PI == "Kozeretska" & Sex == "F"))


e <- extract_random_effects(a, re = "Line") %>% 
  separate(group, c("Line", "Population"), ":") %>%
  rename(Estimate = value, SE = se) 

b <- lmer(TL_micrometers ~ Population + (1|Population:Line), data = filter(droseu$tl, Supervisor.PI == "Kozeretska" & Sex == "F"))

e <- extract_random_effects(b, re = "Population:Line") %>% 
  separate(group, c("Population", "Line"), ":") %>%
  rename(Estimate = value, SE = se) 

inner_join(e, filter(droseu$tl, Supervisor.PI == "Kozeretska" & Sex == "F") %>% group_by(Line) %>% summarize(m = mean(TL_micrometers))) %>%
  ggplot(aes(x = Estimate, y = m)) +
  geom_point(aes(color = Population))


d <- lmer(TL_micrometers ~ 1 + (1|Population:Line), data = filter(droseu$tl, Supervisor.PI == "Schmidt" & Sex == "F"))

e <- extract_random_effects(d, re = "Population:Line") %>% 
  separate(group, c("Population", "Line"), ":") %>%
  rename(Estimate = value, SE = se) 

inner_join(e, filter(droseu$tl, Supervisor.PI == "Schmidt" & Sex == "F") %>% group_by(Line) %>% summarize(m = mean(TL_micrometers))) %>%
  ggplot(aes(x = Estimate, y = m)) +
  geom_point(aes(color = Population))

inner_join(e, pops$by_lat) %>%
  ggplot(aes(x = Estimate, y = Latitude)) +
  geom_point(aes(color = Population))


p <- inner_join(e, pops$by_lat)

dd <- lme4::lmer(TL_micrometers ~ Population + (1|Population:Line), data = filter(droseu$tl, Supervisor.PI == "Schmidt" & Sex == "F"))


ee <- extract_random_coefs(dd, re = "Population:Line") %>% 
  separate(group, c("Population", "Line"), ":") %>%
  rename(Estimate = value, SE = se) 

inner_join(ee, pops$by_lat) %>%
  ggplot(aes(x = Estimate, y = Latitude)) +
  geom_point(aes(color = Population))




eee <- extract_random_effects(dd, re = "Population:Line") %>% 
  separate(group, c("Population", "Line"), ":") %>%
  rename(Estimate = value, SE = se) 

inner_join(eee, pops$by_lat) %>%
  ggplot(aes(x = Estimate, y = Latitude)) +
  geom_point(aes(color = Population))

ggplot(aes(x = ee$Estimate, y = eee$Estimate)) +
  geom_point(aes(color = Population))


fe <- data.frame(stats::coef(summary(dd))) %>%
  dplyr::mutate(term = rownames(.)) %>%
  dplyr::select(term, dplyr::everything())

extract_fixed_effects(d)

fixed_effects  <- extract_fixed_effects(dd) %>%
  dplyr::rename(effect = term,
                se_fe = se,
                value_fe = value)

random_effects <- extract_random_effects(dd)
