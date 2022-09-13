


##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(ggbeeswarm)
library(ggdist)
library(ggpubr)
library(MetBrewer)
#library(lme4)
#library(lsmeans)
#library(afex)
#library(multcomp)



##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")


##### dir
dir.create("Diets")

##### source functions
source("Code/functions.R")

##### load data
droseu <- readRDS("Data/droseu_master_list_2022-05-02.rds")
diets <- read.csv("InfoTables/DrosEU_Diets_Sept22.csv")


######### DIET VARIATION ##########


##### prep the data for plotting
d <- bind_rows(filter(diets, !is.na(PC)) %>% mutate(Group = "All P/C"),
               filter(diets, PC < 1) %>% mutate(Group = "P/C < 1")) 

##### plot the different PC ratios and Diets, regardless of traits

### facet plot
p1 <- select(d, Lab, PC, Diet, Group) %>% distinct() %>%
  ggplot(aes(x = Diet, y = PC, color = Diet)) +
  facet_wrap(. ~ Group, scales = "free") +
  geom_boxplot(outlier.colour = NA) +
  geom_quasirandom(size = 4, alpha = 0.5, width = .3) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none") +
  labs(x = "Diet", y = "P/C", title = "Protein/Carbohydrate ratios (P/C)") +
  theme(plot.title = element_text(size=15))
ggsave(p1, filename = "Diets/DrosEU_Diets_PC_ratios_facets.pdf", width = 6, height = 6)


### same as above but different layout
p2 <- select(d, Lab, PC, Diet, Group) %>% filter(Group == "All P/C") %>% distinct() %>%
  ggplot(aes(x = Diet, y = PC, color = Diet)) +
  geom_boxplot(outlier.colour = NA) +
  geom_quasirandom(size = 4, alpha = 0.5, width = .3) +
  theme_bw(base_size = 18) + #theme_classic(base_size = 24) 
  theme(legend.position = "none") +
  labs(x = "Diet", y = "P/C", title = "All P/C") +
  theme(plot.title = element_text(size=15)) # size=22
p3 <- select(d, Lab, PC, Diet, Group) %>% filter(Group == "P/C < 1") %>% distinct() %>%
  ggplot(aes(x = Diet, y = PC, color = Diet)) +
  geom_boxplot(outlier.colour = NA) +
  geom_quasirandom(size = 4, alpha = 0.5, width = .3) +
  theme_bw(base_size = 18) + 
  theme(legend.position = "none") +
  labs(x = "Diet", y = "P/C", title = "P/C < 1") +
  theme(plot.title = element_text(size=15))
p23 <- ggarrange(p2, p3)
ggsave(p23, filename = "Diets/DrosEU_Diets_PC_ratios.pdf", width = 6, height = 6)



##### plot the different PC ratios and Traits


p4 <- ggplot(data = d, aes(x = PC, y = Trait_long, color = Diet)) +
  geom_point(size = 3, alpha = 0.5) +
  facet_grid(. ~ Group, scales = "free") +
  theme_bw(base_size = 18) +
  labs(y = "Traits", x = "P/C", title = "Protein/Carbohydrate (P/C) ratios") +
  theme(panel.grid.major.y = element_line(size = 0.5)) +
  theme(plot.title = element_text(size=15))
ggsave(p4, filename = "Diets/DrosEU_Diets_PC_ratios_traits_facets.pdf", width = 8, height = 6)



p5 <- ggplot(data = d, aes(x = PC, y = Trait_short, color = Diet)) +
  geom_point(size = 3, alpha = 0.5) +
  theme_bw(base_size = 18) +
  labs(y = "Traits", x = "P/C", title = "All ratios") +
  theme(panel.grid.major.y = element_line(size = 0.5)) +
  theme(plot.title = element_text(size=15))
p6 <- ggplot(data = d, aes(x = PC, y = Trait_short, color = Diet)) +
  geom_point(size = 3, alpha = 0.5) +
  theme_bw(base_size = 18) +
  labs(y = "Traits", x = "P/C", title = "Without extreme ratios") +
  theme(panel.grid.major.y = element_line(size = 0.5)) +
  theme(plot.title = element_text(size=15))
p56 <- ggarrange(p5, p6, common.legend = TRUE)
ggsave(p56, filename = "Diets/DrosEU_Diets_PC_ratios_traits.pdf", width = 6, height = 6)




######### EFFECT OF DIET ON TRAITS ##########

pop_coefs <- read.csv("LinearModelsPop/all_models_pop_estimates.csv") %>%
  inner_join(dplyr::select(diets, Lab, Diet, PC) %>% distinct())

line_coefs <- read.csv("LinearModelsPop/all_models_line_random_coefs.csv") %>%
  inner_join(dplyr::select(diets, Lab, Diet, PC) %>% distinct())

myColors <- met.brewer("Johnson", 9)
names(myColors) <- as.factor(c("AK", "GI", "KA", "MA", "MU", "RE", "UM", "VA", "YE"))
colScale <- scale_colour_manual(name = "Population", values = myColors)


##### TL

### pop level

tl_pop_coefs <- filter(pop_coefs, Trait == "TL")

lab_pc_pop <- dplyr::select(tl_pop_coefs, Lab, PC, Sex) %>% 
  distinct() %>% arrange(PC) %>% mutate(ypos = 0.97*min(tl_pop_coefs$Estimate))

tl_pop_coefs <- mutate(tl_pop_coefs, Lab = factor(Lab, levels = unique(lab_pc_pop$Lab)))

p7 <- ggplot(data = tl_pop_coefs) +
  geom_quasirandom(aes(x = Lab, y = Estimate, color = Population), size = 3, 
                   alpha = 1, width = .3) +
  facet_wrap(. ~ Sex) +
  theme_bw(base_size = 18) +
  colScale +
  labs(y = "Thorax Length", x = "Lab (ordered by increasing P/C)", title = "Lab and Diet effects on Thorax Length - Pop level") +
  theme(panel.grid.major.y = element_line(size = 0.5)) +
  theme(plot.title = element_text(size=15)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  geom_text(data = lab_pc_pop, aes(x = Lab, y = ypos, label = PC)) +
  expand_limits(y=unique(lab_pc_pop$ypos))
ggsave(p7, filename = "Diets/DrosEU_Diets_PC_ratios_TL_pop_facets.pdf", width = 7, height = 6)


### line level

tl_line_coefs <- filter(line_coefs, Trait == "TL")

lab_pc_line <- dplyr::select(tl_line_coefs, Lab, PC, Sex) %>% 
  distinct() %>% arrange(PC) %>% mutate(ypos = 0.97*min(tl_line_coefs$Coef))

tl_line_coefs <- mutate(tl_line_coefs, Lab = factor(Lab, levels = unique(lab_pc_line$Lab)))

p8 <- ggplot(data = tl_line_coefs) +
  geom_quasirandom(aes(x = Lab, y = Coef, color = Population), size = 1, 
                   alpha = 1, width = .3, pch = 1) +
  facet_wrap(. ~ Sex) +
  theme_bw(base_size = 18) +
  colScale +
  labs(y = "Thorax Length", x = "Lab (ordered by increasing P/C)", title = "Lab and Diet effects on Thorax Length - Line level") +
  theme(panel.grid.major.y = element_line(size = 0.5)) +
  theme(plot.title = element_text(size=15)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  geom_text(data = lab_pc_line, aes(x = Lab, y = ypos, label = PC)) +
  expand_limits(y=unique(lab_pc_line$ypos))
ggsave(p8, filename = "Diets/DrosEU_Diets_PC_ratios_TL_line_facets.pdf", width = 7, height = 6)




##### Viab

### pop level

via_pop_coefs <- filter(pop_coefs, Trait == "Via")

lab_pc_pop <- dplyr::select(via_pop_coefs, Lab, PC, Sex) %>% 
  distinct() %>% arrange(PC) %>% mutate(ypos = 0.97*min(via_pop_coefs$Estimate))

via_pop_coefs <- mutate(via_pop_coefs, Lab = factor(Lab, levels = unique(lab_pc_pop$Lab)))

p9 <- ggplot(data = via_pop_coefs) +
  geom_quasirandom(aes(x = Lab, y = Estimate, color = Population), size = 3, 
                   alpha = 1, width = .3) +
  theme_bw(base_size = 18) +
  facet_wrap(Trait ~ .) +
  colScale +
  labs(y = "Viability", x = "Lab (ordered by increasing P/C)", title = "Lab and Diet effects on Viability - Pop level") +
  theme(panel.grid.major.y = element_line(size = 0.5)) +
  theme(plot.title = element_text(size=15)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  geom_text(data = lab_pc_pop, aes(x = Lab, y = ypos, label = PC)) +
  expand_limits(y=unique(lab_pc_pop$ypos))
ggsave(p9, filename = "Diets/DrosEU_Diets_PC_ratios_Via_pop_facets.pdf", width = 7, height = 6)


### line level

via_line_coefs <- filter(line_coefs, Trait == "Via")

lab_pc_line <- dplyr::select(via_line_coefs, Lab, PC, Sex) %>% 
  distinct() %>% arrange(PC) %>% mutate(ypos = 0.97*min(via_line_coefs$Coef))

via_line_coefs <- mutate(via_line_coefs, Lab = factor(Lab, levels = unique(lab_pc_line$Lab)))

p10 <- ggplot(data = via_line_coefs) +
  geom_quasirandom(aes(x = Lab, y = Coef, color = Population), size = 1, 
                   alpha = 1, width = .3, pch = 1) +
  facet_wrap(Trait ~ .) +
  theme_bw(base_size = 18) +
  colScale +
  labs(y = "Viability", x = "Lab (ordered by increasing P/C)", title = "Lab and Diet effects on Viability - Line level") +
  theme(panel.grid.major.y = element_line(size = 0.5)) +
  theme(plot.title = element_text(size=15)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_text(data = lab_pc_line, aes(x = Lab, y = ypos, label = PC)) +
  expand_limits(y=unique(lab_pc_line$ypos))
ggsave(p10, filename = "Diets/DrosEU_Diets_PC_ratios_Via_line_facets.pdf", width = 7, height = 6)





##### DTA

### pop level

dta_pop_coefs <- filter(pop_coefs, Trait == "DT_A")

lab_pc_pop <- dplyr::select(dta_pop_coefs, Lab, PC, Sex) %>% 
  distinct() %>% arrange(PC) %>% mutate(ypos = 0.97*min(dta_pop_coefs$Estimate))

dta_pop_coefs <- mutate(dta_pop_coefs, Lab = factor(Lab, levels = unique(lab_pc_pop$Lab)))

p11 <- ggplot(data = dta_pop_coefs) +
  geom_quasirandom(aes(x = Lab, y = Estimate, color = Population), size = 3, 
                   alpha = 1, width = .3) +
  theme_bw(base_size = 18) +
  facet_wrap(Sex ~ .) +
  colScale +
  labs(y = "Egg-to-adult development time", x = "Lab (ordered by increasing P/C)", title = "Lab and Diet effects on dev. time - Pop level") +
  theme(panel.grid.major.y = element_line(size = 0.5)) +
  theme(plot.title = element_text(size=15)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  geom_text(data = lab_pc_pop, aes(x = Lab, y = ypos, label = round(PC, 2))) +
  expand_limits(y=unique(lab_pc_pop$ypos))
ggsave(p11, filename = "Diets/DrosEU_Diets_PC_ratios_DT_A_pop_facets.pdf", width = 7, height = 6)


### line level

dta_line_coefs <- filter(line_coefs, Trait == "DT_A")

lab_pc_line <- dplyr::select(dta_line_coefs, Lab, PC, Sex) %>% 
  distinct() %>% arrange(PC) %>% mutate(ypos = 0.97*min(dta_line_coefs$Coef))

dta_line_coefs <- mutate(dta_line_coefs, Lab = factor(Lab, levels = unique(lab_pc_line$Lab)))

p12 <- ggplot(data = dta_line_coefs) +
  geom_quasirandom(aes(x = Lab, y = Coef, color = Population), size = 1, 
                   alpha = 1, width = .3, pch = 1) +
  facet_wrap(Sex ~ .) +
  theme_bw(base_size = 18) +
  colScale +
  labs(y = "Viability", x = "Lab (ordered by increasing P/C)", title = "Lab and Diet effects  on dev. time - Line level") +
  theme(panel.grid.major.y = element_line(size = 0.5)) +
  theme(plot.title = element_text(size=15)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_text(data = lab_pc_line, aes(x = Lab, y = ypos, label = round(PC, 2))) +
  expand_limits(y=unique(lab_pc_line$ypos))
ggsave(p12, filename = "Diets/DrosEU_Diets_PC_ratios_DT_A_line_facets.pdf", width = 7, height = 6)




pop_coefs <- read.csv("LinearModelsPop/all_models_pop_estimates.csv") %>%
  filter(!(Trait == "Dia" & Model != "glmer")) %>%
  inner_join(dplyr::select(diets, Lab, Diet, PC) %>% distinct())



pop_coefs_list <- group_split(pop_coefs, Trait)

for (i in 1:length(pop_coefs_list)) {
  
  trait <- pop_coefs_list[[i]] %>% filter(!is.na(PC))
    
  if (length(unique(trait$Lab)) >= 2) {
    
    lab_pc <- dplyr::select(trait, Lab, PC, Sex) %>% 
      distinct() %>% arrange(PC) %>% mutate(ypos = 0.95*min(trait$Estimate))
    
    trait <- mutate(trait, Lab = factor(Lab, levels = unique(lab_pc$Lab)))
    
    p <- ggplot(data = trait) +
      geom_quasirandom(aes(x = Lab, y = Estimate, color = Population), size = 3, width = .3) +
      facet_wrap(Sex ~ .) +
      theme_bw(base_size = 18) +
      colScale +
      labs(y = "Population estimate", x = "Lab (ordered by increasing P/C)", title = paste("Lab and Diet effects on", unique(trait$Trait))) +
      theme(panel.grid.major.y = element_line(size = 0.5)) +
      theme(plot.title = element_text(size=15)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      geom_text(data = lab_pc, aes(x = Lab, y = ypos, label = round(PC, 2))) +
      expand_limits(y=unique(lab_pc$ypos))
    
    ggsave(p, filename = paste("Diets/DrosEU_Diets_PC_ratios_", unique(trait$Trait), "_pop_facets.pdf", sep = ""), width = 7, height = 6) 
    
  }
  
}
    




library(metafor)


estimates <- readRDS("LinearModelsPop/all_models_pop_estimates_list.rds")

### dw

# get effects
dw_effects <- makeEffects(estimates$dw_lmer)

# add diet info
dw_effects <- inner_join(dw_effects, select(diets, Lab, PC) %>% distinct())

dw_F_meta_reg <- rma(yi = Y, sei = SE, data = filter(dw_effects, Sex == "F"), method = "ML", mods = ~ Population)

dw_F_meta_reg_full <- rma(yi = Y, sei = SE, data = filter(dw_effects, Sex == "F"), method = "ML", mods = ~ Population + Lab)

anova(dw_F_meta_reg, dw_F_meta_reg_full)


### via

# get effects
via_effects <- makeEffects(estimates$via_lmer)

# add diet info
via_effects <- inner_join(via_effects, select(diets, Lab, PC) %>% distinct())

via_F_meta_reg <- rma(yi = Y, sei = SE, data = via_effects, method = "ML", mods = ~ Population)

via_F_meta_reg_full <- rma(yi = Y, sei = SE, data = via_effects, method = "ML", mods = ~ Lab)

anova(via_F_meta_reg, via_F_meta_reg_full)



# get effects
sr_effects <- makeEffects(estimates$sr_lmer)

# add diet info
sr_effects <- inner_join(sr_effects, select(diets, Lab, PC) %>% distinct())

sr_M_meta_reg <- rma(yi = Y, sei = SE, data = filter(sr_effects, Sex == "M"), method = "ML", mods = ~ Population)

sr_M_meta_reg_full <- rma(yi = Y, sei = SE, data = filter(sr_effects, Sex == "M"), method = "ML", mods = ~ PC)

anova(via_F_meta_reg, via_F_meta_reg_full)



# get effects
dt_effects <- makeEffects(estimates$dt_lmer)

# add diet info
dt_effects <- inner_join(dt_effects, select(diets, Lab, PC) %>% distinct())

dt_M_meta_reg <- rma(yi = Y, sei = SE, data = filter(sr_effects, Sex == "M"), method = "ML", mods = ~ Population)

sr_M_meta_reg_full <- rma(yi = Y, sei = SE, data = filter(sr_effects, Sex == "M"), method = "ML", mods = ~ PC)

anova(via_F_meta_reg, via_F_meta_reg_full)


    
    
    

