


##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(ggbeeswarm)
library(ggdist)
library(ggpubr)
library(MetBrewer)
library(metafor)



##### set working directory #####
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")


##### make new ouput directory #####
dir.create("Diets", showWarnings = FALSE)

##### source functions #####
source("Code/functions.R")

##### load data #####
droseu <- readRDS("Data/droseu_master_list_2022-05-02.rds")
diets <- read.csv("InfoTables/DrosEU_Diets_Sept22.csv")
random_coefs <- read.csv("LinearModelsPop/all_models_line_random_coefs.csv")






######### EXPLORE DIET VARIATION ##########


##### prep the data for plotting

d <- bind_rows(
  filter(diets, !is.na(PC)) %>% mutate(Group = "All P/C"),
  filter(diets, PC < 1) %>% mutate(Group = "P/C < 1")
)

##### plot the different PC ratios and Diets, regardless of traits

### facet plot

p1 <- select(d, Lab, PC, Diet, Group) %>%
  distinct() %>%
  ggplot(aes(x = Diet, y = PC, color = Diet)) +
  facet_wrap(. ~ Group, scales = "free") +
  geom_boxplot(outlier.colour = NA) +
  geom_quasirandom(size = 4, alpha = 0.5, width = .3) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none") +
  labs(x = "Diet", y = "P/C", title = "Protein/Carbohydrate ratios (P/C)") +
  theme(plot.title = element_text(size = 15))

ggsave(p1,
  filename = "Diets/DrosEU_Diets_PC_ratios_facets.png",
  width = 6,
  height = 6,
  dpi = 120
)

### same as above but different layout

p2 <- select(d, Lab, PC, Diet, Group) %>%
  filter(Group == "All P/C") %>%
  distinct() %>%
  ggplot(aes(x = Diet, y = PC, color = Diet)) +
  geom_boxplot(outlier.colour = NA) +
  geom_quasirandom(size = 4, alpha = 0.5, width = .3) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none") +
  labs(x = "Diet", y = "P/C", title = "All P/C") +
  theme(plot.title = element_text(size = 15))

p3 <- select(d, Lab, PC, Diet, Group) %>%
  filter(Group == "P/C < 1") %>%
  distinct() %>%
  ggplot(aes(x = Diet, y = PC, color = Diet)) +
  geom_boxplot(outlier.colour = NA) +
  geom_quasirandom(size = 4, alpha = 0.5, width = .3) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none") +
  labs(x = "Diet", y = "P/C", title = "P/C < 1") +
  theme(plot.title = element_text(size = 15))

ggsave(ggarrange(p2, p3),
  filename = "Diets/DrosEU_Diets_PC_ratios.png",
  width = 6,
  height = 6,
  dpi = 120
)



##### plot the different PC ratios and Traits

### facet plot

p4 <- ggplot(data = d, aes(x = PC, y = Trait_long, color = Diet)) +
  geom_point(size = 3, alpha = 0.5) +
  facet_grid(. ~ Group, scales = "free") +
  theme_bw(base_size = 18) +
  labs(y = "Traits", x = "P/C", title = "Protein/Carbohydrate (P/C) ratios") +
  theme(panel.grid.major.y = element_line(size = 0.5)) +
  theme(plot.title = element_text(size = 15))

ggsave(p4,
  filename = "Diets/DrosEU_Diets_PC_ratios_traits_facets.png",
  width = 8,
  height = 6,
  dpi = 120
)

### same as above but different layout

p5 <- ggplot(data = d, aes(x = PC, y = Trait_short, color = Diet)) +
  geom_point(size = 3, alpha = 0.5) +
  theme_bw(base_size = 18) +
  labs(y = "Traits", x = "P/C", title = "All ratios") +
  theme(panel.grid.major.y = element_line(size = 0.5)) +
  theme(plot.title = element_text(size = 15))

p6 <- d %>%
  filter(Group == "P/C < 1") %>%
  distinct() %>%
  ggplot(aes(x = PC, y = Trait_short, color = Diet)) +
  geom_point(size = 3, alpha = 0.5) +
  theme_bw(base_size = 18) +
  labs(y = "Traits", x = "P/C", title = "Without extreme ratios") +
  theme(panel.grid.major.y = element_line(size = 0.5)) +
  theme(plot.title = element_text(size = 15))

  ggsave(ggarrange(p5, p6, common.legend = TRUE),
    filename = "Diets/DrosEU_Diets_PC_ratios_traits.png",
    width = 6,
    height = 6,
    dpi = 120
  )







######### IDENTIFY SUBSETS OF TRAITS WITH SIMILAR PCR RATIOS ########

# be consistent with Ewan's PCAs
# nine traits in common between males and females
c9 <- c("ccrt", "csm", "hsm", "dta", "dw", "ls", "sr", "tl", "wa")
# female traits
fmax <- c(c9, "dia", "fec", "pgm")
# female traits plus
fmax_plus <- c(fmax, "via")

# based on diets plots above, the below range covers most of the traits
diet_range <- c(0.094, 0.183)

# keep labs within diet range, remove partial labs and split by trait
diet_list <- filter(d, PC >= diet_range[1] & PC <= diet_range[2] &
  Trait_short %in% fmax_plus & !Lab %in% c("Posnien", "Ritchie")) %>%
  select(Lab, Trait_short, PC) %>%
  arrange(PC, Trait_short) %>%
  distinct() %>%
  group_split(Trait_short)

# randomly sample one lab per trait when there are more than one lab

sample_row <- function(x, n) {
  x[sample(nrow(x), n), ]
}

set.seed(1)
diet_sub <- lapply(diet_list, sample_row, 1) %>%
  bind_rows()

# add in fecundity and add a column for merging with random coefs

diet_plus_fec <- bind_rows(
  diet_sub,
  filter(d, Trait_short == "fec" & Lab == "Fricke") %>%
    select(Lab, Trait_short, PC) %>%
    distinct()
  ) %>%
  arrange(PC) %>%
  mutate(
    Trait_lab = paste(Trait_short, Lab, sep = "_"),
    Trait_lab = str_replace(Trait_lab, "wa", "wa_l"),
    Trait_lab = str_replace(Trait_lab, "dta", "dt_a"),
    Trait_lab = str_replace(Trait_lab, "pgm", "pgm_total")
  )

# keep random coefs for the selected labs, remove lmer Dia

rc_trait_lab <- random_coefs %>%
  mutate(Trait_lab = paste(tolower(Trait), Lab, sep = "_")) %>%
  filter(Trait_lab %in% diet_plus_fec$Trait_lab) %>%
  select(-Trait_lab) %>%
  filter(!(Trait == "Dia" & Model == "lmer_pop"))

write.csv(
  rc_trait_lab,
  "LinearModelsPop/all_models_line_random_coefs_similar_diet.csv",
  row.names = FALSE
)


######### IDENTIFY SUBSETS OF TRAITS WITH SIMILAR PCR RATIOS WITH MAXIMIZING NUMBER OF LABS ########

sel_labs <- filter(diets, PC_control2 == 1) %>%
  mutate(
    Trait_lab = paste(Trait_short, Lab, sep = "_"),
    Trait_lab = str_replace(Trait_lab, "wa", "wa_l"),
    Trait_lab = str_replace(Trait_lab, "dta", "dt_a"),
    Trait_lab = str_replace(Trait_lab, "pgm", "pgm_total")
  )

rc_trait_lab <- random_coefs %>%
  mutate(Trait_lab = paste(tolower(Trait), Lab, sep = "_")) %>%
  filter(Trait_lab %in% sel_labs$Trait_lab) %>%
  select(-Trait_lab) %>%
  filter(!(Trait == "Dia" & Model == "lmer_pop"))

write.csv(
  rc_trait_lab,
  "LinearModelsPop/all_models_line_random_coefs_similar_diet_v2.csv",
  row.names = FALSE
)




######### EFFECT OF LAB AND DIET ON TRAITS ##########

# get data in

pop_coefs <- read.csv("LinearModelsPop/all_models_pop_estimates.csv") %>%
  filter(!(Trait == "Dia" & Model != "glmer")) %>%
  filter(Lab != "Posnien") %>%
  inner_join(dplyr::select(diets, Lab, Diet, PC) %>%
    distinct())


# loop over all traits to plot population estimates for each lab

pop_coefs_list <- group_split(pop_coefs, Trait)

for (i in seq_len(length(pop_coefs_list))) {
  trait <- pop_coefs_list[[i]] %>% filter(!is.na(PC))

  if (length(unique(trait$Lab)) >= 2) {
    lab_pc <- dplyr::select(trait, Lab, PC, Sex) %>%
      distinct() %>%
      arrange(PC) %>%
      mutate(ypos = 0.95 * min(trait$Estimate))

    trait <- mutate(trait, Lab = factor(Lab, levels = unique(lab_pc$Lab)))

    p <- ggplot(data = trait) +
      geom_quasirandom(aes(x = Lab, y = Estimate, color = Population), size = 3, width = .3) +
      facet_wrap(Sex ~ .) +
      theme_bw(base_size = 18) +
      droseu_color_scale +
      labs(
        y = "Population estimate",
        x = "Lab (ordered by increasing P/C)",
        title = paste("Lab and Diet effects on", unique(trait$Trait))
        ) +
      theme(panel.grid.major.y = element_line(size = 0.5)) +
      theme(plot.title = element_text(size = 15)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      geom_text(data = lab_pc, aes(x = Lab, y = ypos, label = round(PC, 2))) +
      expand_limits(y = unique(lab_pc$ypos))

    ggsave(p,
      filename = paste("Diets/DrosEU_Diets_PC_ratios_", unique(trait$Trait),
        "_pop_facets.png", sep = ""),
      width = 7, height = 6, dpi = 120
    )
  }
}


# facet plot for pop estimates and labs

pop_estimates_labs_facets <- pop_coefs %>%
  arrange(PC) %>%
  mutate(
    Group = paste(Trait, Sex),
    Lab = str_replace(Lab, "StamenkovicRadak", "S-R"),
    Lab = factor(Lab, levels = (unique(Lab)))
    ) %>%
  ggplot(aes(x = Lab, y = Estimate, color = Population)) +
    droseu_color_scale +
    geom_quasirandom(size = 2, width = .3) +
    facet_wrap(Group ~ ., scale = "free", ncol = 7) +
    theme_bw(base_size = 14) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(x = "Lab (ordered by increasing P/C ratio)", y = "Population estimate")

ggsave(pop_estimates_labs_facets,
  filename = "Diets/DrosEU_Diets_lab_traits_pop_facets2.png",
  width = 14,
  height = 10,
  dpi = 120
)
 


# facet plot for pop estimates and PC ratios

pop_estimates_pc_facets <- pop_coefs %>%
  arrange(PC) %>%
  mutate(Group = paste(Trait, Sex)) %>%
  ggplot(aes(x = PC, y = Estimate, color = Lab)) +
  geom_point(size = 2) +
  facet_wrap(Group ~ ., scale = "free", ncol = 7) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = "Lab P/C ratio", y = "Population estimate")

ggsave(pop_estimates_pc_facets,
  filename = "Diets/DrosEU_Diets_PC_ratios_traits_pop_facets.png",
  width = 14,
  height = 10,
  dpi = 120
)




# facet plot with Lab and Diet at the same time


















####### META REGRESSION WITH LAB AND DIET #######

dir.create("MetaRegressionDiet", showWarnings = FALSE)

##### get all the line estimates, from all the traits

estimates <- readRDS("LinearModelsPop/all_models_pop_estimates.rds")

partial_labs <- "Posnien"
partial_traits <- c(
  "CCRT", "DT_P", "LA_NDlog2", "LA_Period", "LA_CircPhase",
  "LA_AbsPhase", "LA_Activity"
)

estimates_list <- estimates %>%
  filter(!Lab %in% partial_labs) %>%
  filter(!Trait %in% partial_traits) %>%
  filter(!(Trait == "Dia" & Model != "glmer")) %>%
  group_split(Trait, Sex)



traits <-  unique(estimates$Trait)

directories <- c(
  "ChillComa", "ColdShock", rep("DevelopmentTime", 2), "Diapause", "DryWeight",
  "Fecundity", "HeatShock", "Lifespan", rep("Locomotor", 5), rep("Pigmentation", 4),
  "Starvation", "ThoraxLength", "Viability", rep("WingArea", 2)
)

trait_dirs <- data.frame(Trait = traits, Dir = directories)




##### loop over estimates per trait and sex

r2 <- list()

for (i in seq_len(length(estimates_list))) {

  trait_sex <- estimates_list[[i]]
  trait <- unique(trait_sex$Trait)
  sex <- unique(trait_sex$Sex)
  model <- paste0(unique(trait_sex$Model)[1], "s")
  out_dir <- trait_dirs$Dir[trait_dirs$Trait == trait]
  out_name <- paste(trait, sex, model, "pop_meta_reg_diet.rds", sep = "_")
  out_path <- file.path("MetaRegressionDiet", out_dir, out_name)

  effects <- makeEffects(trait_sex)
  effects_diet <- inner_join(
    effects,
    select(diets, Lab, PC) %>%
      distinct()
  )

  meta_reg <- list()

  meta_reg$pop <- rma(
    yi = Y,
    sei = SE,
    data = effects_diet,
    method = "REML",
    mods = ~ Population
  )

  meta_reg$pc <- rma(
    yi = Y,
    sei = SE,
    data = effects_diet,
    method = "REML",
    mods = ~ PC
  )

  meta_reg$lab <- rma(
    yi = Y,
    sei = SE,
    data = effects_diet,
    method = "REML",
    mods = ~ Lab
  )

  meta_reg$pop_lab <- rma(
    yi = Y,
    sei = SE,
    data = effects_diet,
    method = "REML",
    mods = ~ Population + Lab
  )

  meta_reg$pop_pc <- rma(
    yi = Y,
    sei = SE,
    data = effects_diet,
    method = "REML",
    mods = ~ Population + PC
  )

  dir.create(file.path("MetaRegressionDiet", out_dir), showWarnings = FALSE)
  saveRDS(meta_reg, file = out_path)

  stats_output <- function(x) {
    data.frame(
      Trait = trait,
      Sex = sex,
      Moderator = as.character(x$formula.mods)[2],
      R2 = x$R2,
      P = x$QMp
    )
  }


  output <- lapply(meta_reg, stats_output) %>%
    bind_rows()

  write.csv(output, file = sub(".rds", "_r2.csv", out_path))

  r2[[i]] <- output

}

r2_all <- bind_rows(r2)

r2_facets <- r2_all %>%
  mutate(
    Group = paste(Trait, Sex),
    Moderator = str_replace(Moderator, "Population", "Pop")) %>%
  ggplot(aes(x = Moderator, y = R2, fill = Moderator)) +
    geom_bar(stat = "identity") +
    facet_wrap(Group ~ .) +
    theme_bw(14) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))



ggsave(r2_facets,
  filename = "MetaRegressionDiet/all_models_pop_meta_reg_diets_facets.png",
  width = 14,
  height = 10
)

