



# need to take care of partial data and show a way to display them on the boxplots (DT_A)
# keep y axis scale consistent between males and females - to do



##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)

##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")

##### source functions
source("Code/functions.R")



##### load latest data
droseu <- readRDS("Data/droseu_master_list_2022-05-02.rds")
pops <- readRDS("InfoTables/DrosEU_Populations.rds")

##### load tukey stats
tuk_lmers <- readRDS("LinearModelsPop/all_lmers_pop_tukey_list.rds")
tuk_glmers <- readRDS("LinearModelsPop/all_glmers_pop_tukey_list.rds")

##### create output directory
lmer_plots_dir <- "LinearModelsPopPlots"
dir.create(lmer_plots_dir, showWarnings = FALSE)



##### fix and prepare the data



##### define trait and common variables

trait_var <- c(
  "CCRT_seconds", "ZT_hours_MESA", "ZT_hours_LSPR", "Period_MESA",
  "Period_LSPR", "Rhythmicity_LSPR_amplitude", "Rhythmicity_JTK_p_BH_corrected",
  "CSM_PropDead_ED", "Prop_Max_Stage9", "DT_EggAdult", "DT_EggPupa",
  "DW_micrograms", "NumberOfAdultsEclosed", "TimeDeath_min", "Period",
  "CircPhase", "AbsPhase", "ND", "Activity", "LSL_AgeAtDeath_days",
  "LSM_AgeAtDeath_days", "LSP_AgeAtDeath_days", "PercT4", "PercT5", "PercT6",
  "TotalPerc", "AgeAtDeath_hours", "TL_micrometers",
  "ProportionEggtoAdultSurvival", "CentroidSizeLeft_micrometers",
  "CentroidSizeRight_micrometers"
)

common_var <- c("Supervisor.PI", "Batch", "Population", "Line", "Sex")


##### add extra columns for easier batch processing and split data by single trait variables

droseu_long <- droseu
for (i in seq_along(droseu_long)){
  d <- droseu_long[[i]]
  if (!"Line" %in% colnames(d)) d$Line <- NA
  if (!"Batch" %in% colnames(d)) d$Batch <- NA
  if (!"Sex" %in% colnames(d)) d$Sex <- NA
  d <- d[, colnames(d) %in% c(common_var, trait_var)]
  d <- pivot_longer(d, cols = -all_of(common_var), names_to = "Trait_name", values_to = "Value")
  d$Trait <- names(droseu_long)[i]
  d <- split(d, d$Trait_name)
  droseu_long[[i]] <- d
}

droseu_long <- bind_rows(unlist(droseu_long, recursive = FALSE)) %>%
  rename(Trait_handle = Trait, Lab = Supervisor.PI)

trait_names <- read.csv("InfoTables/trait_names.csv") %>%
  dplyr::select(Trait_handle, Trait, Trait_name_raw, Legend, Title) %>%
  rename(Trait_name = Trait_name_raw)

droseu_long <- inner_join(
  droseu_long,
  trait_names
)

droseu_long <- inner_join(
  droseu_long,
  pops$by_lat
)


droseu_trait <- group_by(droseu_long, Trait, Trait_name, Sex) %>%
  group_split()


##### populations as main groups, labs as subgroups

for (i in seq_along(droseu_trait)) {

  title_text <- paste(unique(droseu_trait[[i]]$Title), unique(droseu_trait[[i]]$Sex), sep = " - ")
  title_text <- sub(" - F", " - Females", title_text)
  title_text <- sub(" - M", " - Males", title_text)

  out_file_png <- paste(
    "p_", unique(droseu_trait[[i]]$Trait), "_",
    unique(droseu_trait[[i]]$Sex), ".png", sep = ""
  )

  p <- ggplot(
    droseu_trait[[i]],
    aes(x = Country_code, y = Value, color = Lab, fill = Country_code)
  ) +
  geom_boxplot(outlier.shape = NA) +
  labs(
    x = "Country (by increasing latitude)",
    y = unique(droseu_trait[[i]]$Legend),
    title = title_text) +
  droseu_fill_scale_country +
  scale_color_grey(start = 0, end = 0) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 22),
    plot.title = element_text(size = 22)
  )

  ggsave(
    print(p),
    filename = file.path(lmer_plots_dir, out_file_png),
    width = 10, height = 6
  )

}


##### same as above but with flipped coords


for (i in seq_along(droseu_trait)) {

  title_text <- paste(unique(droseu_trait[[i]]$Title), unique(droseu_trait[[i]]$Sex), sep = " - ")
  title_text <- sub(" - F", " - Females", title_text)
  title_text <- sub(" - M", " - Males", title_text)
  title_text <- sub(" - NA", "", title_text)

  out_file_png <- paste(
    "p_", unique(droseu_trait[[i]]$Trait), "_",
    unique(droseu_trait[[i]]$Sex), "_flipped.png", sep = ""
  )

  p <- ggplot(
    droseu_trait[[i]],
    aes(x = Country_code, y = Value, color = Lab, fill = Country_code)
  ) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.9)) +
  labs(
    x = "Country (by latitude)",
    y = unique(droseu_trait[[i]]$Legend),
    title = title_text) +
  droseu_fill_scale_country +
    scale_color_grey(start = 0, end = 0) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 22),
      axis.title = element_text(size = 22),
      plot.title = element_text(size = 22)
  ) +
  coord_flip()

  ggsave(
    p,
    filename = file.path(lmer_plots_dir, out_file_png),
    width = 6, height = 8
  )

}