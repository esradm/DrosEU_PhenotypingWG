

rm(list = ls())

##### set wd
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")

##### load packages
library(tidyverse)

##### laod data and fucntions
droseu <- readRDS("Data/droseu_master_list_2022-05-02.rds")
pops <- readRDS("InfoTables/DrosEU_Populations.rds")
traits <- read.csv("InfoTables/trait_names.csv")

##### create output directory structure
dir.create("SummaryTables", showWarnings = FALSE)


##### merge with population info

common_var <- c("Lab", "Diet", "Batch", "Population", "Line", "Sex", "Replicate", "Individual")
trait_var <- traits$Trait_name_raw


##### add extra columns for easier batch processing and split data by single trait variables

droseu_long <- droseu
for (i in seq_along(droseu_long)){
  d <- droseu_long[[i]]
  colnames(d)[colnames(d) == "Supervisor.PI"] <- "Lab"
  colnames(d)[colnames(d) == "ReplicateVial"] <- "Replicate"
  colnames(d)[colnames(d) == "ReplicateCage"] <- "Replicate"
  colnames(d)[colnames(d) == "ReplicateChamber"] <- "Replicate"
  if (!"Line" %in% colnames(d)) d$Line <- NA
  if (!"Batch" %in% colnames(d)) d$Batch <- NA
  if (!"Sex" %in% colnames(d)) d$Sex <- NA
  if (!"ReplicateVial" %in% colnames(d)) d$Replicate <- NA
  if (!"Individual" %in% colnames(d)) d$Individual <- NA
  d <- d[, colnames(d) %in% c(common_var, trait_var)]
  d <- pivot_longer(d, cols = -all_of(common_var), names_to = "Trait_name", values_to = "Value")
  d$Trait <- names(droseu_long)[i]
  droseu_long[[i]] <- d
}

droseu_long <- bind_rows(droseu_long) %>%
  rename(Trait_handle = Trait)


trait_names <- read.csv("InfoTables/trait_names.csv") %>%
  dplyr::select(Trait_handle, Trait, Trait_name_raw, Legend) %>%
  rename(Trait_name = Trait_name_raw)

droseu_long <- inner_join(droseu_long, trait_names)



### summarise for each trait

droseu_trait <- group_by(droseu_long, Trait_name) %>%
  group_split()


droseu_summary <- list()

for (i in seq_along(droseu_trait)) {

  x <- droseu_trait[[i]]

  n_obs <- x %>%
    group_by(Lab, Diet, Sex, Population, Line, ReplicateVial) %>%
    summarise(
      Observation_count = n(),
      Batch_count = length(unique(Batch)),
      Batch_ids = paste(sort(unique(Batch)), collapse = "; ")
    )

  n_reps <- x %>%
    group_by(Lab, Diet, Sex, Population, Line) %>%
    summarise(
      Replicate_count = length(unique(ReplicateVial)),
      Replicate_ids = paste(sort(unique(ReplicateVial)), collapse = "; "),
      Observation_count = n(),
      Batch_count = length(unique(Batch)),
      Batch_ids = paste(sort(unique(Batch)), collapse = "; ")
    )

  n_lines <- x %>%
    group_by(Lab, Diet, Sex, Population) %>%
    summarise(
      Line_count = length(unique(Line)),
      Line_ids = paste(sort(unique(Line)), collapse = "; "),
      Replicate_count = length(unique(ReplicateVial)),
      Replicate_ids = paste(sort(unique(ReplicateVial)), collapse = "; "),
      Observation_count = n(),
      Batch_count = length(unique(Batch)),
      Batch_ids = paste(sort(unique(Batch)), collapse = "; ")
    )

  n_pops <- x %>%
    group_by(Lab, Diet, Sex) %>%
    summarise(
      Population_count = length(unique(Population)),
      Population_ids = paste(sort(unique(Population)), collapse = "; "),
      Line_count = length(unique(Line)),
      Line_ids = paste(sort(unique(Line)), collapse = "; "),
      Replicate_count = length(unique(ReplicateVial)),
      Replicate_ids = paste(sort(unique(ReplicateVial)), collapse = "; "),
      Observation_count = n(),
      Batch_count = length(unique(Batch)),
      Batch_ids = paste(sort(unique(Batch)), collapse = "; ")
    )

  n_sex <- x %>%
    group_by(Lab, Diet) %>%
    summarise(
      Sex_count = length(unique(Sex)),
      Sex_ids = paste(sort(unique(Sex)), collapse = "; "),
      Population_count = length(unique(Population)),
      Population_ids = paste(sort(unique(Population)), collapse = "; "),
      Line_count = length(unique(Line)),
      Line_ids = paste(sort(unique(Line)), collapse = "; "),
      Replicate_count = length(unique(ReplicateVial)),
      Replicate_ids = paste(sort(unique(ReplicateVial)), collapse = "; "),
      Observation_count = n(),
      Batch_count = length(unique(Batch)),
      Batch_ids = paste(sort(unique(Batch)), collapse = "; ")
    )

  by_batch <- x %>%
    group_by(Lab, Diet, Sex, Batch) %>%
    summarise(
      Population_count = length(unique(Population)),
      Population_ids = paste(sort(unique(Population)), collapse = "; "),
      Line_count = length(unique(Line)),
      Line_ids = paste(sort(unique(Line)), collapse = "; "),
      Replicate_count = length(unique(ReplicateVial)),
      Replicate_ids = paste(sort(unique(ReplicateVial)), collapse = "; "),
      Observation_count = n()
    )

  droseu_summary[[i]] <- list(
    obs = n_obs,
    rep = n_reps,
    line = n_lines,
    pop = n_pops,
    sex = n_sex,
    batch = by_batch
  )

}




trait_out <- lapply(droseu_trait, function(x) {
  paste(unique(x$Trait), unique(x$Trait_name), sep = "_")
})


for (i in seq_along(droseu_summary)) {
  file_path <- paste("SummaryTables", trait_out[[i]], sep = "/")
  trait <- droseu_summary[[i]]
  for (j in seq_along(trait)) {
    file_name <- paste(file_path, "_", names(trait)[j], ".csv", sep = "")
    write.csv(trait[[j]], file_name, row.names = FALSE)
  }
}