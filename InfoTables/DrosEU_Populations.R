
##########################################################################
##################### DROSEU PHENOTYPING POPULATIONS #####################
##########################################################################

rm(list = ls())

##### set wd

setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")


##### load packages and source functions

library(tidyverse)

##### load data

updated_coords <- read.csv("InfoTables/DrosEU_Coordinates.csv")


##### populations

pops <- updated_coords %>%
  mutate(
    Population = str_to_upper(str_sub(Location, 1, 2))
  ) %>%
  relocate(Population, Country, Country_code)


ids_by_lat <- arrange(pops, Latitude)$Population
ids_by_lon <- arrange(pops, Longitude)$Population
ids_by_alt <- arrange(pops, Altitude)$Population


pops_by_lat <- mutate(pops, Population = factor(Population, levels = ids_by_lat)) %>%
  arrange(Population) %>%
  mutate(Color = c("#cab2d6", "#ff7f00", "#fdbf6f", "#e31a1c", "#fb9a99",
                  "#33a02c", "#b2df8a", "#1f78b4", "#a6cee3"))

pops_by_lon <- mutate(pops_by_lat, Population = factor(Population, levels = ids_by_lon)) %>%
  arrange(Population)

pops_by_alt <- mutate(pops_by_lat, Population = factor(Population, levels = ids_by_alt)) %>%
  arrange(Population)

pops <- list(by_lat = pops_by_lat, by_lon = pops_by_lon, by_alt = pops_by_alt)
saveRDS(pops, file = "InfoTables/Droseu_Populations.rds")