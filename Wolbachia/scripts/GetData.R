
##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(lme4)
#
##### set working directory
setwd("/Users/martinkapun/Documents/GitHub/DrosEU_PhenotypingWG")

##### source functions
source("Code/functions.R")

##### load data
droseu <- readRDS("Data/droseu_master_list_2022-05-02.rds")

## read Wolbachia dataset
Wolb <- read.table("Data/Wolbachia.txt",
  header = TRUE,
  na.strings = "NA"
)

Diet <- read.table("Wolbachia/data/DrosEU_Diets_Feb22.csv",
  header = TRUE,
  na.strings = "NA",
  sep = ","
)

Diet <- Diet %>%
  group_by(PI.Lab.head) %>%
  summarize(P.C = mean(P.C))

colnames(Diet) <- c("Supervisor.PI", "PC.ratio")

Wolb$Line <- toupper(Wolb$Line)
Wolb.cons <- na.omit(Wolb)

out_dir <- "Wolbachia/data"

############# Viability #############

# create output directory
dir.create(file.path(out_dir, "Viability"), showWarnings = FALSE)

Viability <- droseu$via %>%
  inner_join(Wolb.cons, by = c("Line"))

Viability <- Viability %>%
  inner_join(Diet, by = c("Supervisor.PI"))

write.table(Viability,
  paste(out_dir, "Viability", "Viability.txt", sep = "/"),
  quote = FALSE,
  row.names = FALSE, sep = "\t"
)

############# Fecundity #############

# create output directory
dir.create(file.path(out_dir, "Fecundity"), showWarnings = FALSE)

Fecundity <- droseu$fec %>%
  inner_join(Wolb.cons, by = c("Line"))

Fecundity <- Fecundity %>%
  inner_join(Diet, by = c("Supervisor.PI"))


write.table(Fecundity,
  paste(out_dir, "Fecundity", "Fecundity.txt", sep = "/"),
  quote = FALSE,
  row.names = FALSE, sep = "\t"
)

############# DevelopmentTime_ETP #############

# create output directory
out_dir <- "Wolbachia/data"
dir.create(file.path(out_dir, "DevelopmentTime_ETP"), showWarnings = FALSE)

DevelopmentTime_ETP <- droseu$dtp %>%
  inner_join(Wolb.cons, by = c("Line"))

DevelopmentTime_ETP <- DevelopmentTime_ETP %>%
  inner_join(Diet, by = c("Supervisor.PI"))

write.table(DevelopmentTime_ETP,
  paste(out_dir, "DevelopmentTime_ETP", "DevelopmentTime_ETP.txt", sep = "/"),
  quote = FALSE,
  row.names = FALSE, sep = "\t"
)

############# DevelopmentTime_ETA #############

# create output directory
out_dir <- "Wolbachia/data"
dir.create(file.path(out_dir, "DevelopmentTime_ETA"), showWarnings = FALSE)

DevelopmentTime_ETA <- droseu$dta %>%
  inner_join(Wolb.cons, by = c("Line"))

DevelopmentTime_ETA <- DevelopmentTime_ETA %>%
  inner_join(Diet, by = c("Supervisor.PI"))

write.table(DevelopmentTime_ETA,
  paste(out_dir, "DevelopmentTime_ETA", "DevelopmentTime_ETA.txt", sep = "/"),
  quote = FALSE,
  row.names = FALSE, sep = "\t"
)

############# DryWeight #############

# create output directory
out_dir <- "Wolbachia/data"
dir.create(file.path(out_dir, "DryWeight"), showWarnings = FALSE)

DryWeight <- droseu$dw %>%
  inner_join(Wolb.cons, by = c("Line"))

DryWeight <- DryWeight %>%
  inner_join(Diet, by = c("Supervisor.PI"))


write.table(DryWeight,
  paste(out_dir, "DryWeight", "DryWeight.txt", sep = "/"),
  quote = FALSE,
  row.names = FALSE, sep = "\t"
)

############# ThoraxLength #############

# create output directory
out_dir <- "Wolbachia/data"
dir.create(file.path(out_dir, "ThoraxLength"), showWarnings = FALSE)

ThoraxLength <- droseu$tl %>%
  inner_join(Wolb.cons, by = c("Line"))

ThoraxLength <- ThoraxLength %>%
  inner_join(Diet, by = c("Supervisor.PI"))

write.table(ThoraxLength,
  paste(out_dir, "ThoraxLength", "ThoraxLength.txt", sep = "/"),
  quote = FALSE,
  row.names = FALSE, sep = "\t"
)

############# WingArea #############

# create output directory
out_dir <- "Wolbachia/data"
dir.create(file.path(out_dir, "WingArea"), showWarnings = FALSE)

WingArea <- droseu$wa %>%
  inner_join(Wolb.cons, by = c("Line"))

WingArea <- WingArea %>%
  inner_join(Diet, by = c("Supervisor.PI"))

write.table(WingArea,
  paste(out_dir, "WingArea", "WingArea.txt", sep = "/"),
  quote = FALSE,
  row.names = FALSE, sep = "\t"
)

############# Lifespan #############

# create output directory
out_dir <- "Wolbachia/data"
dir.create(file.path(out_dir, "Lifespan"), showWarnings = FALSE)

Lifespan <- droseu$lsl %>%
  inner_join(Wolb.cons, by = c("Line"))

Lifespan <- Lifespan %>%
  inner_join(Diet, by = c("Supervisor.PI"))

write.table(Lifespan,
  paste(out_dir, "Lifespan", "Lifespan.txt", sep = "/"),
  quote = FALSE,
  row.names = FALSE, sep = "\t"
)

############# ColdShockMortality #############

# create output directory
out_dir <- "Wolbachia/data"
dir.create(file.path(out_dir, "ColdShockMortality"), showWarnings = FALSE)

ColdShockMortality <- droseu$csm %>%
  inner_join(Wolb.cons, by = c("Line"))

ColdShockMortality <- ColdShockMortality %>%
  inner_join(Diet, by = c("Supervisor.PI"))

write.table(ColdShockMortality,
  paste(out_dir, "ColdShockMortality", "ColdShockMortality.txt", sep = "/"),
  quote = FALSE,
  row.names = FALSE, sep = "\t"
)

############# ChillComa #############

# create output directory
out_dir <- "Wolbachia/data"
dir.create(file.path(out_dir, "ChillComa"), showWarnings = FALSE)

ChillComa <- droseu$ccrt %>%
  inner_join(Wolb.cons, by = c("Line"))

ChillComa <- ChillComa %>%
  inner_join(Diet, by = c("Supervisor.PI"))

write.table(ChillComa,
  paste(out_dir, "ChillComa", "ChillComa.txt", sep = "/"),
  quote = FALSE,
  row.names = FALSE, sep = "\t"
)

############# HeatShock #############

# create output directory
out_dir <- "Wolbachia/data"
dir.create(file.path(out_dir, "HeatShock"), showWarnings = FALSE)

HeatShock <- droseu$hsm %>%
  inner_join(Wolb.cons, by = c("Line"))

HeatShock <- HeatShock %>%
  inner_join(Diet, by = c("Supervisor.PI"))

write.table(HeatShock,
  paste(out_dir, "HeatShock", "HeatShock.txt", sep = "/"),
  quote = FALSE,
  row.names = FALSE, sep = "\t"
)

############# Diapause #############

# create output directory
out_dir <- "Wolbachia/data"
dir.create(file.path(out_dir, "Diapause"), showWarnings = FALSE)

Diapause <- droseu$dia %>%
  inner_join(Wolb.cons, by = c("Line"))

Diapause <- Diapause %>% # nolint # nolint
  inner_join(Diet, by = c("Supervisor.PI"))

write.table(Diapause,
  paste(out_dir, "Diapause", "Diapause.txt", sep = "/"),
  quote = FALSE,
  row.names = FALSE, sep = "\t"
)

############# Starvation #############

# create output directory
out_dir <- "Wolbachia/data"
dir.create(file.path(out_dir, "Starvation"), showWarnings = FALSE)

Starvation <- droseu$sr %>%
  inner_join(Wolb.cons, by = c("Line"))

Starvation <- Starvation %>%
  inner_join(Diet, by = c("Supervisor.PI"))

write.table(Starvation,
  paste(out_dir, "Starvation", "Starvation.txt", sep = "/"),
  quote = FALSE,
  row.names = FALSE, sep = "\t"
)

############# Pigmentation #############

# create output directory
out_dir <- "Wolbachia/data"
dir.create(file.path(out_dir, "Pigmentation"), showWarnings = FALSE)

Pigmentation <- droseu$pgm %>%
  inner_join(Wolb.cons, by = c("Line"))

Pigmentation <- Pigmentation %>%
  inner_join(Diet, by = c("Supervisor.PI"))

write.table(Pigmentation,
  paste(out_dir, "Pigmentation", "Pigmentation.txt", sep = "/"),
  quote = FALSE,
  row.names = FALSE, sep = "\t"
)
