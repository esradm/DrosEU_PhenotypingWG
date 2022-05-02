
#############################################################
###################### SURVIVAL CURVES ######################
#############################################################


##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(MetBrewer)


##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")


##### load data
droseu <- readRDS("Data/droseu_master_list_2022-05-02.rds")
pops <- readRDS("InfoTables/DrosEU_Populations.rds")

##### source functions
source("Functions/survival_functions.R")


##### create output directory
surv_dir <- "SurvivalAnalyses"
dir.create(surv_dir, showWarnings = F) 


##### define theme and colors

droseu_theme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black",),  axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16),axis.title.y = element_text(size = 16))


myColors <- met.brewer("Johnson", 9)
names(myColors) <- as.factor(c("AK", "GI", "KA", "MA", "MU", "RE", "UM", "VA", "YE"))
colScale <- scale_colour_manual(name = "Population", values = myColors)





############# LIFESPAN #############

# create output directory
out_dir <- "Lifespan"
dir.create(file.path(surv_dir, out_dir), showWarnings = F) 

# get proportions
LS_prop_pop <- droseu$lsm %>%
  filter(Censor != 1) %>%
  rename(AgeAtDeath = LSM_AgeAtDeath_days) %>%
  group_split(Supervisor.PI, Population, Sex) %>%
  lapply(., survProp, cuts = seq(0, 120, 5)) %>%
  bind_rows()

# Flatt females
LS_F_Flatt <- ggplot(data = filter(LS_prop_pop, Supervisor.PI == "Flatt" & Sex == "F")) + 
  geom_line(aes(y = PropSurv, x = AgeAtDeath, colour = Population), size = 1.5) +
  scale_x_continuous(breaks=seq(0, 120, 20)) +
  labs(title = "LS_F_Flatt", x = "Age at death (days)", y = "Proportion survival") +
  colScale + droseu_theme 

pdf(file.path(surv_dir, out_dir, "p_LS_F_Flatt_survival_curves.pdf"), width=8, height=5)
LS_F_Flatt
dev.off()

# Flatt males
LS_M_Flatt <- ggplot(data = filter(LS_prop_pop, Supervisor.PI == "Flatt" & Sex == "M")) + 
  geom_line(aes(y = PropSurv, x = AgeAtDeath, colour = Population), size = 1.5) +
  scale_x_continuous(breaks=seq(0, 120, 20)) +
  labs(title = "LS_M_Flatt", x = "Age at death (days)", y = "Proportion survival") +
  colScale + droseu_theme 

pdf(file.path(surv_dir, out_dir, "p_LS_M_Flatt_survival_curves.pdf"), width=8, height=5)
LS_M_Flatt
dev.off()


# Parsch females
LS_F_Parsch <- ggplot(data = filter(LS_prop_pop, Supervisor.PI == "Parsch" & Sex == "F")) + 
  geom_line(aes(y = PropSurv, x = AgeAtDeath, colour = Population), size = 1.5) +
  scale_x_continuous(breaks=seq(0, 120, 20)) +
  labs(title = "LS_F_Parsch", x = "Age at death (days)", y = "Proportion survival") +
  colScale + droseu_theme 

pdf(file.path(surv_dir, out_dir, "p_LS_F_Parsch_survival_curves.pdf"), width=8, height=5)
LS_F_Parsch
dev.off()

# Parsch males
LS_M_Parsch <- ggplot(data = filter(LS_prop_pop, Supervisor.PI == "Parsch" & Sex == "M")) + 
  geom_line(aes(y = PropSurv, x = AgeAtDeath, colour = Population), size = 1.5) +
  scale_x_continuous(breaks=seq(0, 120, 20)) +
  labs(title = "LS_M_Parsch", x = "Age at death (days)", y = "Proportion survival") +
  colScale + droseu_theme 

pdf(file.path(surv_dir, out_dir, "p_LS_M_Parsch_survival_curves.pdf"), width=8, height=5)
LS_M_Parsch
dev.off()


# Pasyukova females
LS_F_Pasyukova <- ggplot(data = filter(LS_prop_pop, Supervisor.PI == "Pasyukova" & Sex == "F")) + 
  geom_line(aes(y = PropSurv, x = AgeAtDeath, colour = Population), size = 1.5) +
  scale_x_continuous(breaks=seq(0, 120, 20)) +
  labs(title = "LS_F_Pasyukova", x = "Age at death (days)", y = "Proportion survival") +
  colScale + droseu_theme 

pdf(file.path(surv_dir, out_dir, "p_LS_F_Pasyukova_survival_curves.pdf"), width=8, height=5)
LS_F_Pasyukova
dev.off()

# Pasyukova males
LS_M_Pasyukova <- ggplot(data = filter(LS_prop_pop, Supervisor.PI == "Pasyukova" & Sex == "M")) + 
  geom_line(aes(y = PropSurv, x = AgeAtDeath, colour = Population), size = 1.5) +
  scale_x_continuous(breaks=seq(0, 120, 20)) +
  labs(title = "LS_M_Pasyukova", x = "Age at death (days)", y = "Proportion survival") +
  colScale + droseu_theme 

pdf(file.path(surv_dir, out_dir, "p_LS_M_Pasyukova_survival_curves.pdf"), width=8, height=5)
LS_M_Pasyukova
dev.off()






############# CHILL-COMA #############

# create output directory
out_dir <- "ChillComa"
dir.create(file.path(surv_dir, out_dir), showWarnings = F) 

# get proportions
CCRT_prop_pop <- droseu$ccrt %>%
  filter(Censor != 1) %>%
  rename(AgeAtDeath = CCRT_seconds) %>%
  group_split(Supervisor.PI, Population, Sex) %>%
  lapply(., survProp, cuts = seq(0, 3600, 360)) %>%
  bind_rows()

# Vieira females
CCRT_F_Vieira <- ggplot(data = filter(CCRT_prop_pop, Supervisor.PI == "Vieira" & Sex == "F")) + 
  geom_line(aes(y = 1-PropSurv, x = AgeAtDeath, colour = Population), size = 1.5) +
  scale_x_continuous(breaks = seq(0, 3600, 600)) +
  labs(title = "CCRT_F_Vieira", x = "CCRT (sec)", y = "Proportion recovery") +
  colScale + droseu_theme 

pdf(file.path(surv_dir, out_dir, "p_CCRT_F_Vieira_survival_curves.pdf"), width=8, height=5)
CCRT_F_Vieira
dev.off()

# Vieira males
CCRT_M_Vieira <- ggplot(data = filter(CCRT_prop_pop, Supervisor.PI == "Vieira" & Sex == "M")) + 
  geom_line(aes(y = 1-PropSurv, x = AgeAtDeath, colour = Population), size = 1.5) +
  scale_x_continuous(breaks = seq(0, 3600, 600)) +
  labs(title = "CCRT_M_Vieira", x = "CCRT (sec)", y = "Proportion recovery") +
  colScale + droseu_theme 

pdf(file.path(surv_dir, out_dir, "p_CCRT_M_Vieira_survival_curves.pdf"), width=8, height=5)
CCRT_M_Vieira
dev.off()


# Mensch females
CCRT_F_Mensch <- ggplot(data = filter(CCRT_prop_pop, Supervisor.PI == "Mensch" & Sex == "F")) + 
  geom_line(aes(y = 1-PropSurv, x = AgeAtDeath, colour = Population), size = 1.5) +
  scale_x_continuous(breaks = seq(0, 3600, 600)) +
  labs(title = "CCRT_F_Mensch", x = "CCRT (sec)", y = "Proportion recovery") +
  colScale + droseu_theme 

pdf(file.path(surv_dir, out_dir, "p_CCRT_F_Mensch_survival_curves.pdf"), width=8, height=5)
CCRT_F_Mensch
dev.off()

# Mensch males
CCRT_M_Mensch <- ggplot(data = filter(CCRT_prop_pop, Supervisor.PI == "Mensch" & Sex == "M")) + 
  geom_line(aes(y = 1-PropSurv, x = AgeAtDeath, colour = Population), size = 1.5) +
  scale_x_continuous(breaks = seq(0, 3600, 600)) +
  labs(title = "CCRT_M_Mensch", x = "CCRT (sec)", y = "Proportion recovery") +
  colScale + droseu_theme 

pdf(file.path(surv_dir, out_dir, "p_CCRT_M_Mensch_survival_curves.pdf"), width=8, height=5)
CCRT_M_Mensch
dev.off()



############# STARVATION RESISTANCE #############



#ggplot(filter(droseu$sr, Supervisor.PI == "Gonzalez" & Sex == "M")) + 
#  geom_boxplot(aes(x = Population, y = AgeAtDeath_hours))


# create output directory
out_dir <- "Starvation"
dir.create(file.path(surv_dir, out_dir), showWarnings = F) 

# get proportions
SR_prop_pop <- droseu$sr %>%
  rename(AgeAtDeath = AgeAtDeath_hours) %>%
  group_split(Supervisor.PI, Population, Sex) %>%
  lapply(., survProp, cuts = seq(0, 240, 8)) %>%
  bind_rows()

# Gonzalez females
SR_F_Gonzalez <- ggplot(data = filter(SR_prop_pop, Supervisor.PI == "Gonzalez" & Sex == "F")) + 
  geom_line(aes(y = PropSurv, x = AgeAtDeath, colour = Population), size = 1.5) +
  scale_x_continuous(breaks = seq(0, 240, 24)) +
  labs(title = "SR_F_Gonzalez", x = "Age at death (hours)", y = "Proportion survival") +
  colScale + droseu_theme 

pdf(file.path(surv_dir, out_dir, "p_SR_F_Gonzalez_survival_curves.pdf"), width=8, height=5)
SR_F_Gonzalez
dev.off()

# Gonzalez males
SR_M_Gonzalez <- ggplot(data = filter(SR_prop_pop, Supervisor.PI == "Gonzalez" & Sex == "M")) + 
  geom_line(aes(y = PropSurv, x = AgeAtDeath, colour = Population), size = 1.5) +
  scale_x_continuous(breaks = seq(0, 240, 24)) +
  labs(title = "SR_M_Gonzalez", x = "Age at death (hours)", y = "Proportion survival") +
  colScale + droseu_theme 

pdf(file.path(surv_dir, out_dir, "p_SR_M_Gonzalez_survival_curves.pdf"), width=8, height=5)
SR_M_Gonzalez
dev.off()


# Onder females
SR_F_Onder <- ggplot(data = filter(SR_prop_pop, Supervisor.PI == "Onder" & Sex == "F")) + 
  geom_line(aes(y = PropSurv, x = AgeAtDeath, colour = Population), size = 1.5) +
  scale_x_continuous(breaks = seq(0, 240, 24)) +
  labs(title = "SR_F_Onder", x = "Age at death (hours)", y = "Proportion survival") +
  colScale + droseu_theme 

pdf(file.path(surv_dir, out_dir, "p_SR_F_Onder_survival_curves.pdf"), width=8, height=5)
SR_F_Onder
dev.off()

# Onder males
SR_M_Onder <- ggplot(data = filter(SR_prop_pop, Supervisor.PI == "Onder" & Sex == "M")) + 
  geom_line(aes(y = PropSurv, x = AgeAtDeath, colour = Population), size = 1.5) +
  scale_x_continuous(breaks = seq(0, 240, 24)) +
  labs(title = "SR_M_Onder", x = "Age at death (hours)", y = "Proportion survival") +
  colScale + droseu_theme 

pdf(file.path(surv_dir, out_dir, "p_SR_M_Onder_survival_curves.pdf"), width=8, height=5)
SR_M_Onder
dev.off()


# Pasyukova females
SR_F_Pasyukova <- ggplot(data = filter(SR_prop_pop, Supervisor.PI == "Pasyukova" & Sex == "F")) + 
  geom_line(aes(y = PropSurv, x = AgeAtDeath, colour = Population), size = 1.5) +
  scale_x_continuous(breaks = seq(0, 240, 24)) +
  labs(title = "SR_F_Pasyukova", x = "Age at death (hours)", y = "Proportion survival") +
  colScale + droseu_theme 

pdf(file.path(surv_dir, out_dir, "p_SR_F_Pasyukova_survival_curves.pdf"), width=8, height=5)
SR_F_Pasyukova
dev.off()

# Pasyukova males
SR_M_Pasyukova <- ggplot(data = filter(SR_prop_pop, Supervisor.PI == "Pasyukova" & Sex == "M")) + 
  geom_line(aes(y = PropSurv, x = AgeAtDeath, colour = Population), size = 1.5) +
  scale_x_continuous(breaks = seq(0, 240, 24)) +
  labs(title = "SR_M_Pasyukova", x = "Age at death (hours)", y = "Proportion survival") +
  colScale + droseu_theme 

pdf(file.path(surv_dir, out_dir, "p_SR_M_Pasyukova_survival_curves.pdf"), width=8, height=5)
SR_M_Pasyukova
dev.off()



############# HEAT-SHOCK MORTALITY #############


# create output directory
out_dir <- "HeatShock"
dir.create(file.path(surv_dir, out_dir), showWarnings = F) 

# get proportions
HSM_prop_pop <- droseu$hsm %>%
  filter(Censor != 1) %>%
  rename(AgeAtDeath = TimeDeath_min) %>%
  group_split(Supervisor.PI, Population, Sex) %>%
  lapply(., survProp, cuts = seq(0, 600, 60)) %>%
  bind_rows()

# Parsch females
HSM_F_Parsch <- ggplot(data = filter(HSM_prop_pop, Supervisor.PI == "Parsch" & Sex == "F")) + 
  geom_line(aes(y = PropSurv, x = AgeAtDeath, colour = Population), size = 1.5) +
  scale_x_continuous(breaks = seq(0, 600, 120)) +
  labs(title = "HSM_F_Parsch", x = "Age at death (min)", y = "Proportion survival") +
  colScale + droseu_theme 

pdf(file.path(surv_dir, out_dir, "p_HSM_F_Parsch_survival_curves.pdf"), width=8, height=5)
HSM_F_Parsch
dev.off()

# Parsch males
HSM_M_Parsch <- ggplot(data = filter(HSM_prop_pop, Supervisor.PI == "Parsch" & Sex == "M")) + 
  geom_line(aes(y = PropSurv, x = AgeAtDeath, colour = Population), size = 1.5) +
  scale_x_continuous(breaks = seq(0, 600, 120)) +
  labs(title = "HSM_M_Parsch", x = "Age at death (min)", y = "Proportion survival") +
  colScale + droseu_theme 

pdf(file.path(surv_dir, out_dir, "p_HSM_M_Parsch_survival_curves.pdf"), width=8, height=5)
HSM_M_Parsch
dev.off()


# Vieira females
HSM_F_Vieira <- ggplot(data = filter(HSM_prop_pop, Supervisor.PI == "Vieira" & Sex == "F")) + 
  geom_line(aes(y = PropSurv, x = AgeAtDeath, colour = Population), size = 1.5) +
  scale_x_continuous(breaks = seq(0, 600, 120)) +
  labs(title = "HSM_F_Vieira", x = "Age at death (min)", y = "Proportion survival") +
  colScale + droseu_theme 

pdf(file.path(surv_dir, out_dir, "p_HSM_F_Vieira_survival_curves.pdf"), width=8, height=5)
HSM_F_Vieira
dev.off()

# Vieira males
HSM_M_Vieira <- ggplot(data = filter(HSM_prop_pop, Supervisor.PI == "Vieira" & Sex == "M")) + 
  geom_line(aes(y = PropSurv, x = AgeAtDeath, colour = Population), size = 1.5) +
  scale_x_continuous(breaks = seq(0, 600, 120)) +
  labs(title = "HSM_M_Vieira", x = "Age at death (min)", y = "Proportion survival") +
  colScale + droseu_theme 

pdf(file.path(surv_dir, out_dir, "p_HSM_M_Vieira_survival_curves.pdf"), width=8, height=5)
HSM_M_Vieira
dev.off()




