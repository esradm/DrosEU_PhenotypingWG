


############################################################
###################### HERITABILITIES ######################
############################################################



###### some things to consider




# for LSP and CETS, H2 can be calculated for Population

# what should be the input data when measures have been made at the individual level? Is each individual a line replicate or should we use the average value per replicate vial?

# for between labs H2, how do we pool teh data from all labs together or do we first get an average line value per lab? Then each lab is treated as a replicate measure?

# for between labs H2, instead of using mean Line values, one could consider using Line random coefficients extracted from individual labs linear mixed models

# maybe just stick to within labs H2 in the first place?


##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(lme4)

##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")



################# DEFINE H2 FUNCTION ################# 

# this function is based on a linear mixed model (lmer) in which genotypes (i.e. lines) are included as random-effect variables. It allows for the extraction of genetic and environmental variances used to compute broad sense heritability, H2.
# usage: H2_lmm(d = via_df, phenot = "via", genot = "Line")
# d: data.frame with columns for trait values and genotypes
# phenot: trait values
# genot: genotype ids
# min_rep_genot: minimal genotype replication allowed

H2_lmm <- function(d, phenot, genot, min_rep_genot = 1) {
  require("lme4")
  # check genotypes
  n_unique_genot <- length(unlist(unique(d[,genot])))
  d <- d %>% group_by(.data[[genot]]) %>% filter(n() >= min_rep_genot)
  n_rep_genot <- length(unlist(unique(d[,genot])))
  # make sure that H can be computed by checking that genotypes have been measured more than once (i.e. that genotypes are replicated). If not the case H cannot be computed and a data.frame with NAs is returned.
  if (n_rep_genot < 2 | n_unique_genot == nrow(d)) { 
    data.frame(H2 = NA, Vp = NA, Vg = NA, Ve = NA, Genot = genot, N_unique_genot = n_unique_genot, Min_rep_genot = min_rep_genot, N_rep_genot = n_rep_genot, Warnings = "Not enough genotype replicates")
  } else {
    # build and run the lmer model with genot as random-effect variable
    lmm <- lmer(as.formula(paste0(phenot, "~1 + 1|", genot)), data = d)	
    # extract variance components
    lmm.varcor <- as.data.frame(summary(lmm)$varcor)
    # variance explained by genot
    Vg <- lmm.varcor[lmm.varcor$grp == genot, "vcov"]
    # variance explained by environment (residuals)
    Ve <- lmm.varcor[lmm.varcor$grp == "Residual", "vcov"]
    # total variance
    Vp <- Vg + Ve
    # calculate H2
    H2 <- round(Vg/Vp, 2)
    # output
    w <- summary(lmm)$optinfo$conv$lme4$messages
    data.frame(H2, Vp, Vg, Ve, Genot = genot, N_unique_genot = n_unique_genot, Min_rep_genot = min_rep_genot, N_rep_genot = n_rep_genot, Warnings = ifelse(length(w) != 0, w, "NA"))
  }
}




################# LINE H2 WITHIN LAB WITH RAW DATA #################

# code below allows for automated calculation of H2 for all traits

# get data in
trait_names <- read.csv("InfoTables/trait_names.csv")
droseu <- readRDS("Data/droseu_master_list_2022-05-02.rds")

# add dummy Sex and Line columns if not present for easier batch processing
for (i in 1:length(droseu)) {
  if (!"Sex" %in% colnames(droseu[[i]])) droseu[[i]]$Sex <- "NA"
  if (!"Line" %in% colnames(droseu[[i]])) droseu[[i]]$Line <- "NA"
}

# keep variables of interest and make long table
for (i in 1:length(droseu)) {
  droseu[[i]] <- dplyr::select(droseu[[i]], contains(c("Supervisor.PI", "Population", "Line", "Sex", trait_names$Trait_name_raw)))
  droseu[[i]] <- dplyr::select(droseu[[i]], -contains("Population_"))
  droseu[[i]] <- pivot_longer(droseu[[i]], cols = contains(trait_names$Trait_name_raw), names_to = "Trait_name_raw", values_to = "Value")
}

# update trait names
droseu <- bind_rows(droseu) %>% 
  inner_join(trait_names, by = "Trait_name_raw") %>% 
  dplyr::select(-c("Trait_handle", "Trait_name_raw")) %>%
  rename(Lab = Supervisor.PI)

# keep traits that can be used for Line H2
droseu_line <- filter(droseu, !Trait %in% c("LSM", "LSP", "CETS"))

# split to list
droseu_line_list <- group_split(droseu_line, Trait, Sex, Lab)

# loop over traits to calculate H2
h2_line_within <- list()
for (i in 1:length(droseu_line_list)) {
  trait <- droseu_line_list[[i]][1, c("Trait", "Sex", "Lab")]
  h2_lmm <- H2_lmm(droseu_line_list[[i]] %>% filter(!Value %in% c(-Inf, Inf)), "Value", "Line", 1)
  h2_line_within[[i]] <- bind_cols(trait, h2_lmm) 
}
h2_line_within <- bind_rows(h2_line_within)
  
# save results
write.csv(h2_line_within, "Heritability/H2_line_raw_data_within_labs.csv", row.names = F)




################# LINE H2 BETWEEN LABS WITH RAW DATA #################


# split to list
droseu_line_list2 <- group_split(droseu_line, Trait, Sex)

# loop over traits to calculate H2
h2_line_between <- list()
for (i in 1:length(droseu_line_list2)) {
  trait <- droseu_line_list2[[i]][1, c("Trait", "Sex")]
  h2_lmm <- H2_lmm(droseu_line_list2[[i]] %>% filter(!Value %in% c(-Inf, Inf)), "Value", "Line", 1)
  h2_line_between[[i]] <- bind_cols(trait, h2_lmm) 
}
h2_line_between <- bind_rows(h2_line_between)

# save results
write.csv(h2_line_between, "Heritability/H2_line_raw_data_between_labs.csv", row.names = F)













################# LINE H2 BETWEEN LABS WITH LMERS RANDOM COEFS #################

# get data in
lrc <- read.csv("LinearModelsPop/all_models_line_random_coefs.csv")

# split to list
lrc_list <- group_split(lrc, Trait, Sex)

# calculate H2
trait <- bind_rows(lapply(lrc_list, function(x) x[1,c("Trait", "Sex")]))
h2_lmm <- bind_rows(lapply(lrc_list, H2_lmm, "Coef", "Line", 1))
h2_line_between_rc <- bind_cols(trait, h2_lmm)

# save results
write.csv(h2_line_between_rc, "Heritability/H2_line_random_coefs_between_labs.csv", row.names = F)



















