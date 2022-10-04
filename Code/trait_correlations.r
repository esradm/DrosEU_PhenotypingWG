


library(tidyverse)
library(GGally)

setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")

dir.create("TraitCorrelations/", showWarnings = FALSE)


rc <- read.csv("./MetaAnalyses/all_models_line_meta_compound_random_coefs.csv")
rc_diet <- read.csv("./LinearModelsPop/all_models_line_random_coefs_similar_diet.csv")


# nine traits in common between males and females
c9 <- c("CCRT", "CSM", "HSM", "DT_A", "DW", "LS", "SR", "TL", "WA_L")
fmax <- c(c9, "Dia", "Fec", "Pgm_Total") # female traits
fmaxp <- c(fmax, "Via") # female traits plus


rc_wide <- rc %>%
    dplyr::select(Trait, Population, Line, Sex, Estimate) %>%
    pivot_wider(values_from = Estimate, names_from = c(Trait, Sex))

cor_f9 <- ggpairs(select(rc_wide, matches(c9) & matches("_F")))
ggsave(cor_f9,
    file = "TraitCorrelations/trait_correlations_f9_line_random_coefs.pdf",
    width = 10, height = 10)

cor_m9 <- ggpairs(select(rc_wide, matches(c9) & matches("_M")))
ggsave(cor_m9,
    file = "TraitCorrelations/trait_correlations_m9_line_random_coefs.pdf",
    width = 10, height = 10)

cor_fmax <- ggpairs(select(rc_wide, matches(fmax) & matches("_F")))
ggsave(cor_fmax,
    file = "TraitCorrelations/trait_correlations_fmax_line_random_coefs.pdf",
    width = 10, height = 10)

cor_fmaxp <- ggpairs(select(rc_wide, matches(fmaxp) & !matches("_M")))
ggsave(cor_fmaxp,
    file = "TraitCorrelations/trait_correlations_fmaxp_line_random_coefs.pdf",
    width = 10, height = 10)



rc_diet_wide <- rc_diet %>%
    dplyr::select(Trait, Population, Line, Sex, Coef) %>%
    pivot_wider(values_from = Coef, names_from = c(Trait, Sex))

cor_f9_diet <- ggpairs(select(rc_diet_wide, matches(c9) & matches("_F")))
ggsave(cor_f9_diet,
    file = "TraitCorrelations/trait_correlations_f9_line_random_coefs_diet.pdf",
    width = 10, height = 10)

cor_m9_diet <- ggpairs(select(rc_diet_wide, matches(c9) & matches("_M")))
ggsave(cor_m9_diet,
    file = "TraitCorrelations/trait_correlations_m9_line_random_coefs_diet.pdf",
    width = 10, height = 10)

cor_fmax_diet <- ggpairs(select(rc_diet_wide, matches(fmax) & matches("_F")))
ggsave(cor_fmax_diet,
    file = "TraitCorrelations/trait_correlations_fmax_line_random_coefs_diet.pdf",
    width = 10, height = 10)

cor_fmaxp_diet <- ggpairs(select(rc_diet_wide, matches(fmaxp) & !matches("_M")))
ggsave(cor_fmaxp_diet,
    file = "TraitCorrelations/trait_correlations_fmaxp_line_random_coefs_diet.pdf",
    width = 10, height = 10)