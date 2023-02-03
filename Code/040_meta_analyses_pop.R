
######################################################################
###################### POPULATION META ANALYSES ######################
######################################################################



# CCRT_F_meta.rds contains the output of subgroup meta
# CCRT_F_meta_summary.txt summary of subgroup meta as txt file 
# CCRT_F_meta_pop_compound_estimates.rds compound estimates
# CCRT_F_meta_pop_compound_estimates.txt compound estimates as txt file



##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(meta)
library(ggpubr)
library(MetBrewer)
library(ggrepel)

##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")

##### source functions
source("Code/functions.R")

##### set and create main output directory
meta_dir <- "MetaAnalyses"
dir.create(meta_dir, showWarnings = F) 

##### set input directory containing the linear model outputs
lmer_dir <- "LinearModelsPop"
lmer_pat <- "lmers_pop_model_estimates.rds"



############# PARTIAL DATA ############# 

##### define data that is considered partial and that should be removed from meta analyses

# Posnien in WA and TL 3 lines per pop
# Ritchie in WA and TL 5 lines per pop, we lose TL Males if removed
# Grath in Via and DT_A 10 lines for 3 pops

# would remove Posnien for sure, Ritchie and Grath to be discussed

#partial_labs <- c("Posnien", "Grath", "Ritchie")
partial_labs <- "Posnien"



############# SUBGROUP META ANALYSES ############# 

##### get all the line estimates, from all the traits
pops_coefs <- list.files(path = lmer_dir, recursive = T, full.names = T, pattern = lmer_pat)
print(pops_coefs)

##### remove Dia lmers, keep glmers instead
pops_coefs <- pops_coefs[-grep("Dia_lmers", pops_coefs)]



##### meta loop to get compound estimates per trait

for (f in 1:length(pops_coefs)) {
  # get linear models path
  fpath <- pops_coefs[f]
  # create output directory
  trait_dir <- str_split(fpath, "/", simplify = T)[2]
  trait_dir <- paste(meta_dir, trait_dir, sep = "/")
  dir.create(trait_dir, showWarnings = F) 
  # read linear models output
  m <- readRDS(fpath)
  # remove previoulsy defined partial data
  m_complete <- filter(m, !Lab %in% partial_labs)
  # make effects
  m_effects <- makeEffects(m_complete)
  # get the subtraits
  traits <- unique(m_effects$Trait)
  # loop over subtraits
  for(tr in 1:length(traits)) {
    trait <- filter(m_effects, Trait == traits[tr])
    # get the sexes
    sexes <- unique(trait$Sex)
    # loop over sexes and perform meta
    for (s in 1:length(sexes)) {
      # define output files names
      mod <- ifelse(grepl( "glmers", fpath), "glmers", "lmers") 
      out_sex_rds <- file.path(trait_dir, paste(traits[tr], sexes[s], mod, "pop_meta.rds", sep ="_")) 
      out_sex_txt <- sub(".rds", "_summary.txt", out_sex_rds)
      out_sex_comp_rds <- sub(".rds", "_compound_estimates.rds", out_sex_rds)
      out_sex_comp_txt <- sub(".rds", ".txt", out_sex_comp_rds)
      # meta analysis
      msex <- metagen(data = filter(trait, Sex == sexes[s]), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE, method.tau = "REML")
      # subgroup analysis
      mres <- update.meta(msex, subgroup = Population, tau.common = FALSE)
      # save meta object and get meta summary for each sex
      saveRDS(mres, file = out_sex_rds)
      capture.output(summary(mres), file = out_sex_txt) 
      # extract compound estimates from the meta object to a data.frame
      comp_estimates <- data.frame(Models = mod, 
                                   Trait = traits[tr], 
                                   Sex = sexes[s],
                                   Population = str_sub(mres$bylevs, 1, 2), 
                                   Estimate = mres$TE.random.w, 
                                   SE = mres$seTE.random.w, 
                                   LLEst = mres$lower.random.w, 
                                   ULEst = mres$upper.random.w,
                                   Q = mres$Q.b.random, 
                                   P = mres$pval.Q.b.random,
                                   N_lab = mres$k.w,
                                   N_lab_av = mean(mres$k.w))
      # save compound estimates object as rds and txt
      saveRDS(comp_estimates, file = out_sex_comp_rds)
      write.table(comp_estimates, out_sex_comp_txt, row.names = F, quote = F, sep = "\t")
    }
  }
}


############# COMBINE META ANALYSES OUTPUTS ############# 

##### as a global list
all_metas_pop_list <- rdsBatchReaderToList(path = meta_dir, recursive = T, full.names = T, pattern = "lmers_pop_meta_compound_estimates.rds")
saveRDS(all_metas_pop_list, file.path(file = meta_dir, "all_models_pop_meta_compound_estimates_list.rds"))

##### as collapsed list
all_metas_pop <- bind_rows(all_metas_pop_list)
saveRDS(all_metas_pop, file = file.path(meta_dir, "all_models_pop_meta_compound_estimates.rds"))
write.csv(all_metas_pop, file = file.path(meta_dir, "all_models_pop_meta_compound_estimates.csv"), row.names = F)




############# PVALUES AND MULTIPLE TESTING CORRECTION ############# 

##### read in meta results
all_metas_pop_list <- rdsBatchReaderToList(path = meta_dir, recursive = T, full.names = T, pattern = "lmers_pop_meta_compound_estimates.rds")

##### define function to extract p values and other statistcics
getPvalue <- function(x) {
    x %>% mutate(Min_lab = min(N_lab), Max_lab = max(N_lab)) %>% 
    select(Models, Trait, Sex, Q, P, Min_lab, Max_lab) %>%
    distinct()
}

##### apply the function to all meta results
all_metas_pvalues <- lapply(all_metas_pop_list, getPvalue) %>%
  bind_rows()


##### remove non informative meta analyses for DT_P, LA and TL males

all_metas_pvalues <- all_metas_pvalues %>% filter(Min_lab > 1)

##### correct p values for multiple testing

all_metas_pvalues_adj <- all_metas_pvalues %>%
  mutate(P_bonf = p.adjust(P, "bonferroni"),
         P_bh = p.adjust(P, "BH"))


##### define number of traits and apply Bonferroni correction
#n_traits <- 13 # unique traits we run meta on
#n_traits <- 22 # unique traits we run meta on + sex
#n_traits <- 17 - 1 # unique traits + subtraits we run meta on
#n_traits <- 27 - 1 # unique traits + subtraits we run meta on + sex
# TL_M is also not meta analysed, hence - 1


##### output pvalues
saveRDS(all_metas_pvalues_adj, file.path(meta_dir, "all_models_pop_meta_pvalues.rds"))
write.csv(all_metas_pvalues_adj, file.path(meta_dir, "all_models_pop_meta_pvalues.csv"), row.names = F)


############# Q AND P VALUES PLOT ############# 

bh_thresh <- sum(sort(all_metas_pvalues_adj$P) < 0.05) / nrow(all_metas_pvalues_adj) * 0.05
bonf_thresh <- 0.05 / nrow(all_metas_pvalues_adj)
  
  
pvalue_plot <- all_metas_pvalues_adj %>%
  ggplot(aes(x = Q, y = -log10(P))) +
  geom_point(size = 4, alpha = 0.3, pch = 19) +
  theme_bw(16) +
  geom_hline(yintercept = -log10(bonf_thresh), linetype = 2, size = 0.5, col = "red") +
  geom_hline(yintercept = -log10(bh_thresh), linetype = 2, size = 0.5) +
  labs(title = "Meta analyses Q and P values", x = "Q value", y = "-log10(Pvalue)") +
  geom_label_repel(data = mutate(all_metas_pvalues_adj, Label = ifelse(P_bh >= 0.05, NA, paste(Trait, Sex))), aes(x = Q, y = -log10(P), label = Label), size = 3, box.padding = 0.45, label.padding = 0.45, point.padding = 0, segment.color = 'grey50', max.overlaps = 10, min.segment.length = 0, seed = 1, force = 20) +
  annotate("text", x = max(all_metas_pvalues_adj$Q), y = -log10(bonf_thresh) + 0.2, label = "Bonferroni threshold", color = "red", hjust = 1, vjust = 0) +
  annotate("text", x = max(all_metas_pvalues_adj$Q), y = -log10(bh_thresh) + 0.2, label = "BH threshold", hjust = 1, vjust = 0) +
  theme(plot.title = element_text(size = 16))

ggsave(pvalue_plot, filename = file.path(meta_dir, "all_models_pop_meta_pvalues.png"), height = 5, width = 5, dpi = 120)





############# FACET PLOT ############# 


##### read in population info and pop compound estimates
pops <- readRDS("InfoTables/DrosEU_Populations.rds")
meta_pops <- readRDS(file.path(meta_dir, "all_models_pop_meta_compound_estimates.rds")) 

##### define variables for plotting
meta_country <- inner_join(meta_pops, pops$by_lat) %>%
  mutate(Group = paste(Trait, Sex, Models)) %>%
  group_by(Group) %>% 
  mutate(y = as.numeric(Country_code),
         xpos = min(LLEst),
         ypos_qp = 11.2,
         ypos_n = 10.2,
         N_plot = paste("N ==", round(N_lab_av, 1))) %>%
  filter(N_lab_av >= 2)

##### get model pvalues, add Group and format statistics
meta_pops_pvalues <- readRDS(file.path(meta_dir, "all_models_pop_meta_pvalues.rds")) %>% 
  mutate(Group = paste(Trait, Sex, Models),
         Q_plot = paste0("italic(Q) == ", round(Q, 2)),
         P_bonf_plot = ifelse(P_bonf < 0.001, "italic(p) < 0.001", paste0("italic(p) == ", round(P_bonf, 3))))


##### prepare text plotting
stats_text <- inner_join(select(meta_pops_pvalues, Group, Q_plot, P_bonf_plot), select(meta_country, Group, N_plot, xpos, ypos_qp, ypos_n) %>% distinct())

##### facet plot
p_meta_CI_facet <- ggplot(data = meta_country, aes(x = Estimate, y = y, color = Country_code)) +
  facet_wrap(Group ~ ., scales = "free_x", ncol = 5) +
  theme_bw() +
  geom_point(shape = 15, size = 2) +
  geom_errorbarh(aes(xmin = LLEst, xmax = ULEst), height = 0) +
  droseu_color_scale_country +
  scale_y_continuous(name = "Population", breaks = 1:max(meta_country$y), labels = unique(meta_country$Country_code)) +
  labs(x = "Population summary effect", y = "Population", title = "Subgroup meta analyses results", subtitle = "Populations summary effects with 95% CI, Q and P values and average number of labs (N)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text(data = stats_text, aes(x = xpos, y = ypos_qp, label = paste(Q_plot, P_bonf_plot, sep = "~~")), color = "black", parse = T, hjust = 0, size = 3, vjust = 1) +
  geom_text(data = stats_text, aes(x = xpos, y = ypos_n, label = N_plot), color = "black", parse = T, hjust = 0, size = 3, vjust = 1) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.position = "none")

##### save facet plot
ggsave(p_meta_CI_facet, filename = "MetaAnalyses/all_models_pop_meta_summary_effect.pdf", width = 8.27, height = 11.69)
ggsave(p_meta_CI_facet, filename = "MetaAnalyses/all_models_pop_meta_summary_effect.png", width = 8.27, height = 11.69, dpi = 120)



############# INDVIDUAL PLOTS, LOW RES FOR HTML ############# 

# re usage of facet plot stats_text

##### list all metas
metas_pop_trait <- list.files(path = meta_dir, recursive = T, full.names = T, pattern = "lmers_pop_meta_compound_estimates.rds")

for (i in 1:length(metas_pop_trait)){
  fpath <- metas_pop_trait[i]
  # set output names
  p_out_pdf <- sub("compound_estimates.rds", "summary_effect.pdf", fpath)
  p_out_png <- sub(".pdf", ".png", p_out_pdf)
  # read compound estimates in and define variables
  m <- readRDS(fpath) %>%
    inner_join(pops$by_lat) %>%
    mutate(Group = paste(Trait, Sex, Models),
           y = as.numeric(Country_code))
  # get stats defined before
  m <- inner_join(m, stats_text)
  # plot relevant metas
  if (nrow(m) > 0) {
    title_text <- paste(unique(m$Trait), unique(m$Sex), "summary effect with 95% CI")
    qvalue <- unique(m$Q_plot)
    pvalue <- unique(m$P_bonf_plot)
    nvalue <- unique(m$N_plot)
    # plot
    p_meta_CI <- ggplot(data = m, aes(x = Estimate, y = y, color = Country_code)) +
      theme_bw(16) +
      geom_point(shape = 15, size = 6) +
      geom_errorbarh(aes(xmin = LLEst, xmax = ULEst), height = 0, lwd = 1.5) +
      droseu_color_scale_country +
      scale_y_continuous(
        name = "Population", breaks = 1:max(m$y),
        labels = unique(m$Country_code)
      ) +
      labs(x = "Population summary effect", y = "Population", title = title_text) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      geom_text(aes(x = unique(xpos), y = unique(ypos_n)), label = paste(qvalue, pvalue, nvalue, sep = "~~"), color = "black", parse = T, hjust = 0, size = 5, vjust = 1) +
      theme(plot.title = element_text(size = 16),
            legend.position = "none")
    # save plot
    ggsave(p_meta_CI, filename = p_out_pdf, width = 5, height = 5)
    ggsave(p_meta_CI, filename = p_out_png, width = 5, height = 5, dpi = 120)
  }
}


############# INDVIDUAL PLOTS FOR PAPER FIGURES ############# 

# re usage of facet plot stats_text


trait_names <- read.csv("InfoTables/trait_names.csv") %>%
  dplyr::select(Trait_handle, Trait, Trait_name_raw, Legend, Title) %>%
  rename(Trait_name = Trait_name_raw)


##### list all metas
metas_pop_trait <- list.files(path = meta_dir, recursive = T, full.names = T, pattern = "lmers_pop_meta_compound_estimates.rds")

for (i in 1:length(metas_pop_trait)){
  fpath <- metas_pop_trait[i]
  # set output names
  p_out_pdf <- sub("compound_estimates.rds", "summary_effect_high_res.pdf", fpath)
  p_out_png <- sub(".pdf", ".png", p_out_pdf)
  # read compound estimates in and define variables
  m <- readRDS(fpath) %>%
    inner_join(pops$by_lat) %>%
    mutate(Group = paste(Trait, Sex, Models),
           y = as.numeric(Country_code))
  # get stats defined before
  m <- inner_join(m, stats_text)
  m <- inner_join(m, trait_names)
  # plot relevant metas
  if (nrow(m) > 0) {
    title_text <- paste(unique(m$Title), unique(m$Sex), sep = " - ")
    title_text <- sub(" - F", " - Females", title_text)
    title_text <- sub(" - M", " - Males", title_text)
    title_text <- sub(" - NA", "", title_text)
    sub_title_text <- bquote(italic("Q") == .(m$Q) ~ italic("p") == .(m$P) ~ N == .(m$N_lab_av))
    qvalue <- unique(m$Q_plot)
    pvalue <- unique(m$P_bonf_plot)
    nvalue <- unique(m$N_plot)
    # plot
    p_meta_CI <- ggplot(data = m, aes(x = Estimate, y = Country_code, color = Country_code)) +
      theme_classic() +
      geom_point(shape = 15, size = 6) +
      geom_errorbarh(aes(xmin = LLEst, xmax = ULEst), height = 0, lwd = 1.5) +
      droseu_color_scale_country +
      labs(
        x = "Summary effect with 95% CI", y = "Population",
        title = title_text, subtitle = sub_title_text
      ) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(legend.position = "none",
        axis.text = element_text(size = 22),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 22),
        plot.subtitle = element_text(size = 18)
      )
    # save plot
    ggsave(p_meta_CI, filename = p_out_pdf, width = 6, height = 8)
    ggsave(p_meta_CI, filename = p_out_png, width = 6, height = 8)
  }
}
