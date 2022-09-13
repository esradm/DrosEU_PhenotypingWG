
#################################################################################
################# META ANALYSES TO GET COMPOUND LINE ESTIMATES ##################
#################################################################################

 




##### quick plot for tarit correlations

library(GGally)
corF <- ggpairs(select(compound_effects_wide, contains(c("_F", "NA"))))
ggsave(corF, file = "~/Desktop/trait_cor_F.pdf", width = 14, height = 14)

corM <- ggpairs(select(compound_effects_wide, contains(c("_M", "NA"))))
ggsave(corM, file = "~/Desktop/trait_cor_M.pdf", width = 14, height = 14)







for (f in 1:length(lines_coefs)) {
  
  m <- readRDS(lines_coefs[f]) %>% rename(Estimate = Coef)
  m_effects <- makeEffects(m)
  
  out_rds <- sub("line_random", "line_compound_random", lines_coefs[f])
  out_txt <- sub(".rds", ".txt", out_rds)
  
  traits <- unique(m_effects$Trait)
  meta_trait_res <- list()
  meta_trait_out <- list()
  
  for(tr in 1:length(traits)) {
    trait <- filter(m_effects, Trait == traits[tr])
    sexes <- unique(trait$Sex)
    meta_sex_res <- list()
    meta_sex_out <- list()
    
    for (s in 1:length(sexes)) {
      msex <- metagen(data = filter(trait, Sex == sexes[s]), TE = Y, seTE = SE, studlab = Study, fixed = FALSE, random = TRUE, method.tau = "REML")
      mres <- update.meta(msex, subgroup = Line, tau.common = FALSE)
      mout <- data.frame(Trait = traits[tr], 
                         Population = str_sub(mres$bylevs, 1, 2), 
                         Line = mres$bylevs, Sex = sexes[s], 
                         Value = mres$TE.random.w, 
                         SE = mres$seTE.random.w, 
                         LLM = mres$lower.random.w, 
                         ULM = mres$upper.random.w, 
                         N_lab = mres$k.w)
      meta_sex_res[[s]] <- mres
      mnane <- str_match(out_rds, '([^/]+)(?:/[^/]+){0}$')[,1]
      names(meta_sex_res)[s] <- sub(".rds", paste0("_", traits[tr], "_", sexes[s]), mnane)
      meta_sex_out[[s]] <- mout
    }
    meta_trait_out[[tr]] <- bind_rows(meta_sex_out)
    meta_trait_res[[tr]] <- meta_sex_res
  }
  meta_trait_out <- bind_rows(meta_trait_out)
  meta_trait_res <- unlist(meta_trait_res, recursive=FALSE)
  
  saveRDS(meta_trait_out, out_rds)
  write.table(meta_trait_out, out_txt, row.names = F, quote = F, sep = "\t")
  capture.output(meta_trait_res, file = sub(".txt", "_meta.txt", out_txt)) 
  
}





