

##### DEFINE FUNCTIONS
# these functions have been written based on the "introduction to meta-analysis" book

# computes study weights and summary effect under a fixed effect model.
# study: vector of studies IDs
# Y: vector of studies effect sizes
# V: vector of studies effect sizes variances
# returns a data.frame with all the computed statistics
summaryEffectFixedComp <- function(study, Y, V){
  W = 1 / V # studies weights
  W_rel = W / sum(W) * 100 # studies absolute weights
  M = sum(W * Y) / sum(W) # summary effect
  VM = 1 / sum(W) # estimated variance for the summary effect
  SEM = sqrt(VM) # estimated SE for the summary effect 
  LLM <- M - 1.96 * SEM # 95% lower limit for the summary effect 
  ULM <- M + 1.96 * SEM # 95% upper limit for the summary effect 
  Z <- M / SEM # Z-value for p-value computation to test null hypothesis (mean effect is 0)
  p1 <- 1 - pnorm(Z) # one-tailed test p-value
  p2 <- 2 * p1 # two-tailed test p-value
  data.frame(Study = study, Y, V, W, W_rel, M, VM, SEM, LLM, ULM, Z, p1, p2)
}

# computes study weights and summary effect under a random-effects model.
# study: vector of studies IDs
# Y: vector of studies effect size
# V: vector of within-study variances
# T2: single value representing the between-study variance
# star notation stands for random-effects version of the computed statistics
# returns a data.frame with all the computed statistics
summaryEffectRandComp <- function(study, Y, V, T2){
  Wstar = 1 / (V + T2) # studies weights
  Wstar_rel = Wstar / sum(Wstar) * 100 # studies absolute weights
  Mstar = sum(Wstar * Y) / sum(Wstar) # summary effect
  VMstar = 1 / sum(Wstar) # estimated variance for the summary effect
  SEMstar = sqrt(VMstar) # estimated SE for the summary effect 
  LLMstar <- Mstar - 1.96 * SEMstar # 95% lower limit for the summary effect 
  ULMstar <- Mstar + 1.96 * SEMstar # 95% upper limit for the summary effect 
  Zstar <- Mstar / SEMstar # Z-value for p-value computation to test null hypothesis (mean effect is 0)
  p1 <- 1 - pnorm(Zstar) # one-tailed test p-value
  p2 <- 2 * p1 # two-tailed test p-value
  data.frame(Study = study, Y, V, T2, Wstar, Wstar_rel, Mstar, VMstar, 
             SEMstar, LLMstar, ULMstar, Zstar, p1, p2)
}

# computes the tau^2 statistics that represents the between-study variance
# also computes different measures of heterogeneity and tests whether heterogeneity is significant
# Y: vector of studies effect size
# W: vector of studies weight
# returns a data.frame with all the computed statistics
T2Comp <- function(Y, V) {
  W = 1 / V
  Q = sum(W*Y^2) - sum(W*Y)^2 / sum(W) # the observed WSS, reflects total dispersion
  df = length(Y) - 1 # degrees of freedom, the expected WSS
  excess = Q - df # excess variation attributed to differences in true effects
  p = pchisq(Q, df = df, lower.tail = FALSE) # test if the heterogeneity is statistically significant
  C = sum(W) - sum(W^2) / sum(W) # puts T2 in the same metric as the effect size
  T2 = excess / C # estimate of variance between studies, variance of the true effect size
  if (T2 < 0) T2 = 0
  T = sqrt(T2) # estimate of SD of the true effect size
  # below are the confidence intervals for T2
  sw1 = sum(W)
  sw2 = sum(W^2)
  sw3 = sum(W^3)
  A = df + 2 * (sw1 - sw2/sw1) * T2 + (sw2 - 2 * (sw3/sw1) + sw2^2/sw1^2) * T2^2
  VT2 = 2*(A/C^2)
  SET2 = sqrt(VT2)
  if (Q > df + 1) B = 0.5 * (log(Q) - log(df)) / (sqrt(2*Q) - sqrt(2*df-1))
  if (Q <= df + 1) B = sqrt(1 / (2 * (df-1) * (1 - (1 / (3 * (df-1)^2)))))
  L = exp(0.5 * log(Q/df) - 1.96 * B)
  U = exp(0.5 * log(Q/df) + 1.96 * B)
  LLT2 =  df * (L^2 - 1) / C
  if (LLT2 < 0) LLT2 = 0
  ULT2 =  df * (U^2 - 1) / C
  if (ULT2 < 0) ULT2 = 0
  LLT = sqrt(LLT2)
  ULT = sqrt(ULT2)
  data.frame(Q, df, excess, p, C, T2, T, A, VT2, SET2, B, L, U, LLT2, ULT2, LLT, ULT)
}


# computes the I^2 statistics that represents of kind of signal-to-noise ratio, a ratio of true heterogeneity to total variation in observed effects
# Q: the observed WSS
# df: degrees of freedom
# excess: the expected WSS
# returns a data.frame with all the computed statistics
I2Comp <- function(Q, df, excess) {
  I2 = excess / Q * 100
  I2 = (Q-df)/Q
  if (I2 < 0) I2 = 0
  if (Q > df + 1) B = 0.5 * (log(Q) - log(df)) / (sqrt(2*Q) - sqrt(2*df-1))
  if (Q <= df + 1) B = sqrt(1 / (2 * (df-1) * (1 - (1 / (3 * (df-1)^2)))))
  L = exp(0.5 * log(Q/df) - 1.96 * B)
  U = exp(0.5 * log(Q/df) + 1.96 * B)
  LLI2 = (L^2-1) / L^2 * 100
  if (LLI2 < 0) LLI2 = 0
  ULI2 = (U^2-1) / U^2 * 100
  if (ULI2 < 0) ULI2 = 0
  data.frame(Q, df, excess, I2, B, L, U, LLI2, ULI2)
}




metaAnalysisRandomModel <- function(studies.effects) {
  ### 1) calculate tau2 for each subgroup
  tau2 <- studies.effects %>%
    group_split(Population) %>%
    lapply(function(x) T2Comp(Y = x$Y, V = x$V)) %>%
    bind_rows() %>%
    mutate(Population = unique(studies.effects$Population)) %>%
    relocate(Population)
  ### 2) random effect model for each group. Caculates study weigths and summary effects with the ramdom model. Also calculates confidence intervals for summary effects
  summary_effects_random <- inner_join(studies.effects, dplyr::select(tau2, Population, T2)) %>%
    group_split(Population) %>%
    lapply(function(x) summaryEffectRandComp(study = x$Study, Y = x$Y, V = x$V, T2 = unique(x$T2))) %>%
    bind_rows() %>%
    mutate(Population = studies.effects$Population) %>%
    relocate(Population)
  ### 3) Effects comparison between sub groups
  ## 3.1) First method, Q-test based on analysis of variance
  # 3.1.1) calculate Qstar, the dispersion within each subgroup taking into account between study variance (T2). It can be done with the T2Comp function and the previously estimated T2 estimates)
  Qstar <- group_split(summary_effects_random, Population) %>%
    lapply(function(x) T2Comp(Y = x$Y, V = x$V + x$T2)) %>%
    bind_rows() %>%
    mutate(Population = unique(studies.effects$Population)) %>%
    relocate(Population) %>%
    dplyr::select(Population, Q) %>%
    dplyr::rename(Qstar = Q)
  # 3.1.2) random effect model for all studies
  summary_effects_random_all <- summaryEffectRandComp(study = summary_effects_random$Study, Y = summary_effects_random$Y, V = summary_effects_random$V, T2 = summary_effects_random$T2)
  # 3.1.3) Q for all studies
  Qstar_all <- T2Comp(Y = summary_effects_random_all$Y, 
                      V = summary_effects_random_all$V + summary_effects_random_all$T2) %>%
    mutate(Population = "combined") %>%
    dplyr::select(Population, Q) %>%
    dplyr::rename(Qstar = Q)
  Qstars <- bind_rows(Qstar, Qstar_all)
  # 3.1.4) put the data together
  Qstars <- bind_cols(
    bind_rows(summary_effects_random %>% 
                dplyr::select(Mstar, VMstar, SEMstar, LLMstar, ULMstar, Zstar, p2) %>% 
                distinct(),
              summary_effects_random_all %>% dplyr::select(Mstar, VMstar, SEMstar, LLMstar, ULMstar, Zstar, p2) %>%
                distinct()),
    bind_rows(Qstar, Qstar_all)) %>%
    relocate(Population)
  # 3.1.5) analysis of variance
  Qtest_var <- data.frame(Qwithin = Qstars$Qstar[Qstars$Population != "combined"] %>% sum) %>%
    mutate(Qbetween = Qstars$Qstar[Qstars$Population == "combined"] - Qwithin,
           df = nrow(Qstars)-2,
           p = pchisq(Qbetween, df, lower.tail = FALSE))
  ## 3.2) Second method, Q-test for heterogeneity
  Qtest_het <- T2Comp(Y = unique(summary_effects_random$Mstar), 
                      V = unique(summary_effects_random$VMstar)) %>%
    dplyr::select(Q, df, excess, p)
  
  
  list(summary_effects_random = inner_join(studies.effects, summary_effects_random), tau2 = tau2, Qstars = Qstars, Qtest_var = Qtest_var, Qtest_het = Qtest_het)

}





test_data_for_meta <- bind_cols(
  Population = c(c(rep("A", 5), rep("B", 5))),
  Study = c("Thornhill", "Kendall", "Vandamm", "Leonard", "Professor", "Jefferies", "Fremont", "Doyle", "Stella", "Thorwal"),
  Y = c(0.110, 0.224, 0.338, 0.451, 0.480, 0.440, 0.492, 0.651, 0.710, 0.740),
  V = c(0.01, 0.03, 0.02, 0.01, 0.01, 0.015, 0.020, 0.015, 0.025, 0.012))




