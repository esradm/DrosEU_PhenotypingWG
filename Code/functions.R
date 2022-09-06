# functions
# source("Code/functions.R")



# list.files wrapper to read in RDS files an turn them into list elements whose names are extracted from file paths.
rdsBatchReaderToList <- function(...) {
  temp = list.files(...)
  tnames <- str_split(temp, "/", simplify = T)[,3]
  tnames <- str_replace(tnames, ".rds", "")
  tlist <- lapply(temp, readRDS)
  names(tlist) <- tnames
  return(tlist)
}

# compute pairwise diff between Population using emmeans and multcomp::cld.
compTukeyCLD <- function(x) {
  require(emmeans)
  require(multcomp)
  em <- emmeans(x, pairwise ~ Population, mode = "asymp", adjust = "tukey")
  let <- cld(em, Letters = letters, alpha = 0.05) %>%
    as.data.frame %>% mutate(Population = as.factor(Population)) %>%
    arrange(Population) %>% dplyr::rename(cld = .group)
  list(contrasts = em, letters = let)
}


# plot linear model residuals. Takes file path of a RDS list object as input and produces a qqplot and a histogram for each element of the list.
plotResiduals <- function(f) {
  p <- str_match(f, '(.*[^/]+)(?:/[^/]+){1}$')[,2]
  dir.create(file.path(p, "by_lab_models_residuals"), showWarnings = F)
  m <- readRDS(f)
  n <- names(m)
  f_out <- file.path(p, "by_lab_models_residuals", n)
  qq_out_png <- paste0(f_out, "_qq_plot_residuals.png")
  hist_out_png <- paste0(f_out, "_hist_residuals.png")
  for (j in 1:length(n)) {
    png(filename = qq_out_png[j], height = 2100, width = 2100, res = 300)
    qqnorm(resid(m[[j]]), main = n[j])
    qqline(resid(m[[j]]))
    dev.off()
    png(filename = hist_out_png[j], height = 2100, width = 2100, res = 300)
    hist(resid(m[[j]]), main = n[j], xlab = "Residuals")
    dev.off() }
}

# output models statistics (summary, anova if lm, tukey if Pop) by trait
outputModelsStats <- function(f) {
  s_out_rds <- sub(".rds", "_summary.rds", f)
  s_out_txt <- sub(".rds", "_summary.txt", f)
  a_out_rds <- sub(".rds", "_anova.rds", f)
  a_out_txt <- sub(".rds", "_anova.txt", f)
  t_out_rds <- sub(".rds", "_tukey.rds", f)
  t_out_txt <- sub(".rds", "_tukey.txt", f)
  m <- readRDS(f)
  s <- lapply(m, summary)
  saveRDS(s, file = s_out_rds)
  capture.output(s, file = s_out_txt)
  if (class(m[[1]]) != "glmerMod") {
    a <- lapply(m, anova)
    saveRDS(a, file = a_out_rds)
    capture.output(a, file = a_out_txt) }
  if (grepl("lmers_pop.rds", f)) {
    t <- lapply(m, compTukeyCLD)
    saveRDS(t, file = t_out_rds)
    capture.output(t, file = t_out_txt) }
}

# output models statistics (summary, anova if lm, tukey if Pop) by trait and lab
outputModelsStatsLab <- function(f) {
  p <- str_match(f, '(.*[^/]+)(?:/[^/]+){1}$')[,2]
  dir.create(file.path(p, "by_lab_txt_output"), showWarnings = F)
  dir.create(file.path(p, "by_lab_rds_output"), showWarnings = F)
  m <- readRDS(f)
  n <- names(m)
  f_out_txt <- file.path(p, "by_lab_txt_output", n)
  f_out_rds <- file.path(p, "by_lab_rds_output", n)
  s_out_rds <- paste0(f_out_rds, "_summary.rds")
  s_out_txt <- paste0(f_out_txt, "_summary.txt")
  a_out_rds <- paste0(f_out_rds, "_anova.rds")
  a_out_txt <- paste0(f_out_txt, "_anova.txt")
  t_out_rds <- paste0(f_out_rds, "_tukey.rds")
  t_out_txt <- paste0(f_out_txt, "_tukey.txt")
  for (j in 1:length(n)) {
    s <- summary(m[[j]])
    saveRDS(s, file = s_out_rds[j])
    capture.output(s, file = s_out_txt[j])
    if (class(m[[j]]) != "glmerMod") {
      a <- anova(m[[j]])
      saveRDS(a, file = a_out_rds[j])
      capture.output(a, file = a_out_txt[j]) }
    if (grepl("lmers_pop.rds", f)) {
      t <- compTukeyCLD(m[[j]])
      saveRDS(t, file = t_out_rds[j])
      capture.output(t, file = t_out_txt[j]) }
  }
}

# get model estimates (fitted values) and their SE using emmeans package
getEstSE <- function(x, groups = "Population") {
  require(emmeans)
  df <- as.data.frame(emmeans(x, groups, mode = "asymp"))[,1:3]
  colnames(df)[2] <- "Estimate"
  df
}

# combine model estimates by trait
combineEstSE <- function(x) {
  for (lab in 1:length(x)) {
    info <- str_split(names(x)[lab], "_") %>% unlist
    if (length(info) >= 5) info[1] <- paste(info[1], info[2], sep = "_")
    info[1] <- sub("_F", "", info[1])
    info[1] <- sub("_M", "", info[1])
    es <- x[[lab]] %>%
      mutate(Model = info[length(info)-1], Predictor = info[length(info)], Trait = info[1], Lab = info[length(info)-2], Sex = info[length(info)-3]) %>% relocate(Model, Predictor, Trait, Lab, Sex)
    if (!unique(es$Sex) %in% c("F", "M")) es$Sex <- "NA"
    if (length(unique(es$Trait)) == 1 & unique(es$Trait) %in% c("Dia", "Fec")) es$Sex <- "F"
    if (length(unique(es$Trait)) == 1 & grepl("Pgm_", unique(es$Trait))) es$Sex <- "F"
    if (length(unique(es$Trait)) == 1 & grepl("LA_", unique(es$Trait))) es$Sex <- "B"
    x[[lab]] <- es }
  bind_rows(x) }


# combine model P values
combinePValues <- function(x) {
  pvals <- list()
  for (a in 1:length(x)) {
    info <- str_split(names(x)[a], "_") %>% unlist
    if (length(info) >= 5) info[1] <- paste(info[1], info[2], sep = "_")
    info[1] <- sub("_F", "", info[1])
    info[1] <- sub("_M", "", info[1])
    info <- data.frame(Lab = info[length(info)-2],
                       Trait = info[1],
                       Sex = info[length(info)-3],
                       Model = info[length(info)-1],
                       Predictor = info[length(info)])
    info$P <- ifelse(info$Model == "lmer", x[[a]]$P[1], x[[a]]$P[2])
    if (!info$Sex %in% c("F", "M")) info$Sex <- "NA"
    if (info$Trait %in% c("Dia", "Fec")) info$Sex <- "F"
    if (grepl("Pgm_", info$Trait)) info$Sex <- "F"
    if (grepl("LA_", info$Trait)) info$Sex <- "B"
    pvals[[a]] <- info }
  bind_rows(pvals) }

  # combine model P values for Wolbachia
  combinePValuesWolb <- function(x) {
    pvals <- list()
    for (a in 1:length(x)) {
      info <- str_split(names(x)[a], "_") %>% unlist
      if (length(info) >= 5) info[1] <- paste(info[1], info[2], sep = "_")
      info[1] <- sub("_F", "", info[1])
      info[1] <- sub("_M", "", info[1])
      info <- data.frame(Lab = info[length(info)-2],
                         Trait = info[1],
                         Sex = info[length(info)-3],
                         Model = info[length(info)-1],
                         Predictor = "Wolbachia")
      ## skip glmers for now
      if(info$Model == "lmer"){
      info$P <- x[[a]]$P[2]
      if (!info$Sex %in% c("F", "M")) info$Sex <- "NA"
      if (info$Trait %in% c("Dia", "Fec")) info$Sex <- "F"
      if (grepl("Pgm_", info$Trait)) info$Sex <- "F"
      if (grepl("LA_", info$Trait)) info$Sex <- "NA"
      pvals[[a]] <- info }
    }
    bind_rows(pvals) }



# combine estimates and format them to run meta analyses. Populations are ordered by latitude
makeEffects <- function(x) {
  x <- mutate(x, Population = factor(Population, levels = c("YE","RE","GI","MU","MA","UM","KA","VA","AK")), Lab = as.factor(Lab), V = SE^2, Study = paste(Population, Lab, sep = "_"))
  x <- relocate(x, Trait, Population, Sex, Lab, Study) %>%
    arrange(Population) %>% dplyr::rename(Y = Estimate)
  if ("Line" %in% colnames(x)) {
    x <- relocate(x, Trait, Population, Line, Sex, Lab, Study) %>%
      arrange(Population, Line) }
  x }
