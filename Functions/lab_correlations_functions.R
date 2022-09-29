



# get the line overlap between labs for a given trait. Returns a list with a summary table, the ids of line that overlap, and the ids of lines used in each lab.
getLineOverlap <- function(x) {  
  n_unique <- length(unique(x$Line))
  byPI <- split(x, f = x$Supervisor.PI)
  unique_byPI <- lapply(byPI, function(x) as.character(unique(x$Line)))
  n_unique_byPI <- unlist(lapply(unique_byPI, length))
  overlap <- Reduce(intersect, unique_byPI)
  n_overlap <- length(overlap)
  line_summary <- data.frame(PI = names(n_unique_byPI), 
                             n_lines = n_unique_byPI, 
                             n_overlap, n_unique)
  rownames(line_summary) <- NULL
  list(line_summary = line_summary, overlap = overlap,
       unique_byPI = unique_byPI)
}



# extract model estimates and SE
getEstSE <- function(x, groups = "Population") {
  require(emmeans)
  df <- as.data.frame(emmeans(x, groups, mode = "asymp"))[,1:3]
  colnames(df)[2] <- "Estimate"
  df
}



# extract model estimates and SE
getEstStdErr2 <- function(x, groups = "Population") {
  require(lsmeans)
  df <- as.data.frame(lsmeans(x, groups, mode = "asymp"))[,1:3]
  colnames(df)[2] <- "Estimate"
  df
}

# combine fitted values and SE extracted from different models 
combineEst <- function(x) {
  for(i in (1:length(x))) {
    info <- str_split(names(x)[i], "_") %>% unlist
    if (length(info) == 5) info[1] <- paste(info[1], info[2], sep = "_")
    x[[i]] <- x[[i]] %>% 
      mutate(Trait = info[1], Lab = info[length(info)-1], Sex = info[length(info)-2]) %>%
      relocate(Trait, Lab, Sex) } 
  bind_rows(x) }


# combine fitted values and SE extracted from different models
combineEst2 <- function(x) {
  for(i in (1:length(x))) {
    info <- str_split(names(x)[i], "_") %>% unlist
    if (length(info) == 6) info[1] <- paste(info[1], info[2], sep = "_")
    x[[i]] <- x[[i]] %>% 
      mutate(Trait = info[1], Lab = info[length(info)-2], Sex = info[length(info)-3]) %>%
      relocate(Trait, Lab, Sex) 
    if (!unique(x[[i]]$Sex) %in% c("F", "M")) x[[i]]$Sex <- "F" 
    } 
  bind_rows(x) }

# combine fitted values and SE extracted from different models
combineEst3 <- function(x) {
  for(i in (1:length(x))) {
    info <- str_split(names(x)[i], "_") %>% unlist
    if (length(info) >=5) info[1] <- paste(info[1], info[2], sep = "_")
    info[1] <- sub("_F", "", info[1])
    info[1] <- sub("_M", "", info[1])
    x[[i]] <- x[[i]] %>% 
      mutate(Trait = info[1], Lab = info[length(info)-2], Sex = info[length(info)-3]) %>%
      relocate(Trait, Lab, Sex) 
    if (!unique(x[[i]]$Sex) %in% c("F", "M")) x[[i]]$Sex <- "NA" 
    if (length(unique(x[[i]]$Trait)) == 1 & unique(x[[i]]$Trait) %in% c("Dia", "Fec")) x[[i]]$Sex <- "F"
    if (length(unique(x[[i]]$Trait)) == 1 & grepl("Pgm", unique(x[[i]]$Trait))) x[[i]]$Sex <- "F"
    if (length(unique(x[[i]]$Trait)) == 1 & unique(x[[i]]$Trait) %in% c("LA")) x[[i]]$Sex <- "M"
  } 
  bind_rows(x) }


# extract model estimates and SE
getEstStdErr <-  function(x = lmer_model) {
  data.frame(Estimate = coef(summary(x))[,1], 
             SE = coef(summary(x))[,2])
}

# combine fitted values and SE extracted from different models 
combineFitted <- function(labs = rep(c("Flatt", "Parsch", "Pasyukova"), each = 2), 
                          sex = rep(c("F", "M"), 3), 
                          models = list(fla_fem_lmer1, fla_mal_lmer1, par_fem_lmer3, 
                                        par_mal_lmer3, pas_fem_lmer3, pas_mal_lmer3)){
  models_estimates <- list()
  for (i in 1:length(models)){
    models_estimates[[i]] <- bind_cols(Lab = labs[i], Sex = sex[i], 
                                       getEstStdErr(models[[i]]) %>% 
                                         rownames_to_column("Population"))}
  bind_rows(models_estimates) %>% 
    mutate(Population = str_replace(Population, "Population", ""), Value = "Fitted") %>%
    return()
}


# turn long table to wide for correlations matrices
longToWide <- function(fitted_values) {
  fitted_values %>% 
    dplyr::select(-c(SE)) %>%
    pivot_wider(names_from = c(Lab, Sex), names_sep = "_", values_from = Estimate) %>%
    relocate(contains("_F"), .after = Population) %>%
    mutate(Population = factor(Population, levels = c("YE","RE","GI","MU","MA","UM","KA","VA","AK"))) %>% arrange(Population) %>% 
    dplyr::select(!Population)
}



# turn long table to wide for correlations matrices
longToWide2 <- function(estimates) {
  ltw <- dplyr::select(estimates, -c(SE, Trait, Model, Predictor)) %>%
    pivot_wider(names_from = c(Lab, Sex), names_sep = "_", values_from = Estimate) %>%
    #mutate(Population = factor(Population, levels = c("YE","RE","GI","MU","MA","UM","KA","VA","AK"))) %>%
    arrange(Population) %>% 
    dplyr::select(!Population) 
  if ("Line" %in% colnames(ltw)) { dplyr::select(ltw, !Line) }
}





# modify colnames so they can be used as plot labels
prepScatterPlotMatrix <- function(x, sex = c("F", "M", "NA")) {
  sex <- paste0("_", sex)
  x <- dplyr::select(x, contains(sex))
  if (sex == "_F") names(x) <- str_replace(names(x), "_F", "\nFemales")
  if (sex == "_M") names(x) <- str_replace(names(x), "_M", "\nMales")
  if (sex == "_NA") names(x) <- str_replace(names(x), "_NA", "") 
  x
}

# modify colnames so they can be used as plot labels
prepScatterPlotMatrix2 <- function(x) {
  sex <- str_split(names(x), "_", simplify = T)[,2] %>% unique()
  if (sex == "F") names(x) <- str_replace(names(x), "_F", "\nFemales")
  if (sex == "M") names(x) <- str_replace(names(x), "_M", "\nMales")
  if (sex == "NA") names(x) <- str_replace(names(x), "_NA", "")
  x
}

# custom pair function
pairs2 <- function (x, labels, panel = points, ..., lower.panel = panel, 
                    upper.panel = panel, diag.panel = NULL, text.panel = textPanel, 
                    label.pos = 0.5 + has.diag/3, cex.labels = NULL, font.labels = 1, 
                    row1attop = TRUE, gap = 1, ax.labels = TRUE, ax.ticks = TRUE) {
  textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
                                                               y, txt, cex = cex, font = font)
  localAxis <- function(side, x, y, xpd, bg, col = NULL, main, 
                        oma, ...) {
    if (side%%2 == 1) 
      Axis(x, side = side, xpd = NA, ...)
    else Axis(y, side = side, xpd = NA, ...)
  }
  localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
  localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
  localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
  localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
  dots <- list(...)
  nmdots <- names(dots)
  if (!is.matrix(x)) {
    x <- as.data.frame(x)
    for (i in seq_along(names(x))) {
      if (is.factor(x[[i]]) || is.logical(x[[i]])) 
        x[[i]] <- as.numeric(x[[i]])
      if (!is.numeric(unclass(x[[i]]))) 
        stop("non-numeric argument to 'pairs'")
    }
  }
  else if (!is.numeric(x)) 
    stop("non-numeric argument to 'pairs'")
  panel <- match.fun(panel)
  if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
    lower.panel <- match.fun(lower.panel)
  if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
    upper.panel <- match.fun(upper.panel)
  if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) 
    diag.panel <- match.fun(diag.panel)
  if (row1attop) {
    tmp <- lower.panel
    lower.panel <- upper.panel
    upper.panel <- tmp
    tmp <- has.lower
    has.lower <- has.upper
    has.upper <- tmp
  }
  nc <- ncol(x)
  if (nc < 2) 
    stop("only one column in the argument to 'pairs'")
  has.labs <- TRUE
  if (missing(labels)) {
    labels <- colnames(x)
    if (is.null(labels)) 
      labels <- paste("var", 1L:nc)
  }
  else if (is.null(labels)) 
    has.labs <- FALSE
  oma <- if ("oma" %in% nmdots) 
    dots$oma
  else NULL
  main <- if ("main" %in% nmdots) 
    dots$main
  else NULL
  if (is.null(oma)) {
    oma <- c(4, 4, 4, 4)
    if (!is.null(main)) 
      oma[3L] <- 6
  }
  opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
  on.exit(par(opar))
  dev.hold()
  on.exit(dev.flush(), add = TRUE)
  for (i in if (row1attop) 
    1L:nc
    else nc:1L) for (j in 1L:nc) {
      localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
                type = "n", ...)
      if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
        box()
        # edited here...
        #           if (i == 1 && (!(j%%2) || !has.upper || !has.lower)) 
        #           localAxis(1 + 2 * row1attop, x[, j], x[, i], 
        #                       ...)
        # draw x-axis
        if (i == nc & j != nc) 
          localAxis(1, x[, j], x[, i], labels = ax.labels, tick = ax.ticks, ...)
        # draw y-axis
        if (j == 1 & i != 1) 
          localAxis(2, x[, j], x[, i], labels = ax.labels, tick = ax.ticks, ...)
        #           if (j == nc && (i%%2 || !has.upper || !has.lower)) 
        #             localAxis(4, x[, j], x[, i], ...)
        mfg <- par("mfg")
        if (i == j) {
          if (has.diag) 
            localDiagPanel(as.vector(x[, i]), ...)
          if (has.labs) {
            par(usr = c(0, 1, 0, 1))
            if (is.null(cex.labels)) {
              l.wid <- strwidth(labels, "user")
              cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
            }
            text.panel(0.5, label.pos, labels[i], cex = cex.labels, 
                       font = font.labels)
          }
        }
        else if (i < j) 
          localLowerPanel(as.vector(x[, j]), as.vector(x[, 
                                                         i]), ...)
        else localUpperPanel(as.vector(x[, j]), as.vector(x[, 
                                                            i]), ...)
        if (any(par("mfg") != mfg)) 
          stop("the 'panel' function made a new plot")
      }
      else par(new = FALSE)
    }
  if (!is.null(main)) {
    font.main <- if ("font.main" %in% nmdots) 
      dots$font.main
    else par("font.main")
    cex.main <- if ("cex.main" %in% nmdots) 
      dots$cex.main
    else par("cex.main")
    mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
  }
  invisible(NULL)
}

# define the appearance of the correlation coefficients in pairs2
panCor <- function(x, y, digits=2, prefix="", cex.cor){
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y, use = "complete.obs"), 2)
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  #if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  #text(0.5, 0.5, txt, cex = cex.cor * r)
  text(0.5, 0.5, txt, cex = 2.5)
  cex.lab <- 2
}



# define the appearance of the scatter plots in pairs2
panScatterPlot <- function(x, y, bg.col = met.brewer("Johnson", 9)){
  if (length(x) > 9) { 
    bg.col = "black"
    cex = 1.4 } 
  else {
    cex = 1.8
  }
  points(x, y, pch = 21, cex = cex, bg = bg.col)
  abline(lm(y~x))
}

# run longToWide and prepScatterPlotMatrix and plot the matrix with pairs2
# requires longToWide2(), prepScatterPlotMatrix2(), pairs2(), panScatterPlot(), panCor()
scatterPlotMatrix <- function(x, sex) {
  require(MetBrewer)
  x <- longToWide2(x)
  if (!missing(sex)) { x <- prepScatterPlotMatrix2(x, sex = sex) }
  else { x <- prepScatterPlotMatrix2(x) }
  par(las = 2)
  par(cex.axis = 1.6)
  pairs2(x, lower.panel = panScatterPlot, upper.panel = panCor, oma = c(6.5,4.5,2,4), ax.labels = TRUE, ax.ticks = TRUE, gap = 1)
  if (nrow(x) <= 9) {
    legend("bottom", xjust = 0.5, inset = -ncol(x)*0.012, legend = c("AK", "GI", "KA", "MA", "MU", "RE", "UM", "VA", "YE"), pch = 21, pt.cex = 2, pt.bg = met.brewer("Johnson", 9), horiz = T, cex = 0.8, bty = "n", xpd = T) }
}








