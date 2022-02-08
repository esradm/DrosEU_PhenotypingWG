
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
    dplyr::select(-c(SE, Value)) %>%
    pivot_wider(names_from = c(Lab, Sex), names_sep = "_", values_from = Estimate) %>%
    relocate(contains("_F"), .after = Population) %>%
    mutate(Population = factor(Population, levels = c("YE","RE","GI","MU","MA","UM","KA","VA","AK"))) %>% arrange(Population) %>% 
    dplyr::select(!Population)
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
  r <- abs(cor(x, y, use = "complete.obs"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
  cex.lab <- 2
}

# define the appearance of the scatter plots in pairs2
panScatterPlot <- function(x, y, bg.col = c("#cab2d6", "#ff7f00", "#fdbf6f", "#e31a1c", "#fb9a99", "#33a02c", "#b2df8a", "#1f78b4", "#a6cee3")){
  points(x, y, pch = 21, cex = 1.8, bg = bg.col)
  abline(lm(y~x))
}

# run longToWide and prepScatterPlotMatrix and plot the matrix with pairs2
scatterPlotMatrix <- function(x, sex) {
  x <- longToWide(x)
  if (!missing(sex)) { x <- prepScatterPlotMatrix(x, sex = sex) }
  else { x <- prepScatterPlotMatrix(x) }
  par(las = 2)
  pairs2(x, lower.panel = panScatterPlot, upper.panel = panCor, oma = c(6,4,2,4), ax.labels = TRUE, ax.ticks = TRUE, gap = 1)
  legend("bottom", xjust = 0.5, inset = -ncol(x)*0.012, legend = c("YE","RE","GI","MU","MA","UM","KA","VA","AK"), pch = 21, pt.bg = c("#cab2d6", "#ff7f00", "#fdbf6f", "#e31a1c", "#fb9a99", "#33a02c", "#b2df8a", "#1f78b4", "#a6cee3"), horiz = T, cex = 0.8, bty = "n", xpd = T)}

