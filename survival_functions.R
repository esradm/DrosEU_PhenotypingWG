



survProp <- function(x, cuts = seq(0, 120, 10)) {
  cols_to_keep <- c("Supervisor.PI", "Population", "Line", 
                    "Sex", "ReplicateVial", "ReplicateCage")
  prop <- x[1, colnames(x) %in% cols_to_keep]
  prop <- rbind(prop, prop[rep(1, length(cuts)-1), ])
  prop$AgeAtDeath = cuts
  prop$PropSurv = NA
  ninds <- nrow(x)
  for (i in 1:nrow(prop)){
    dead <- sum(x$AgeAtDeath <= prop$AgeAtDeath[i])
    prop$PropSurv[i] <- 1 - (dead / ninds)
  }
  return(prop)
}


