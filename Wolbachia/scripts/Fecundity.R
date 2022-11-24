
library(afex)
library(tidyverse)
library(sjPlot)
library(sjmisc)
library(car)
library(see)

## set working directory
# setwd("/media/inter/mkapun/projects/DrosEU_PhenotypingWG/")
setwd("/Users/martinkapun/Documents/GitHub/DrosEU_PhenotypingWG/")

## Function to plot lineplots with Standard Deviations using ggplot
data_summary <- function(x) {
    m <- mean(x)
    se <- sd(x)
    # se <- sd(x)/sqrt(length(x))
    ymin <- m - se
    ymax <- m + se
    return(c(y = m, ymin = ymin, ymax = ymax))
}

## Function to plot lineplots with Standard Errors using ggplot
data_summary_se <- function(x) {
    m <- mean(x)
    # se <- sd(x)
    se <- sd(x) / sqrt(length(x))
    ymin <- m - se
    ymax <- m + se
    return(c(y = m, ymin = ymin, ymax = ymax))
}

## make output directory
dir.create("Wolbachia/results")

PH <- "Fecundity"
TRt <- "NumberOfAdultsEclosed"

## create output directory
out_dir <- "Wolbachia/results"
dir.create(file.path(out_dir, PH), showWarnings = F)

## start output file for stats
# sink(paste0(out_dir, "/", PH, "/", TRt, ".txt"))
# print(paste0("__________________", TRt, "_______________"))

## Read Raw Data Table
DATA <- read.table(paste0("Wolbachia/data/", PH, "/", PH, ".txt"),
    header = T,
    comment.char = "",
    sep = "\t"
)

## Define as factors
DATA$Wolbachia <- as.factor(DATA$Wolbachia)
DATA$Country <- as.factor(DATA$Country)
DATA$Lab <- as.factor(DATA$Supervisor.PI)

## Exclude Finland and Russia, since all lines are Wolbachia positive there.
DATA <- DATA %>%
    filter(!DATA$Country %in% c("Finland", "Russia"))

# prepare final table for printing

print("______ Table with line counts per country_________")
print(tmp2)

## reorder countries by longitude
DATA$Country <- with(DATA, reorder(Country, Longitude))


## include line connecting the labs
DATA.p <- ggplot(DATA, aes(x = Lab, y = Trait, col = Wolbachia)) +
    stat_summary(
        fun.data = data_summary,
        aes(
            group = Wolbachia,
            color = Wolbachia
        ),
        geom = "line",
        size = 0.2,
        show.legend = F
    ) +
    stat_summary(
        fun.data = data_summary,
        aes(
            colour = Wolbachia,
            pch = Wolbachia
        ),
        size = 1
    ) +
    theme_bw() +
    facet_grid(~Line, scales = "free_y") +
    theme(axis.title.y = element_text(size = 20, angle = 90)) +
    theme(axis.title.x = element_text(size = 20, angle = 00)) +
    theme(axis.text = element_text(size = 10)) +
    theme(legend.text = element_text(size = 20)) +
    theme(legend.title = element_text(size = 20)) +
    theme(strip.text = element_text(size = 20)) +
    ylab(TRt) +
    xlab("Laboratory")

ggsave(paste0(out_dir, "/", PH, "/", TRt, "_byLine.pdf"),
    DATA.p,
    width = 60,
    height = 5,
    limitsize = FALSE
)

ggsave(paste0(out_dir, "/", PH, "/", TRt, "_byLine.png"),
    DATA.p,
    width = 60,
    height = 5,
    limitsize = FALSE
)


sink(paste0(out_dir, "/", PH, "/", TRt, "_byLine.txt"))

# Test for the effect of Wolbachia on the Trait variance
tmp <- DATA %>%
    group_by(Line, Wolbachia, Country, Lab) %>%
    summarise(Mean = mean(Trait)) %>%
    spread(Lab, Mean) %>%
    mutate(SD = abs(Billeter - Fricke))

Test <- lm(SD ~ Country * Wolbachia,
    data = tmp
)

print(Anova(Test,
    type = 3
))

print("Test by Lab")
for (i in c("Billeter", "Fricke")) {
    tmp <- DATA %>%
        filter(DATA$Supervisor.PI == i)

    Test <- lmer(NumberOfAdultsEclosed ~ Country * Wolbachia + (1 | Line:Country), data = tmp)
    print(i)
    print(Anova(Test,
        type = 3
    ))
}


sink()
