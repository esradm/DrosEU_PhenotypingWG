

##### clean workspace
rm(list = ls())

##### libraries

library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(ggrepel)
library(cowplot)
library(ggpubr)


##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG_new/")

source("Code/functions.R")


pops <- read.csv("InfoTables/DrosEU_Coordinates.csv")
labs <- read.csv("InfoTables/DrosEU_cities.csv")



world <- ne_countries(scale = "small", returnclass = "sf")

pops_plot <- ggplot(data = world) +
    geom_sf(color = "grey30", fill = "grey80", size = 0.25) +
    coord_sf(xlim = c(-12, 38), ylim = c(30, 65), expand = FALSE) +
    annotation_scale(location = "bl", text_cex = 0.5, width_hint = 0.5) +
    geom_point(
        data = pops, aes(x = Longitude, y = Latitude, fill = Country_code),
        size = 2, color = "black", shape = 21
    ) +
    droseu_fill_scale_country +
    theme_classic(8) +
    geom_label_repel(
        data = pops, aes(x = Longitude, y = Latitude, label = Country_code), size = 2,
        box.padding = 0.26, min.segment.length = 0.2, point.padding = 0.2, seed = 1
    ) +
    xlab("Longitude") +
    ylab("Latitude") +
    #ggtitle("Sampling locations", subtitle = "9 populations") +
    theme(
        legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        panel.grid.major.y = element_line(size = 0.5)
    )


labs_plot <- ggplot(data = world) +
    geom_sf(color = "grey30", fill = "grey80", size = 0.25) +
    coord_sf(xlim = c(-85, 46), ylim = c(-45, 63), expand = FALSE) +
    geom_point(
        data = labs, aes(x = Longitude, y = Latitude), size = 1.5,
        fill = "red", color = "black", shape = 21
    ) +
    theme_classic(8) +
    xlab("Longitude") +
    ylab("Latitude") +
    #ggtitle("Contributing labs", subtitle = "17 countries, 26 research groups") +
    theme(
        legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        panel.grid.major.y = element_line(size = 0.5)
    )



plots <- plot_grid(
    pops_plot, labs_plot,
    labels = c("A", "B"), label_size = 12, rel_widths = c(1, 1.196),
    label_y = 0.95
)

ggsave(plots,
  filename = "Figures/figure1_v5.png",
  width = 6.3, height = 3.6
)

