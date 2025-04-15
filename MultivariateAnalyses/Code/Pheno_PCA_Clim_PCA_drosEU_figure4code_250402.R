rm(list=ls(all=TRUE))


### libraries
library(data.table)
library(ggplot2)
library(patchwork)
library(ggpubr)

### F9max - 30 day
### Load FMax Trait PC

workingDir = "C://Users//ewanh//Dropbox//Barcelona_IBE//DrosEU";
setwd(workingDir); 
getwd();

#load("F9_drosEU.RData")
load("FmaxP_drosEU.RData")

### trait PCA loadings
# F9_traitPC <- as.data.table(pca2$var$coord)
# F9_traitPC[,trait:=rownames(pca2$var$coord)]
# 
# F9_traitPC[,ordDim1:=rank(Dim.1)]
# F9_traitPC[,ordDim2:=rank(Dim.2)]

Fmp_traitPC <- as.data.table(pca4$var$coord)
Fmp_traitPC[,trait:=rownames(pca4$var$coord)]

Fmp_traitPC[,ordDim1:=rank(Dim.1)]
Fmp_traitPC[,ordDim2:=rank(Dim.2)]

Fmp_traitPC$trait<-c("Chill-coma recovery time","Cold shock mortality",
                    "Egg-to-adult development time","Diapause","Dry weight",
                    "Fecundity","Heat shock mortality","Life span",
                    "Total pigmentation",
                    "Starvation resistance","Thorax length","Wing area - Left",
                    "Viability")


### Load M9 Trait PC
# load("M9_drosEU.RData")
# 
# ### trait PCA loadings
# M9_traitPC <- as.data.table(pca1$var$coord)
# M9_traitPC[,trait:=rownames(pca1$var$coord)]
# 
# M9_traitPC[,ordDim1:=rank(Dim.1)]
# M9_traitPC[,ordDim2:=rank(Dim.2)]

### Load Day 30
load("all_traits_30d_WS.RData")
d30_F9_data <- F9_data
d30_Fmp_data <- F9maxP_data
d30_M9_data <- M9_data

### env PCA loadings
d30_envPC <- as.data.table(d30_pca$var$coord)
d30_envPC[,trait:=rownames(d30_pca$var$coord)]
d30_envPC[,ordDim1:=rank(Dim.1)]
d30_envPC[,ordDim2:=rank(Dim.2)]

d30_envPC$trait<-c("Longwave irradiance", "Frost days","Precipitation","Specific humidity",
                   "Relative humidity","Mean temp (2m)","Dew/frost point (2m)",
                   "Wet bulb temp (2m)","Max temp (2m)","Min temp (2m)","Temp range (2m)",
                   "Mean temp (earth skin)","Max temp (earth skin)","Min temp (earth skin)")

### load y30
load("all_traits_30y_WS.RData")
y30_F9_data <- F9_data
y30_Fmp_data <- F9maxP_data
y30_M9_data <- M9_data

### env PCA loadings
y30_envPC <- as.data.table(all_ann_clim_pca$var$coord)
y30_envPC[,trait:=rownames(all_ann_clim_pca$var$coord)]
y30_envPC[,ordDim1:=rank(Dim.1)]
y30_envPC[,ordDim2:=rank(Dim.2)]

y30_envPC$trait<-c("Longwave irradiance", "Frost days","Precipitation","Specific humidity",
                   "Relative humidity","Mean temp (2m)","Dew/frost point (2m)",
                   "Wet bulb temp (2m)","Max temp (2m)","Min temp (2m)","Temp range (2m)",
                   "Mean temp (earth skin)","Max temp (earth skin)","Min temp (earth skin)")

### general layout
layout <- "
  ABBB
  ABBB
  ABBB
  #CCC"


#####################
# Colour palette
#    Country   Color
#1   Austria #E38800
#2   Germany #F6C200
#4    Russia #095888
#7   Finland #A00E00
#8   Ukraine #0086A8
#13  Denmark #D04E00
#16 Portugal #7BA354
#18   Turkey #132B69
#71    Spain #B82E00

palette =c("#E38800","#D04E00","#A00E00","#F6C200","#7BA354","#095888","#B82E00","#132B69","#0086A8")


##################################
###
### Various objects for final plot
###

#load("all_traits_30y_WS.RData")

Fmp_traitPC_Sig <- subset(Fmp_traitPC, Dim.2 > 0.4 | Dim.2 < -0.4)
y30_envPC_sig <- subset(y30_envPC, Dim.2 > 0.65 | Dim.2 < -0.65)

Fmp_traitPC_Sig[,ordDim2:=rank(Dim.2)]
y30_envPC_sig[,ordDim2:=rank(Dim.2)]


flip_Fmp_30y <- 1
Fmp_30y_traitPC_load <- ggplot(data=Fmp_traitPC_Sig) +
  geom_segment( aes(x=0, xend=Dim.2, y=ordDim2, yend=ordDim2), arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
                linewidth =0.7, 
                lineend = "round",
                linejoin = "mitre")+
  geom_text(data=Fmp_traitPC_Sig[Dim.2<0], aes(y=ordDim2, x=.05, label=trait), size=4.5, angle=90, hjust=0) +
  geom_text(data=Fmp_traitPC_Sig[Dim.2>0], aes(y=ordDim2, x=-.05, label=trait), size=4.5, angle=90, hjust=1) +
  coord_flip() + xlim(-.8, .8) + ylim(0,5.5) +
  xlab("Loading") + ylab("") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=12),
        axis.title.y =element_text(size = 14),
        axis.ticks.y=element_line(linewidth=0.5),
        axis.ticks.length=unit(0.20,"cm"),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey70"))

Fmp_30y_phenoPC_load <- ggplot(data=y30_envPC_sig) +
  geom_segment(aes(x=0, xend=flip_Fmp_30y*Dim.2, y=ordDim2, yend=ordDim2), arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               linewidth =0.7, 
               lineend = "round",
               linejoin = "mitre")+
  geom_text(data=y30_envPC_sig[Dim.2<0], aes(y=ordDim2, x=flip_Fmp_30y*.05, label=trait), size=4.5, hjust=0) +
  geom_text(data=y30_envPC_sig[Dim.2>0], aes(y=ordDim2, x=flip_Fmp_30y*-.05, label=trait), size=4.5, hjust=1) +
  xlim(-0.9, 0.9)  + ylim(0, 4.5) +
  xlab("Loading") + ylab("") +
  theme_bw()+
  theme(axis.text.y=element_blank(),
                     axis.ticks.y=element_blank(),
                     axis.text.x=element_text(size=12),
                     axis.title.x =element_text(size = 14),
                     axis.ticks.x=element_line(linewidth=0.5),
                     axis.ticks.length=unit(0.20,"cm"),
                     panel.grid = element_blank(),
                     panel.border = element_rect(colour = "grey70"))

Fmp_30y_trait_env_plot <- ggplot(data=y30_Fmp_data, aes(x=flip_Fmp_30y*PC2_clim, y=PC2_F9maxP, color=Country)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_manual(values = palette) +
  ylim(-4.5, 4.5) +
  xlab("Climate PC2") + ylab("Phenotype PC2") +
  ggtitle("F13 PC2 vs 30 year PC2") +
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4))+
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4))+
  theme_bw()+
  theme(text = element_text(size = 16), 
                      #plot.title = element_blank(),
                      axis.title = element_text(size = 14),
                      axis.text = element_text(size = 14),
                      panel.grid.minor = element_blank())

F13_pc2_y30_pc2_mega <-
  Fmp_30y_traitPC_load + Fmp_30y_trait_env_plot + Fmp_30y_phenoPC_load +
  plot_layout(design=layout, 
              heights = c(2.5,1),
              widths = c(1,2.2))

#################
###
### For sup mat: F9 year

### trait PCA loadings

load("F9_drosEU.RData")
y30_F9_data <- F9_data

F9_traitPC <- as.data.table(pca2$var$coord)
F9_traitPC[,trait:=rownames(pca2$var$coord)]

F9_traitPC$trait<-c("Chill-coma recovery time","Cold shock mortality",
                    "Egg-to-adult development time","Dry weight",
                    "Heat shock mortality","Life span",
                    "Starvation resistance","Thorax length","Wing area - Left")


# 
#F9_traitPC[,ordDim1:=rank(Dim.1)]
#F9_traitPC[,ordDim2:=rank(Dim.2)]

F9_traitPC_Sig <- subset(F9_traitPC, Dim.2 > 0.4 | Dim.2 < -0.4)
y30_envPC_sig <- subset(y30_envPC, Dim.2 > 0.65 | Dim.2 < -0.65)

F9_traitPC_Sig[,ordDim2:=rank(Dim.2)]
y30_envPC_sig[,ordDim2:=rank(Dim.2)]



flip_F9_30y <- 1
F9_30y_traitPC_load <- ggplot(data=F9_traitPC_Sig) +
  geom_segment( aes(x=0, xend=Dim.2, y=ordDim2, yend=ordDim2), arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
                linewidth =0.7, 
                lineend = "round",
                linejoin = "mitre")+
  geom_text(data=F9_traitPC_Sig[Dim.2<0], aes(y=ordDim2, x=.05, label=trait), size=4.5, angle=90, hjust=0) +
  geom_text(data=F9_traitPC_Sig[Dim.2>0], aes(y=ordDim2, x=-.05, label=trait), size=4.5, angle=90, hjust=1) +
  coord_flip() + xlim(-.8, .8) + ylim(0,5.5) +
  xlab("Loading") + ylab("") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=12),
        axis.title.y =element_text(size = 14),
        axis.ticks.y=element_line(linewidth=0.5),
        axis.ticks.length=unit(0.20,"cm"),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey70"))

F9_30y_phenoPC_load <- ggplot(data=y30_envPC_sig) +
  geom_segment(aes(x=0, xend=flip_F9_30y*Dim.2, y=ordDim2, yend=ordDim2), arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               linewidth =0.7, 
               lineend = "round",
               linejoin = "mitre")+
  geom_text(data=y30_envPC_sig[Dim.2<0], aes(y=ordDim2, x=flip_F9_30y*.05, label=trait), size=4.5, hjust=0) +
  geom_text(data=y30_envPC_sig[Dim.2>0], aes(y=ordDim2, x=flip_F9_30y*-.05, label=trait), size=4.5, hjust=1) +
  xlim(-0.9, 0.9)  + ylim(0, 4.5) +
  xlab("Loading") + ylab("") +
  theme_bw()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=12),
        axis.title.x =element_text(size = 14),
        axis.ticks.x=element_line(linewidth=0.5),
        axis.ticks.length=unit(0.20,"cm"),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey70"))

F9_30y_trait_env_plot <- ggplot(data=y30_F9_data, aes(x=flip_F9_30y*PC2_clim, y=PC2_F9, color=Country)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_manual(values = palette) +
  ylim(-4.5, 4.5) +
  xlab("Climate PC2") + ylab("Phenotype PC2") +
  ggtitle("F9 PC2 vs 30 year PC2") +
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4))+
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4))+
  theme_bw()+
  theme(text = element_text(size = 16), 
        #plot.title = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        panel.grid.minor = element_blank())

F9_pc2_y30_pc2_mega <-
  F9_30y_traitPC_load + F9_30y_trait_env_plot + F9_30y_phenoPC_load +
  plot_layout(design=layout, 
              heights = c(2.5,1),
              widths = c(1,2.2))


#################
###
### For sup mat: M9 year

### trait PCA loadings

load("M9_drosEU.RData")
y30_M9_data <- M9_data

M9_traitPC <- as.data.table(pca1$var$coord)
M9_traitPC[,trait:=rownames(pca1$var$coord)]

M9_traitPC$trait<-c("Chill-coma recovery time","Cold shock mortality",
                    "Egg-to-adult development time","Dry weight",
                    "Heat shock mortality","Life span",
                    "Starvation resistance","Thorax length","Wing area - Left")

# 
#M9_traitPC[,ordDim1:=rank(Dim.1)]
#M9_traitPC[,ordDim2:=rank(Dim.2)]

M9_traitPC_Sig <- subset(M9_traitPC, Dim.2 > 0.4 | Dim.2 < -0.4)
#y30_envPC_sig <- subset(y30_envPC, Dim.2 > 0.65 | Dim.2 < -0.65)

M9_traitPC_Sig[,ordDim2:=rank(Dim.2)]
#y30_envPC_sig[,ordDim2:=rank(Dim.2)]


flip_M9_30y <- 1
M9_30y_traitPC_load <- ggplot(data=M9_traitPC_Sig) +
  geom_segment( aes(x=0, xend=Dim.2, y=ordDim2, yend=ordDim2), arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
                linewidth =0.7, 
                lineend = "round",
                linejoin = "mitre")+
  geom_text(data=M9_traitPC_Sig[Dim.2<0], aes(y=ordDim2, x=.05, label=trait), size=4.5, angle=90, hjust=0) +
  geom_text(data=M9_traitPC_Sig[Dim.2>0], aes(y=ordDim2, x=-.05, label=trait), size=4.5, angle=90, hjust=1) +
  coord_flip() + xlim(-.8, .8) + ylim(0,5.5) +
  xlab("Loading") + ylab("") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=12),
        axis.title.y =element_text(size = 14),
        axis.ticks.y=element_line(linewidth=0.5),
        axis.ticks.length=unit(0.20,"cm"),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey70"))

M9_30y_phenoPC_load <- ggplot(data=y30_envPC_sig) +
  geom_segment(aes(x=0, xend=flip_M9_30y*Dim.2, y=ordDim2, yend=ordDim2), arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               linewidth =0.7, 
               lineend = "round",
               linejoin = "mitre")+
  geom_text(data=y30_envPC_sig[Dim.2<0], aes(y=ordDim2, x=flip_M9_30y*.05, label=trait), size=4.5, hjust=0) +
  geom_text(data=y30_envPC_sig[Dim.2>0], aes(y=ordDim2, x=flip_M9_30y*-.05, label=trait), size=4.5, hjust=1) +
  xlim(-0.9, 0.9)  + ylim(0, 4.5) +
  xlab("Loading") + ylab("") +
  theme_bw()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=12),
        axis.title.x =element_text(size = 14),
        axis.ticks.x=element_line(linewidth=0.5),
        axis.ticks.length=unit(0.20,"cm"),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey70"))

M9_30y_trait_env_plot <- ggplot(data=y30_M9_data, aes(x=flip_M9_30y*PC2_clim, y=PC2_M9, color=Country)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_manual(values = palette) +
  ylim(-4.5, 4.5) +
  xlab("Climate PC2") + ylab("Phenotype PC2") +
  ggtitle("M9 PC2 vs 30 year PC2") +
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4))+
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4))+
  theme_bw()+
  theme(text = element_text(size = 16), 
        #plot.title = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        panel.grid.minor = element_blank())

M9_pc2_y30_pc2_mega <-
  M9_30y_traitPC_load + M9_30y_trait_env_plot + M9_30y_phenoPC_load +
  plot_layout(design=layout, 
              heights = c(2.5,1),
              widths = c(1,2.2))

########################################################################
###
############################################################
############################################################
###
### 30 day plots
###
#########################################################################
##################################
###
### F13
###

#load("all_traits_30d_WS.RData")

Fmp_traitPC_Sig <- subset(Fmp_traitPC, Dim.2 > 0.4 | Dim.2 < -0.4)
d30_envPC_sig <- subset(d30_envPC, Dim.2 > 0.65 | Dim.2 < -0.65)

Fmp_traitPC_Sig[,ordDim2:=rank(Dim.2)]
d30_envPC_sig[,ordDim2:=rank(Dim.2)]


flip_Fmp_30d <- 1
Fmp_30d_traitPC_load <- ggplot(data=Fmp_traitPC_Sig) +
  geom_segment( aes(x=0, xend=Dim.2, y=ordDim2, yend=ordDim2), arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
                linewidth =0.7, 
                lineend = "round",
                linejoin = "mitre")+
  geom_text(data=Fmp_traitPC_Sig[Dim.2<0], aes(y=ordDim2, x=.05, label=trait), size=4.5, angle=90, hjust=0) +
  geom_text(data=Fmp_traitPC_Sig[Dim.2>0], aes(y=ordDim2, x=-.05, label=trait), size=4.5, angle=90, hjust=1) +
  coord_flip() + xlim(-.8, .8) + ylim(0,5.5) +
  xlab("Loading") + ylab("") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=12),
        axis.title.y =element_text(size = 14),
        axis.ticks.y=element_line(linewidth=0.5),
        axis.ticks.length=unit(0.20,"cm"),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey70"))

Fmp_30d_phenoPC_load <- ggplot(data=d30_envPC_sig) +
  geom_segment(aes(x=0, xend=flip_Fmp_30d*Dim.2, y=ordDim2, yend=ordDim2), arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               linewidth =0.7, 
               lineend = "round",
               linejoin = "mitre")+
  geom_text(data=d30_envPC_sig[Dim.2<0], aes(y=ordDim2, x=flip_Fmp_30d*.05, label=trait), size=4.5, hjust=0) +
  geom_text(data=d30_envPC_sig[Dim.2>0], aes(y=ordDim2, x=flip_Fmp_30d*-.05, label=trait), size=4.5, hjust=1) +
  xlim(-0.9, 0.9)  + ylim(0, 4.5) +
  xlab("Loading") + ylab("") +
  theme_bw()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=12),
        axis.title.x =element_text(size = 14),
        axis.ticks.x=element_line(linewidth=0.5),
        axis.ticks.length=unit(0.20,"cm"),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey70"))

Fmp_30d_trait_env_plot <- ggplot(data=d30_Fmp_data, aes(x=flip_Fmp_30d*PC2_clim, y=PC2_F9maxP, color=Country)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_manual(values = palette) +
  ylim(-4.5, 4.5) +
  xlab("Climate PC2") + ylab("Phenotype PC2") +
  ggtitle("F13 PC2 vs 30 day PC2") +
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4))+
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4))+
  theme_bw()+
  theme(text = element_text(size = 16), 
        #plot.title = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        panel.grid.minor = element_blank())

F13_pc2_d30_pc2_mega <-
  Fmp_30d_traitPC_load + Fmp_30d_trait_env_plot + Fmp_30d_phenoPC_load +
  plot_layout(design=layout, 
              heights = c(2.5,1),
              widths = c(1,2.2))

#################
###
### For sup mat: F9 year

### trait PCA loadings

load("F9_drosEU.RData")
d30_F9_data <- F9_data

F9_traitPC <- as.data.table(pca2$var$coord)
F9_traitPC[,trait:=rownames(pca2$var$coord)]


F9_traitPC$trait<-c("Chill-coma recovery time","Cold shock mortality",
                    "Egg-to-adult development time","Dry weight",
                    "Heat shock mortality","Life span",
                    "Starvation resistance","Thorax length","Wing area - Left")


# 
#F9_traitPC[,ordDim1:=rank(Dim.1)]
#F9_traitPC[,ordDim2:=rank(Dim.2)]

F9_traitPC_Sig <- subset(F9_traitPC, Dim.2 > 0.4 | Dim.2 < -0.4)
d30_envPC_sig <- subset(d30_envPC, Dim.2 > 0.65 | Dim.2 < -0.65)

F9_traitPC_Sig[,ordDim2:=rank(Dim.2)]
d30_envPC_sig[,ordDim2:=rank(Dim.2)]


flip_F9_30d <- 1
F9_30d_traitPC_load <- ggplot(data=F9_traitPC_Sig) +
  geom_segment( aes(x=0, xend=Dim.2, y=ordDim2, yend=ordDim2), arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
                linewidth =0.7, 
                lineend = "round",
                linejoin = "mitre")+
  geom_text(data=F9_traitPC_Sig[Dim.2<0], aes(y=ordDim2, x=.05, label=trait), size=4.5, angle=90, hjust=0) +
  geom_text(data=F9_traitPC_Sig[Dim.2>0], aes(y=ordDim2, x=-.05, label=trait), size=4.5, angle=90, hjust=1) +
  coord_flip() + xlim(-.8, .8) + ylim(0,5.5) +
  xlab("Loading") + ylab("") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=12),
        axis.title.y =element_text(size = 14),
        axis.ticks.y=element_line(linewidth=0.5),
        axis.ticks.length=unit(0.20,"cm"),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey70"))

F9_30d_phenoPC_load <- ggplot(data=d30_envPC_sig) +
  geom_segment(aes(x=0, xend=flip_F9_30d*Dim.2, y=ordDim2, yend=ordDim2), arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               linewidth =0.7, 
               lineend = "round",
               linejoin = "mitre")+
  geom_text(data=d30_envPC_sig[Dim.2<0], aes(y=ordDim2, x=flip_F9_30d*.05, label=trait), size=4.5, hjust=0) +
  geom_text(data=d30_envPC_sig[Dim.2>0], aes(y=ordDim2, x=flip_F9_30d*-.05, label=trait), size=4.5, hjust=1) +
  xlim(-0.9, 0.9)  + ylim(0, 4.5) +
  xlab("Loading") + ylab("") +
  theme_bw()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=12),
        axis.title.x =element_text(size = 14),
        axis.ticks.x=element_line(linewidth=0.5),
        axis.ticks.length=unit(0.20,"cm"),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey70"))

F9_30d_trait_env_plot <- ggplot(data=d30_F9_data, aes(x=flip_F9_30d*PC2_clim, y=PC2_F9, color=Country)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_manual(values = palette) +
  ylim(-4.5, 4.5) +
  xlab("Climate PC2") + ylab("Phenotype PC2") +
  ggtitle("F9 PC2 vs 30 day PC2") +
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4))+
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4))+
  theme_bw()+
  theme(text = element_text(size = 16), 
        #plot.title = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        panel.grid.minor = element_blank())

F9_pc2_d30_pc2_mega <-
  F9_30d_traitPC_load + F9_30d_trait_env_plot + F9_30d_phenoPC_load +
  plot_layout(design=layout, 
              heights = c(2.5,1),
              widths = c(1,2.2))


#################
###
### For sup mat: M9 year

### trait PCA loadings

load("M9_drosEU.RData")
d30_M9_data <- M9_data

M9_traitPC <- as.data.table(pca1$var$coord)
M9_traitPC[,trait:=rownames(pca1$var$coord)]

M9_traitPC$trait<-c("Chill-coma recovery time","Cold shock mortality",
                    "Egg-to-adult development time","Dry weight",
                    "Heat shock mortality","Life span",
                    "Starvation resistance","Thorax length","Wing area - Left")

# 
#M9_traitPC[,ordDim1:=rank(Dim.1)]
#M9_traitPC[,ordDim2:=rank(Dim.2)]

M9_traitPC_Sig <- subset(M9_traitPC, Dim.2 > 0.4 | Dim.2 < -0.4)
d30_envPC_sig <- subset(d30_envPC, Dim.2 > 0.65 | Dim.2 < -0.65)

M9_traitPC_Sig[,ordDim2:=rank(Dim.2)]
d30_envPC_sig[,ordDim2:=rank(Dim.2)]


flip_M9_30d <- 1
M9_30d_traitPC_load <- ggplot(data=M9_traitPC_Sig) +
  geom_segment( aes(x=0, xend=Dim.2, y=ordDim2, yend=ordDim2), arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
                linewidth =0.7, 
                lineend = "round",
                linejoin = "mitre")+
  geom_text(data=M9_traitPC_Sig[Dim.2<0], aes(y=ordDim2, x=.05, label=trait), size=4.5, angle=90, hjust=0) +
  geom_text(data=M9_traitPC_Sig[Dim.2>0], aes(y=ordDim2, x=-.05, label=trait), size=4.5, angle=90, hjust=1) +
  coord_flip() + xlim(-.8, .8) + ylim(0,5.5) +
  xlab("Loading") + ylab("") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=12),
        axis.title.y =element_text(size = 14),
        axis.ticks.y=element_line(linewidth=0.5),
        axis.ticks.length=unit(0.20,"cm"),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey70"))

M9_30d_phenoPC_load <- ggplot(data=d30_envPC_sig) +
  geom_segment(aes(x=0, xend=flip_M9_30d*Dim.2, y=ordDim2, yend=ordDim2), arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               linewidth =0.7, 
               lineend = "round",
               linejoin = "mitre")+
  geom_text(data=d30_envPC_sig[Dim.2<0], aes(y=ordDim2, x=flip_M9_30d*.05, label=trait), size=4.5, hjust=0) +
  geom_text(data=d30_envPC_sig[Dim.2>0], aes(y=ordDim2, x=flip_M9_30d*-.05, label=trait), size=4.5, hjust=1) +
  xlim(-0.9, 0.9)  + ylim(0, 4.5) +
  xlab("Loading") + ylab("") +
  theme_bw()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=12),
        axis.title.x =element_text(size = 14),
        axis.ticks.x=element_line(linewidth=0.5),
        axis.ticks.length=unit(0.20,"cm"),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey70"))

M9_30d_trait_env_plot <- ggplot(data=d30_M9_data, aes(x=flip_M9_30d*PC2_clim, y=PC2_M9, color=Country)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_manual(values = palette) +
  ylim(-4.5, 4.5) +
  xlab("Climate PC2") + ylab("Phenotype PC2") +
  ggtitle("M9 PC2 vs 30 day PC2") +
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4))+
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4))+
  theme_bw()+
  theme(text = element_text(size = 16), 
        #plot.title = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        panel.grid.minor = element_blank())

M9_pc2_d30_pc2_mega <-
  M9_30d_traitPC_load + M9_30d_trait_env_plot + M9_30d_phenoPC_load +
  plot_layout(design=layout, 
              heights = c(2.5,1),
              widths = c(1,2.2))

library(cowplot)
obj <- plot_grid(F13_pc2_d30_pc2_mega,
          F9_pc2_d30_pc2_mega,
          M9_pc2_d30_pc2_mega,
          F13_pc2_y30_pc2_mega,
          F9_pc2_y30_pc2_mega,
          M9_pc2_y30_pc2_mega, nrow = 2)

#############################
###
### main fig for MS


Fmp_traitPC_Sig <- subset(Fmp_traitPC, Dim.2 > 0.4 | Dim.2 < -0.4)
y30_envPC_sig <- subset(y30_envPC, Dim.2 > 0.65 | Dim.2 < -0.65)

Fmp_traitPC_Sig[,ordDim2:=rank(Dim.2)]
y30_envPC_sig[,ordDim2:=rank(Dim.2)]


flip_Fmp_30y <- 1
Fmp_30y_traitPC_load <- ggplot(data=Fmp_traitPC_Sig) +
  geom_segment( aes(x=0, xend=Dim.2, y=ordDim2, yend=ordDim2), arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
                linewidth =0.7, 
                lineend = "round",
                linejoin = "mitre")+
  geom_text(data=Fmp_traitPC_Sig[Dim.2<0], aes(y=ordDim2, x=.05, label=trait), size=4.5, angle=90, hjust=0) +
  geom_text(data=Fmp_traitPC_Sig[Dim.2>0], aes(y=ordDim2, x=-.05, label=trait), size=4.5, angle=90, hjust=1) +
  coord_flip() + xlim(-.8, .8) + ylim(0,5.5) +
  xlab("Loading") + ylab("") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=12),
        axis.title.y =element_text(size = 14),
        axis.ticks.y=element_line(linewidth=0.5),
        axis.ticks.length=unit(0.20,"cm"),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey70"))

Fmp_30y_phenoPC_load <- ggplot(data=y30_envPC_sig) +
  geom_segment(aes(x=0, xend=flip_Fmp_30y*Dim.2, y=ordDim2, yend=ordDim2), arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               linewidth =0.7, 
               lineend = "round",
               linejoin = "mitre")+
  geom_text(data=y30_envPC_sig[Dim.2<0], aes(y=ordDim2, x=flip_Fmp_30y*.05, label=trait), size=4.5, hjust=0) +
  geom_text(data=y30_envPC_sig[Dim.2>0], aes(y=ordDim2, x=flip_Fmp_30y*-.05, label=trait), size=4.5, hjust=1) +
  xlim(-0.9, 0.9)  + ylim(0, 4.5) +
  xlab("Loading") + ylab("") +
  theme_bw()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=12),
        axis.title.x =element_text(size = 14),
        axis.ticks.x=element_line(linewidth=0.5),
        axis.ticks.length=unit(0.20,"cm"),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey70"))

Fmp_30y_trait_env_plot <- ggplot(data=y30_Fmp_data, aes(x=flip_Fmp_30y*PC2_clim, y=PC2_F9maxP, color=Country)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_manual(values = palette) +
  ylim(-4.5, 4.5) +
  xlab("Climate PC2") + ylab("Phenotype PC2") +
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4))+
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4))+
  theme_bw()+
  theme(text = element_text(size = 16), 
        #plot.title = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        panel.grid.minor = element_blank())

F13_pc2_y30_pc2_mega2 <-
  Fmp_30y_traitPC_load + Fmp_30y_trait_env_plot + Fmp_30y_phenoPC_load +
  plot_layout(design=layout, 
              heights = c(2.5,1),
              widths = c(1,2.2))

#################

