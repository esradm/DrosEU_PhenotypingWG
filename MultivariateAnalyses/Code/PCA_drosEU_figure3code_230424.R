rm(list=ls(all=TRUE))
library(dplyr)
library(FactoMineR)
library(factoextra)
library(ggforce)
library(cowplot)
library(tidyverse)
library(psych)
library(corrplot)
library(data.table)

# Load the data set into R
masterall <- read.csv("C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/all_models_line_meta_compound_random_coefs_wide_220916.csv")

##################################
###
### Preparing the masterall data
###

# Make variables for Country and the sexes
table(masterall$Population)
# AK GI KA MA MU RE UM VA YE 
# 22 15 20 20 20 17 19 20 20
Country <- c(rep("Turkey", 20),rep("Portugal", 17),rep("Spain", 15),rep("Germany", 20),rep("Austria", 20),rep("Ukraine", 19),rep("Denmark", 20),rep("Russia", 20),rep("Finland", 22))
table(masterall$Population)
# this is simply the number of observations
SexM <-c(rep("M", 173))
SexF <-c(rep("F", 173))

####
#### Formatting the data for the different versions of the PCA
####
#################################
# M9
masterallM9 <- masterall[,c("Population","Line","CCRT_M","CSM_M","DT_A_M","DW_M","HSM_M","LS_M","SR_M","TL_M","WA_L_M")]
dataall1 <- na.omit(cbind(Country,masterallM9))
# F9 
masterallF9 <- masterall[,c("Population","Line","CCRT_F","CSM_F","DT_A_F","DW_F","HSM_F","LS_F","SR_F","TL_F","WA_L_F")]
dataall2 <- na.omit(cbind(Country,masterallF9))

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

#####################
# PCA 1 = M9 (columns 1-3 are qualitative)

pca1<-PCA(dataall1,scale.unit=TRUE,graph=FALSE,quali.sup=c(1:3),axes=c(1,2))

# Here are the individual coordinates for PC1 output with the qualitative vars
PCA1_ind<-cbind(dataall1[,c(1:3)], pca1$ind$coord)
#
p1a<-fviz(pca1, title = "Male PCA - 9 traits (M9) PC1 vs PC2",
          element = "ind",
          habillage = as.factor(dataall1$Country),
          geom = c("point","text"),
          label= "none",
          pointsize = 2,
          pointshape = 18,
          palette = palette,
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          invisible="quali",
          ylim=c(3.4,-3.6),
          xlim=c(-4.3,4.8))

# labels for countries have to be added manually. I can do this later for this figure if required
p1aa<-p1a+
  scale_y_continuous(breaks = c(4, 2, 0, -2, -4))+
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4))+
  theme(text = element_text(size = 16), 
                plot.title = element_blank(),
                legend.position = "none",
                axis.title = element_text(size = 14),
                axis.text = element_text(size = 14),
        panel.grid.minor = element_blank())+
  annotate("text", x = 1.9, y = 0.7, label = "AT", cex = 5, colour = "#E38800") + 
  annotate("text", x = -0.5, y = -0.2, label = "DE", cex = 5, colour = "#F6C200") +
  annotate("text", x = 1.75, y = 0.1, label = "RU", cex = 5, colour = "#095888") +
  annotate("text", x = -1.05, y = -1.2, label = "FI", cex = 5, colour = "#A00E00") +
  annotate("text", x = 1.9, y = -0.65, label = "UA", cex = 5, colour = "#0086A8") +
  annotate("text", x = 1, y = -1.8, label = "DK", cex = 5, colour = "#D04E00") +
  annotate("text", x = 2.3, y = 2.2, label = "PT", cex = 5, colour = "#7BA354") +
  annotate("text", x = -2.4, y = 1.8, label = "TR", cex = 5, colour = "#132B69") +
  annotate("text", x = -0.8, y = 1.4, label = "ES", cex = 5, colour = "#B82E00") 

##############

M9_traitPC <- as.data.table(pca1$var$coord)
M9_traitPC[,trait:=rownames(pca1$var$coord)]

M9_traitPC$trait<-c("Chill-coma recovery time","Cold shock mortality",
                    "Egg-to-adult development time","Dry weight",
                    "Heat shock mortality","Life span",
                    "Starvation resistance","Thorax length","Wing area - Left")

M9_Dim1_sig <- subset(M9_traitPC, Dim.1 > 0.5 | Dim.1 < -0.5)
M9_Dim2_sig <- subset(M9_traitPC, Dim.2 > 0.5 | Dim.2 < -0.5)

M9_Dim1_sig[,ordDim1:=rank(Dim.1)]
M9_Dim2_sig[,ordDim2:=rank(Dim.2)]

xloading_theme <- theme(axis.text.y=element_blank(),
                      axis.ticks.y=element_blank(),
                      axis.text.x=element_text(size=12),
                      axis.title.x =element_text(size = 14),
                      axis.ticks.x=element_line(linewidth=0.5),
                      axis.ticks.length=unit(0.20,"cm"),
                      panel.grid = element_blank(),
                      panel.border = element_rect(colour = "grey70"))

yloading_theme <- theme(axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(),
                        axis.text.y=element_text(size=12),
                        axis.title.y =element_text(size = 14),
                        axis.ticks.y=element_line(linewidth=0.5),
                        axis.ticks.length=unit(0.20,"cm"),
                        panel.grid = element_blank(),
                        panel.border = element_rect(colour = "grey70"))

flip_M9PC1 <- 1
M9_traitPC1_load <- ggplot(data=M9_Dim1_sig) +
  geom_segment( aes(x=0, xend=flip_M9PC1*Dim.1, y=ordDim1, yend=ordDim1), arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
                linewidth =0.7, lineend = "round", linejoin = "mitre") +
  geom_text(data=M9_Dim1_sig[Dim.1< 0], aes(y=ordDim1, x=.05, label=trait), size=4.5, angle=0, hjust=0) +
  geom_text(data=M9_Dim1_sig[Dim.1> 0], aes(y=ordDim1, x=-.05, label=trait), size=4.5, angle=0, hjust=1) +
  xlim(-.85, .85) + ylim(0, 5.5) +
  xlab("Loading") + ylab("") +
  theme_bw() +
  xloading_theme


flip_M9PC2 <- 1
M9_traitPC2_load <- ggplot(data=M9_Dim2_sig) +
  geom_segment( aes(x=0, xend=flip_M9PC2*Dim.2, y=ordDim2, yend=ordDim2), arrow = arrow(length = unit(0.2, "cm"),type = "closed"),
                linewidth =0.7, lineend = "round", linejoin = "mitre") +
  geom_text(data=M9_Dim2_sig[Dim.2< 0], aes(y=ordDim2, x=.05, label=trait), size=4.5, angle=90, hjust=1) +
  geom_text(data=M9_Dim2_sig[Dim.2> 0], aes(y=ordDim2, x=-.05, label=trait), size=4.5, angle=90, hjust=0) +
  coord_flip() + xlim(.75, -.75) + ylim(4.5, 0) +
  xlab("Loading") + ylab("") +
  theme_bw() +
  yloading_theme

BLANK<-ggplot() + theme_void()

plot1<-plot_grid(M9_traitPC2_load,p1aa,BLANK,M9_traitPC1_load,
          rel_widths = c(1,4), rel_heights = c(3.5,1))


#####################
# PCA 2 = F9 (columns 1-3 are qualitative)

pca2<-PCA(dataall2,scale.unit=TRUE,graph=FALSE,quali.sup=c(1:3),axes=c(1,2))

# Here are the individual coordinates for PC2 output with the qualitative vars
PCA2_ind<-cbind(dataall2[,c(1:3)], pca2$ind$coord)

p2a<-fviz(pca2, title = "Female PCA - 9 traits (F9) PC1 vs PC2",
          element = "ind",
          habillage = as.factor(dataall2$Country),
          geom = c("point","text"),
          label= "none",
          pointsize = 2,
          pointshape = 18,
          palette = palette,
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          invisible="quali",
          xlim=c(-4.3,4.8),
          ylim=c(-3.4,3.6))

p2aa<-p2a+
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4))+
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4))+
  theme(text = element_text(size = 16), 
                plot.title = element_blank(),
                legend.position = "none",
                axis.title = element_text(size = 14),
                axis.text = element_text(size = 14),
                panel.grid.minor = element_blank())+
  annotate("text", x = 1.5, y = -0.1, label = "AT", cex = 5, colour = "#E38800") + 
  annotate("text", x = 1.2, y = 0.3, label = "DE", cex = 5, colour = "#F6C200") +
  annotate("text", x = 1.7, y = 1.0, label = "RU", cex = 5, colour = "#095888") +
  annotate("text", x = -1.7, y = 1.0, label = "FI", cex = 5, colour = "#A00E00") +
  annotate("text", x = 0.9, y = -1.5, label = "UA", cex = 5, colour = "#0086A8") +
  annotate("text", x = -0.3, y = 1.75, label = "DK", cex = 5, colour = "#D04E00") +
  annotate("text", x = 2.5, y = -1, label = "PT", cex = 5, colour = "#7BA354") +
  annotate("text", x = -1.5, y = -1.6, label = "TR", cex = 5, colour = "#132B69") +
  annotate("text", x = -1.7, y = -0.2, label = "ES", cex = 5, colour = "#B82E00") 
# labels can be added later

####

F9_traitPC <- as.data.table(pca2$var$coord)
F9_traitPC[,trait:=rownames(pca2$var$coord)]

F9_traitPC$trait<-c("Chill-coma recovery time","Cold shock mortality",
                    "Egg-to-adult development time","Dry weight",
                    "Heat shock mortality","Life span",
                    "Starvation resistance","Thorax length","Wing area - Left")


F9_Dim1_sig <- subset(F9_traitPC, Dim.1 > 0.5 | Dim.1 < -0.5)
F9_Dim2_sig <- subset(F9_traitPC, Dim.2 > 0.4 | Dim.2 < -0.4)

F9_Dim1_sig[,ordDim1:=rank(Dim.1)]
F9_Dim2_sig[,ordDim2:=rank(Dim.2)]


flip_F9PC1 <- 1
F9_traitPC1_load <- ggplot(data=F9_Dim1_sig) +
  geom_segment( aes(x=0, xend=flip_F9PC1*Dim.1, y=ordDim1, yend=ordDim1), arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
                linewidth =0.7, lineend = "round", linejoin = "mitre") +
  geom_text(data=F9_Dim1_sig[Dim.1<0], aes(y=ordDim1, x=.05, label=trait), size=4.5, angle=0, hjust=0) +
  geom_text(data=F9_Dim1_sig[Dim.1>0], aes(y=ordDim1, x=-.05, label=trait), size=4.5, angle=0, hjust=1) +
  xlim(-.85, .85) + ylim(-0.5, 5) +
  xlab("Loading") + ylab("") +
  theme_bw() +
  xloading_theme

flip_F9PC2 <- 1
F9_traitPC2_load <- ggplot(data=F9_Dim2_sig) +
  geom_segment( aes(x=0, xend=flip_F9PC2*Dim.2, y=ordDim2, yend=ordDim2), arrow = arrow(length = unit(0.2, "cm"),type = "closed"),
                linewidth =0.7, lineend = "round", linejoin = "mitre") +
  geom_text(data=F9_Dim2_sig[Dim.2<0], aes(y=ordDim2, x=.05, label=trait), size=4.5, angle=90, hjust=0) +
  geom_text(data=F9_Dim2_sig[Dim.2>0], aes(y=ordDim2, x=-.05, label=trait), size=4.5, angle=90, hjust=1) +
  coord_flip() + xlim(-.75, .75) + ylim(-0.5, 4) +
  xlab("Loading") + ylab("") +
  theme_bw() +
  yloading_theme


plot2<-plot_grid(F9_traitPC2_load,p2aa,BLANK,F9_traitPC1_load,
          rel_widths = c(1,4), rel_heights = c(3.5,1))

plot_fin<-plot_grid(BLANK,BLANK,plot1,plot2,
                    nrow = 2,
                    rel_heights = c(1,20),
                    labels = c('A', 'B','',''),
                    label_size = 18)

