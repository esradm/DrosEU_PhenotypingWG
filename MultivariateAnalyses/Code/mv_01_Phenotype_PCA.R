rm(list=ls(all=TRUE))
library(dplyr)
library(FactoMineR)
library(factoextra)
library(ggforce)
library(cowplot)

workingDir = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/MultivariateAnalyses/Data/";
setwd(workingDir); 
getwd();

master <- read.csv("all_models_compound_coefs_forMultivariateAnalysis.csv")

# Make variables for Country and the sexes
Country <- c(rep("Turkey", 20),rep("Portugal", 17),rep("Spain", 15),rep("Germany", 20),rep("Austria", 20),rep("Ukraine", 19),rep("Denmark", 20),rep("Russia", 20),rep("Finland", 22))
SexM <-c(rep("M", 173))
SexF <-c(rep("F", 173))

####
#### Formatting the data for the different versions of the PCA
####
#################################
# M9
masterM9 <- master[,c("Population","Line","CCRT_M","CSM_M","DT_A_M","DW_M","HSM_M","LS_M","SR_M","TL_M","WA_L_M")]
data1 <- na.omit(cbind(Country,masterM9))
# F9 
masterF9 <- master[,c("Population","Line","CCRT_F","CSM_F","DT_A_F","DW_F","HSM_F","LS_F","SR_F","TL_F","WA_L_F")]
data2 <- na.omit(cbind(Country,masterF9))
# Fmax (not included in final analyses)
# masterFmax <- master[,c("Population","Line","CCRT_F","CSM_F","DT_A_F","Dia_F","DW_F","Fec_F","HSM_F","LS_F","Pgm_Total_F","SR_F","TL_F","WA_L_F")]
# data3 <- na.omit(cbind(Country,masterFmax))
# Fmaxplus
masterFmaxP <- master[,c("Population","Line","CCRT_F","CSM_F","DT_A_F","Dia_F","DW_F","Fec_F","HSM_F","LS_F","Pgm_Total_F","SR_F","TL_F","WA_L_F", "Via_NA")]
data4 <- na.omit(cbind(Country,masterFmaxP))

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
pca1<-PCA(data1,scale.unit=TRUE,graph=FALSE,quali.sup=c(1:3),axes=c(1,2))
pca1$eig
pca1$var$coord
fviz_contrib(pca1, choice = "var", axes = c(1,3), top = 10)
fviz_contrib(pca1, choice = "var", axes = 1, top = 10)
fviz_contrib(pca1, choice = "var", axes = 2, top = 10)
fviz_contrib(pca1, choice = "var", axes = 3, top = 10)

# Here are the individual coordinates for PC1 output with the qualitative vars
PCA1_ind<-cbind(data1[,c(1:3)], pca1$ind$coord)
p1<-plot(pca1,choix="var",axes=c(1,2),cex=1.4)
q1<-plot(pca1,choix="var",axes=c(2,3),cex=1.4)
# simple ellipses with chosen colours
plotellipses(pca1, keepvar = c(1), axes = c(1, 2),label = "quali", level = 0.95, palette=palette)
plotellipses(pca1, keepvar = c(1), axes = c(1, 3),label = "quali", level = 0.95, palette=palette)
plotellipses(pca1, keepvar = c(1), axes = c(2, 3),label = "quali", level = 0.95, palette=palette)

p1a<-fviz(pca1, title = "Male PCA - 9 traits (M9) PC1 vs PC2",
            element = "ind",
                     habillage =  as.factor(data1$Country),
                     geom = c("point","text"),
                     label= "quali",
                     pointsize = 2,
                     pointshape = 18,
                     palette = palette,
                     addEllipses = TRUE, # Concentration ellipses
                     ellipse.type="confidence",
                     legend.title = "Treatment",invisible="quali")

# labels for countries have to be added manually.
p1aa<-p1a+theme(text = element_text(size = 16), 
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14))+
  annotate("text", x = 2.5, y = 0.7, label = "Austria", cex = 5, colour = "#E38800") + 
  annotate("text", x = -1.3, y = -0.2, label = "Germany", cex = 5, colour = "#F6C200") +
  annotate("text", x = 2.4, y = 0.3, label = "Russia", cex = 5, colour = "#095888") +
  annotate("text", x = -1.7, y = -1.2, label = "Finland", cex = 5, colour = "#A00E00") +
  annotate("text", x = 2.5, y = -0.5, label = "Ukraine", cex = 5, colour = "#0086A8") +
  annotate("text", x = 1.7, y = -1.8, label = "Denmark", cex = 5, colour = "#D04E00") +
  annotate("text", x = 3.0, y = 2.2, label = "Portugal", cex = 5, colour = "#7BA354") +
  annotate("text", x = -2.2, y = 1.8, label = "Turkey", cex = 5, colour = "#132B69") +
  annotate("text", x = -0.5, y = 1.4, label = "Spain", cex = 5, colour = "#B82E00") 

plot_grid(p1,p1aa)

#####################
# PCA 2 = F9 (columns 1-3 are qualitative)
pca2<-PCA(data2,scale.unit=TRUE,graph=FALSE,quali.sup=c(1:3),axes=c(1,2))
pca2$eig
pca2$var$coord
# Here are the individual coordinates for PC2 output with the qualitative vars
PCA2_ind<-cbind(data2[,c(1:3)], pca2$ind$coord)
p2<-plot(pca2,choix="var",axes=c(1,2),cex=1.4)
q2<-plot(pca2,choix="var",axes=c(1,3),cex=1.4)
# simple ellipses with chosen colours
plotellipses(pca2, keepvar = c(1), axes = c(1, 2),label = "quali", level = 0.95, palette=palette)
plotellipses(pca2, keepvar = c(1), axes = c(1, 3),label = "quali", level = 0.95, palette=palette)
p2a<-fviz(pca2, title = "Female PCA - 9 traits (F9) PC1 vs PC2",
          element = "ind",
          habillage =  as.factor(data2$Country),
          geom = c("point","text"),
          label= "quali",
          pointsize = 2,
          pointshape = 18,
          palette = palette,
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          legend.title = "Treatment",invisible="quali")

p2aa<-p2a+theme(text = element_text(size = 16), 
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14))+
  annotate("text", x = 1.6, y = 0.1, label = "Austria", cex = 5, colour = "#E38800") + 
  annotate("text", x = 1.7, y = 0.4, label = "Germany", cex = 5, colour = "#F6C200") +
  annotate("text", x = 1.9, y = 1.0, label = "Russia", cex = 5, colour = "#095888") +
  annotate("text", x = -2.1, y = 1.0, label = "Finland", cex = 5, colour = "#A00E00") +
  annotate("text", x = 0.9, y = -1.5, label = "Ukraine", cex = 5, colour = "#0086A8") +
  annotate("text", x = 0.1, y = 1.8, label = "Denmark", cex = 5, colour = "#D04E00") +
  annotate("text", x = 2.8, y = -1.2, label = "Portugal", cex = 5, colour = "#7BA354") +
  annotate("text", x = -1.8, y = -1.6, label = "Turkey", cex = 5, colour = "#132B69") +
  annotate("text", x = -1.9, y = -0.2, label = "Spain", cex = 5, colour = "#B82E00") 

p2b<-fviz(pca2, title = "Female PCA - 9 traits (F9) PC1 vs PC3",
          element = "ind", axes = c(1, 3),
          habillage =  as.factor(data2$Country),
          geom = c("point","text"),
          label= "quali",
          pointsize = 2,
          pointshape = 18,
          palette = palette,
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          legend.title = "Treatment",invisible="quali")

p2bb<-p2b+theme(text = element_text(size = 16), 
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14))+
  annotate("text", x = 1.7, y = 0.7, label = "Austria", cex = 5, colour = "#E38800") + 
  annotate("text", x = -0.4, y = 0.9, label = "Ger.", cex = 5, colour = "#F6C200") +
  annotate("text", x = 0.8, y = -0.6, label = "Rus.", cex = 5, colour = "#095888") +
  annotate("text", x = -1.1, y = 1.3, label = "Finland", cex = 5, colour = "#A00E00") +
  annotate("text", x = 0.9, y = 2.3, label = "Ukraine", cex = 5, colour = "#0086A8") +
  annotate("text", x = 0.9, y = -1.3, label = "Denmark", cex = 5, colour = "#D04E00") +
  annotate("text", x = 3.1, y = -0.5, label = "Portugal", cex = 5, colour = "#7BA354") +
  annotate("text", x = -2.0, y = -0.2, label = "Turkey", cex = 5, colour = "#132B69") +
  annotate("text", x = -2.0, y = -0.8, label = "Spain", cex = 5, colour = "#B82E00") 

plot_grid(p2,p2aa, q2, p2bb)

#####################
# # PCA 3 = Fmax  (columns 1-3 are qualitative)
# pca3<-PCA(data3,scale.unit=TRUE,graph=FALSE,quali.sup=c(1:3),axes=c(1,2))
# pca3$eig
# pca3$var$coord
# # Here are the individual coordinates for PC2 output with the qualitative vars
# PCA3_ind<-cbind(data3[,c(1:3)], pca3$ind$coord)
# p3<-plot(pca3,choix="var",axes=c(1,2),cex=1.4)
# q3<-plot(pca3,choix="var",axes=c(1,3),cex=1.4)
# # simple ellipses with chosen colours
# plotellipses(pca3, keepvar = c(1), axes = c(1, 2),label = "quali", level = 0.95, palette=palette)
# plotellipses(pca3, keepvar = c(1), axes = c(1, 2),label = "quali", level = 0.95, palette=palette)
# p3a<-fviz(pca3, title = "Female PCA - 12 traits (Fmax) PC1 vs PC2",
#           element = "ind",
#           habillage =  as.factor(data3$Country),
#           geom = c("point","text"),
#           label= "quali",
#           pointsize = 2,
#           pointshape = 18,
#           palette = palette,
#           addEllipses = TRUE, # Concentration ellipses
#           ellipse.type="confidence",
#           legend.title = "Treatment",invisible="quali")
# 
# p3aa<-p3a+theme(text = element_text(size = 16), 
#           axis.title = element_text(size = 14),
#           axis.text = element_text(size = 14))+
#   annotate("text", x = 1.8, y = 0.3, label = "Austria", cex = 5, colour = "#E38800") + 
#   annotate("text", x = -0.8, y = 0.2, label = "Ger.", cex = 5, colour = "#F6C200") +
#   annotate("text", x = 1.4, y = 1.3, label = "Russia", cex = 5, colour = "#095888") +
#   annotate("text", x = -2.1, y = 1.2, label = "Finland", cex = 5, colour = "#A00E00") +
#   annotate("text", x = 0.4, y = -1.3, label = "Ukr.", cex = 5, colour = "#0086A8") +
#   annotate("text", x = -0.2, y = 1.9, label = "Denmark", cex = 5, colour = "#D04E00") +
#   annotate("text", x = 3.1, y = -0.5, label = "Portugal", cex = 5, colour = "#7BA354") +
#   annotate("text", x = -1.8, y = -1.8, label = "Turkey", cex = 5, colour = "#132B69") +
#   annotate("text", x = -2.0, y = -0.2, label = "Spain", cex = 5, colour = "#B82E00") 
# 
# p3b<-fviz(pca3, title = "Female PCA - 12 traits (Fmax) PC1 vs PC3",
#           element = "ind", axes = c(1, 3),
#           habillage =  as.factor(data3$Country),
#           geom = c("point","text"),
#           label= "quali",
#           pointsize = 2,
#           pointshape = 18,
#           palette = palette,
#           addEllipses = TRUE, # Concentration ellipses
#           ellipse.type="confidence",
#           legend.title = "Treatment",invisible="quali")
# 
# p3bb<-p3b+theme(text = element_text(size = 16), 
#           axis.title = element_text(size = 14),
#           axis.text = element_text(size = 14))+
#   annotate("text", x = 2.0, y = 0.7, label = "Austria", cex = 5, colour = "#E38800") + 
#   annotate("text", x = -0.6, y =1.0, label = "Germany", cex = 5, colour = "#F6C200") +
#   annotate("text", x = 0.8, y = -1.2, label = "Rus.", cex = 5, colour = "#095888") +
#   annotate("text", x = -1.8, y = 0.7, label = "Finland", cex = 5, colour = "#A00E00") +
#   annotate("text", x = 0.9, y = 2.7, label = "Ukraine", cex = 5, colour = "#0086A8") +
#   annotate("text", x = 0.1, y = -1.2, label = "Den.", cex = 5, colour = "#D04E00") +
#   annotate("text", x = 3.1, y = -0.5, label = "Portugal", cex = 5, colour = "#7BA354") +
#   annotate("text", x = -2.2, y = -0.4, label = "Turkey", cex = 5, colour = "#132B69") +
#   annotate("text", x = -1.4, y = -2.0, label = "Spain", cex = 5, colour = "#B82E00") 

plot_grid(p3,p3aa, q3, p3bb)

#####################
# PCA = Fmaxplus (columns 1-3 are qualitative)
pca4<-PCA(data4,scale.unit=TRUE,graph=FALSE,quali.sup=c(1:3),axes=c(1,2))
pca4$eig
pca4$var$coord
# Here are the individual coordinates for PC2 output with the qualitative vars
PCA4_ind<-cbind(data4[,c(1:3)], pca4$ind$coord)
p4<-plot(pca4,choix="var",axes=c(1,2),cex=1.4)
q4<-plot(pca4,choix="var",axes=c(1,3),cex=1.4)
# simple ellipses with chosen colours
plotellipses(pca4, keepvar = c(1), axes = c(1, 2),label = "quali", level = 0.95, palette=palette)
plotellipses(pca4, keepvar = c(1), axes = c(1, 3),label = "quali", level = 0.95, palette=palette)
p4a<-fviz(pca4, title = "Female PCA - 13 traits (Fmax Plus) PC1 vs PC2",
          element = "ind",
          habillage =  as.factor(data4$Country),
          geom = c("point","text"),
          label= "quali",
          pointsize = 2,
          pointshape = 18,
          palette = palette,
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          legend.title = "Treatment",invisible="quali")

p4aa<-p4a+theme(text = element_text(size = 16), 
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14))+
  annotate("text", x = 2.2, y = 0.3, label = "Austria", cex = 5, colour = "#E38800") + 
  annotate("text", x = -0.8, y = 0.2, label = "Ger.", cex = 5, colour = "#F6C200") +
  annotate("text", x = -0.2, y = -0.2, label = "Russia", cex = 5, colour = "#095888") +
  annotate("text", x = -2.4, y = 1.2, label = "Finland", cex = 5, colour = "#A00E00") +
  annotate("text", x = 1.4, y = 1.4, label = "Ukraine", cex = 5, colour = "#0086A8") +
  annotate("text", x = -0.2, y = 1.9, label = "Denmark", cex = 5, colour = "#D04E00") +
  annotate("text", x = 3.3, y = -0.5, label = "Portugal", cex = 5, colour = "#7BA354") +
  annotate("text", x = -1.8, y = -1.8, label = "Turkey", cex = 5, colour = "#132B69") +
  annotate("text", x = -1.9, y = -0.5, label = "Spain", cex = 5, colour = "#B82E00") 

p4b<-fviz(pca4, title = "Female PCA - 13 traits (Fmax Plus) PC1 vs PC3",
          element = "ind", axes = c(1,3),
          habillage =  as.factor(data4$Country),
          geom = c("point","text"),
          label= "quali",
          pointsize = 2,
          pointshape = 18,
          palette = palette,
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          legend.title = "Treatment",invisible="quali")

p4bb<-p4b+theme(text = element_text(size = 16), 
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14))+
  annotate("text", x = 2.0, y = 0.7, label = "Austria", cex = 5, colour = "#E38800") + 
  annotate("text", x = -0.6, y =1.0, label = "Germany", cex = 5, colour = "#F6C200") +
  annotate("text", x = 1.7, y = -1.0, label = "Russia", cex = 5, colour = "#095888") +
  annotate("text", x = -2.4, y = -0.4, label = "Finland", cex = 5, colour = "#A00E00") +
  annotate("text", x = 0.9, y = 2.7, label = "Ukraine", cex = 5, colour = "#0086A8") +
  annotate("text", x = 0.6, y = -1.5, label = "Denmark", cex = 5, colour = "#D04E00") +
  annotate("text", x = 3.1, y = -0.5, label = "Portugal", cex = 5, colour = "#7BA354") +
  annotate("text", x = -1.8, y = 0.7, label = "Turkey", cex = 5, colour = "#132B69") +
  annotate("text", x = -0.9, y = -1.8, label = "Spain", cex = 5, colour = "#B82E00") 

plot_grid(p4,p4aa, q4, p4bb)

#######
##
## Final plots for manuscript comparing F9 and M9 with loadings plotted on x and y axes
## 
## Note that this reuses the p1a (M9) and p2a (F9) objects from before

# Custom ggplot2 themes are used to plot the loadings
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

# labels for countries added manually.
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

# extract loadings
M9_traitPC <- as.data.table(pca1$var$coord)
# order phenotypes alphabetically
M9_traitPC[,trait:=rownames(pca1$var$coord)]
# use longer version of phenotpye name
M9_traitPC$trait<-c("Chill-coma recovery time","Cold shock mortality",
                    "Egg-to-adult development time","Dry weight",
                    "Heat shock mortality","Life span",
                    "Starvation resistance","Thorax length","Wing area - Left")

# only use traits with loadings above a certain value
M9_Dim1_sig <- subset(M9_traitPC, Dim.1 > 0.5 | Dim.1 < -0.5)
M9_Dim2_sig <- subset(M9_traitPC, Dim.2 > 0.5 | Dim.2 < -0.5)
# order these loadings
M9_Dim1_sig[,ordDim1:=rank(Dim.1)]
M9_Dim2_sig[,ordDim2:=rank(Dim.2)]

# Loadings for PC1
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

# Loadings for PC2
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

# blank plot used as a spacer
BLANK<-ggplot() + theme_void()

# composite plot
plot1<-plot_grid(M9_traitPC2_load,p1aa,BLANK,M9_traitPC1_load,
                 rel_widths = c(1,4), rel_heights = c(3.5,1))

#############
# PCA 2 = F9 

# reusing th p2a object created before, this time with simpler labels
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

# Extract loadings
F9_traitPC <- as.data.table(pca2$var$coord)
# Order by trait name (alphbetical)
F9_traitPC[,trait:=rownames(pca2$var$coord)]
# Use longer version of phenotype name
F9_traitPC$trait<-c("Chill-coma recovery time","Cold shock mortality",
                    "Egg-to-adult development time","Dry weight",
                    "Heat shock mortality","Life span",
                    "Starvation resistance","Thorax length","Wing area - Left")

# only use traits with loadings above a certain value
F9_Dim1_sig <- subset(F9_traitPC, Dim.1 > 0.5 | Dim.1 < -0.5)
F9_Dim2_sig <- subset(F9_traitPC, Dim.2 > 0.4 | Dim.2 < -0.4)

# Order these loadings
F9_Dim1_sig[,ordDim1:=rank(Dim.1)]
F9_Dim2_sig[,ordDim2:=rank(Dim.2)]

# Loadings for PC1
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

# Loadings for PC2
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

# composite plot F9 data
plot2<-plot_grid(F9_traitPC2_load,p2aa,BLANK,F9_traitPC1_load,
                 rel_widths = c(1,4), rel_heights = c(3.5,1))

##############################
# Final composite plot for MS

plot_fin<-plot_grid(BLANK,BLANK,plot1,plot2,
                    nrow = 2,
                    rel_heights = c(1,20),
                    labels = c('A', 'B','',''),
                    label_size = 18)

# Save PCA coordinates in csv files for scripts mv_03A and mv_03B
write.csv(PCA1_ind, file = "M9_drosEU_PCcoords.csv", row.names = F)
write.csv(PCA2_ind, file = "F9_drosEU_PCcoords.csv", row.names = F)
#write.csv(PCA3_ind, file = "Fmax_drosEU_PCcoords.csv", row.names = F)
write.csv(PCA4_ind, file = "FmaxP_drosEU_PCcoords.csv", row.names = F)

# Save PCA objects for use in script mv_03D
save(pca1, file = "M9_drosEU.RData")
save(pca2, file = "F9_drosEU.RData")
#save(pca3, file = "Fmax_drosEU.RData")
save(pca4, file = "FmaxP_drosEU.RData")

####### OPTIONAL: Save PCA objects as R data and/or write individual coordinates to csv output
# save.image(file = "PCA_results_drosEU.RData")
