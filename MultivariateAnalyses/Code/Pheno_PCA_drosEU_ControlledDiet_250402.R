rm(list=ls(all=TRUE))
library(dplyr)
library(FactoMineR)
library(factoextra)
library(ggforce)
library(cowplot)
library(tidyverse)
library(psych)
library(corrplot)

# Load the two data sets into R. They are structured rather differently, so need to be treated a little differently prior to running the PCAs

masterall <- read.csv("C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/all_models_line_meta_compound_random_coefs_wide_220916.csv")
masterctl <- read.csv("C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/all_models_line_random_coefs_similar_diet_v2.csv")

library("Hmisc")

# What are the traits in the data (this is the case for both data sets)
unique(masterctl$Trait)
# Levels: CCRT CSM Dia DT_A DW Fec HSM LS Pgm_Total SR TL Via WA_L

# Analysing the sexes separately. For males this is easy, but there are a few different ways of doing this for females:
# These will be the PCAs analysed:
# M9.       9 male traits:
#           "CCRT_M","CSM_M","DT_A_M","DW_M","HSM_M","LS_M","SR_M","TL_M","WA_L_M"
# F9.       9 female traits corresponding to these male traits:
#           "CCRT_F","CSM_F","DT_A_F","DW_F","HSM_F","LS_F","SR_F","TL_F","WA_L_F"
# Fmax.     All 12 possible female traits:
#           "CCRT_F","CSM_F","DT_A_F","Dia_F","DW_F","Fec_F","HSM_F","LS_F","Pgm_Total_F","SR_F","TL_F","WA_L_F"
# Fmaxplus. All 12 possible female traits plus viability:
#           "CCRT_F","CSM_F","DT_A_F","Dia_F","DW_F","Fec_F","HSM_F","LS_F","Pgm_Total_F","SR_F","TL_F","WA_L_F", "Via_NA"

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

# Female 13 traits
masterallFmaxP <- masterall[,c("Population","Line","CCRT_F","CSM_F","Dia_F","DT_A_F","DW_F","Fec_F","HSM_F","LS_F","Pgm_Total_F","SR_F","TL_F","WA_L_F", "Via_NA")]
dataall_F13 <- na.omit(cbind(Country,masterallFmaxP))

# Female 9 traits 
masterallF9 <- masterall[,c("Population","Line","CCRT_F","CSM_F","DT_A_F","DW_F","HSM_F","LS_F","SR_F","TL_F","WA_L_F")]
dataall_F9 <- na.omit(cbind(Country,masterallF9))

# Male 9 traits
masterallM9 <- masterall[,c("Population","Line","CCRT_M","CSM_M","DT_A_M","DW_M","HSM_M","LS_M","SR_M","TL_M","WA_L_M")]
dataall_M9 <- na.omit(cbind(Country,masterallM9))

# Male 8 traits -  for the purposes of comparing with the P:C diet controlled data set, which is missing the TL trait
masterallM8 <- masterall[,c("Population","Line","CCRT_M","CSM_M","DT_A_M","DW_M","HSM_M","LS_M","SR_M","WA_L_M")]
dataall_M8 <- na.omit(cbind(Country,masterallM8))

###################

Corr_FmaxPP <- masterall[,c("CCRT_F","CSM_F","DT_A_F","Dia_F","DT_P_NA","DW_F","Fec_F","HSM_F","LS_F","Pgm_Total_F","SR_F","TL_F","WA_L_F", "WA_R_F", "Via_NA")]

res2 <- rcorr(as.matrix(Corr_FmaxPP))
res2$r
res2$P

##################################
###
### Preparing the masterctl data
###

masterOrd<-masterctl[order(masterctl$Population),]

table(masterOrd$Population)
# for v2
# AK  GI  KA  MA  MU  RE  UM  VA  YE 
# 406 290 388 387 393 323 348 394 393 

masterOrd$Country <- as.factor(c(rep("Finland", 406),rep("Spain", 290),rep("Denmark", 388),rep("Austria", 387),rep("Germany", 393),rep("Portugal", 323),rep("Ukraine", 348),rep("Russia", 394),rep("Turkey", 393)))

# Table showing which labs have been selected for which traits
table(masterOrd$Lab, masterOrd$Trait)

# get separate male and female data
masterf <- subset(masterOrd, Sex == "F")
masterm <- subset(masterOrd, Sex == "M")

#######################################################
### For Females
# Create small data frames of the coefficient values
dfCCRT_F <- subset(masterf, Trait=="CCRT")[,c(9,5:7)]
dfCSM_F <- subset(masterf, Trait=="CSM")[,c(9,5:7)]
dfDia_F <- subset(masterf, Trait=="Dia")[,c(9,5:7)]
dfDT_A_F <- subset(masterf, Trait=="DT_A")[,c(9,5:7)]
dfDW_F <- subset(masterf, Trait=="DW")[,c(9,5:7)]
dfFec_F <- subset(masterf, Trait=="Fec")[,c(9,5:7)]
dfHSM_F <- subset(masterf, Trait=="HSM")[,c(9,5:7)]
dfLS_F <- subset(masterf, Trait=="LS")[,c(9,5:7)]
dfPgm_Total_F <- subset(masterf, Trait=="Pgm_Total")[,c(9,5:7)]
dfSR_F <- subset(masterf, Trait=="SR")[,c(9,5:7)]
dfTL_F <- subset(masterf, Trait=="TL")[,c(9,5:7)]
dfVia_NA <- subset(masterOrd, Trait=="Via")[,c(9,5:7)]
dfWA_L_F <- subset(masterf, Trait=="WA_L")[,c(9,5:7)]

# put all data frames into list
F_df_list <- list(dfCCRT_F, dfCSM_F, dfDia_F, dfDT_A_F, dfDW_F, dfFec_F, dfHSM_F, dfLS_F, dfPgm_Total_F, dfSR_F, dfTL_F, dfWA_L_F, dfVia_NA)
# merge all list by the qualitative variables
Fdata<- na.omit(F_df_list %>% reduce(full_join, by=c('Country','Population','Line')))
# rename
colnames(Fdata)<-c('Country','Population', 'Line', 'CCRT_F', 'CSM_F', 'Dia_F', 'DT_A_F', 'DW_F', 'Fec_F', 'HSM_F', 'LS_F', 'Pgm_Total_F', 'SR_F', 'TL_F', 'WA_L_F', 'Via_NA')

#######################################################
### For Males
# Create small data frames of the coefficient values
dfCCRT_M <- subset(masterm, Trait=="CCRT")[,c(9,5:7)]
dfCSM_M <- subset(masterm, Trait=="CSM")[,c(9,5:7)]
dfDT_A_M <- subset(masterm, Trait=="DT_A")[,c(9,5:7)]
dfDW_M <- subset(masterm, Trait=="DW")[,c(9,5:7)]
dfHSM_M <- subset(masterm, Trait=="HSM")[,c(9,5:7)]
dfLS_M <- subset(masterm, Trait=="LS")[,c(9,5:7)]
dfSR_M <- subset(masterm, Trait=="SR")[,c(9,5:7)]
# dfTL_M <- subset(masterm, Trait=="TL")[,c(9,5:7)]
# There is no TL data for males...
dfWA_L_M <- subset(masterm, Trait=="WA_L")[,c(9,5:7)]

# put all data frames into list
M_df_list <- list(dfCCRT_M, dfCSM_M, dfDT_A_M, dfDW_M, dfHSM_M, dfLS_M, dfSR_M, dfWA_L_M)
# merge all list by the qualitative variables
Mdata<- na.omit(M_df_list %>% reduce(full_join, by=c('Country','Population','Line')))
# rename
colnames(Mdata)<-c('Country','Population', 'Line', 'CCRT_M', 'CSM_M', 'DT_A_M', 'DW_M', 'HSM_M', 'LS_M', 'SR_M', 'WA_L_M')

#######################################################
#### Formatting the data for the different versions of the controlled diet PCA
####
# Fmaxplus
datactl_F13 <- Fdata
# F9 
datactl_F9 <- Fdata[,c("Country","Population","Line","CCRT_F","CSM_F","DT_A_F","DW_F","HSM_F","LS_F","SR_F","TL_F","WA_L_F")]
# M8
datactl_M8 <- Mdata


#####
###
### Save the various data files ready for PCA

### All data
saveRDS(dataall_F13, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/pca_dataall_F13.rds")
saveRDS(dataall_F9, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/pca_dataall_F9.rds")
saveRDS(dataall_M9, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/pca_dataall_M9.rds")
saveRDS(dataall_M8, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/pca_dataall_M8.rds")

saveRDS(datactl_F13, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/pca_datactl_F13.rds")
saveRDS(datactl_F9, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/pca_datactl_F9.rds")
saveRDS(datactl_M8, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/pca_datactl_M8.rds")


#####################
# Colour palette
#    Country  Color
#   Austria #E38800
#   Germany #F6C200
#    Russia #095888
#   Finland #A00E00
#   Ukraine #0086A8
#   Denmark #D04E00
#  Portugal #7BA354
#    Turkey #132B69
#     Spain #B82E00

palette =c("#E38800","#D04E00","#A00E00","#F6C200","#7BA354","#095888","#B82E00","#132B69","#0086A8")


###########################
##
## Female 13 traits PCA
##
###########################

pca_F13<-PCA(dataall_F13,scale.unit=TRUE,graph=FALSE,quali.sup=c(1:3),axes=c(1,2))
pcactl_F13<-PCA(datactl_F13,scale.unit=TRUE,graph=FALSE,quali.sup=c(1:3),axes=c(1,2))
colnames(pcactl_F13$var$coord) <- c("ctl_dim1", "ctl_dim2", "ctl_dim3", "ctl_dim4", "ctl_dim5")  
# Tucker's congruence between PCs of the all data vs control data
factcor_F13<-factor.congruence(pca_F13$var$coord, pcactl_F13$var$coord,2)
corrplot(factcor_F13, method="number")

pca_F13$eig
pca_F13$var$coord
fviz_eig(pca_F13)

# Here are the individual coordinates for the qualitative vars
pca_F13_ind<-cbind(dataall_F13[,c(1:3)], pca_F13$ind$coord)
p_F13_v1<-plot(pca_F13,choix="var",axes=c(1,2),cex=1.4) + 
  theme(panel.grid = element_blank())
p_F13_v2<-plot(pca_F13,choix="var",axes=c(3,4),cex=1.4) + 
  theme(panel.grid = element_blank())

# plot ellipses with chosen colours
F13_c1<-fviz(pca_F13, title = "F13 PC1 vs PC2",
          element = "ind",
          habillage =  as.factor(dataall_F13$Country),
          geom = c("point","text"),
          label= "quali",
          pointsize = 2,
          pointshape = 18,
          palette = palette,
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          legend.title = "Treatment",invisible="quali")

p_F13_c1<-F13_c1+theme(text = element_text(size = 16), 
                axis.title = element_text(size = 14),
                axis.text = element_text(size = 14))+
  annotate("text", x = 0.6, y = -0.7, label = "AT", cex = 5, colour = "#E38800") + 
  annotate("text", x = -0.8, y = 0.2, label = "DE", cex = 5, colour = "#F6C200") +
  annotate("text", x = -0.2, y = -0.2, label = "RU", cex = 5, colour = "#095888") +
  annotate("text", x = -2.1, y = 0.8, label = "FI", cex = 5, colour = "#A00E00") +
  annotate("text", x = 1.4, y = 1.4, label = "UA", cex = 5, colour = "#0086A8") +
  annotate("text", x = -0.5, y = 1.9, label = "DK", cex = 5, colour = "#D04E00") +
  annotate("text", x = 2.8, y = -0.5, label = "PT", cex = 5, colour = "#7BA354") +
  annotate("text", x = -1.8, y = -1.8, label = "TR", cex = 5, colour = "#132B69") +
  annotate("text", x = -1.55, y = -0.4, label = "ES", cex = 5, colour = "#B82E00") +
  theme(legend.position = "none")

F13_c2<-fviz(pca_F13, title = "F13 PC3 vs PC4",
          element = "ind", axes = c(3,4),
          habillage =  as.factor(dataall_F13$Country),
          geom = c("point","text"),
          label= "quali",
          pointsize = 2,
          pointshape = 18,
          palette = palette,
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          legend.title = "Treatment",invisible="quali")

p_F13_c2<-F13_c2+theme(text = element_text(size = 16), 
                axis.title = element_text(size = 14),
                axis.text = element_text(size = 14))+
  annotate("text", x = 1.3, y = -1.5, label = "AT", cex = 5, colour = "#E38800") + 
  annotate("text", x = 0.25, y = 1.2, label = "DE", cex = 5, colour = "#F6C200") +
  annotate("text", x = -1.3, y = 1.15, label = "RU", cex = 5, colour = "#095888") +
  annotate("text", x = -0.4, y = 1.2, label = "FI", cex = 5, colour = "#A00E00") +
  annotate("text", x = 2.7, y = 0.7, label = "UA", cex = 5, colour = "#0086A8") +
  annotate("text", x = -0.75, y = -1.1, label = "DK", cex = 5, colour = "#D04E00") +
  annotate("text", x = -0.25, y = -1.2, label = "PT", cex = 5, colour = "#7BA354") +
  annotate("text", x = 1.1, y = 0.0, label = "TR", cex = 5, colour = "#132B69") +
  annotate("text", x = -1.9, y = -0.3, label = "ES", cex = 5, colour = "#B82E00") +
  theme(legend.position = "none")

plot_grid(p_F13_v1,p_F13_c1)
plot_grid(p_F13_v2,p_F13_c2)


###########################
##
## Female 9 traits PCA
##
###########################

pca_F9<-PCA(dataall_F9,scale.unit=TRUE,graph=FALSE,quali.sup=c(1:3),axes=c(1,2))
pcactl_F9<-PCA(datactl_F9,scale.unit=TRUE,graph=FALSE,quali.sup=c(1:3),axes=c(1,2))
fviz_contrib(pcactl_F9, choice = "var", axes = 1, top = 10)
fviz_contrib(pcactl_F9, choice = "var", axes = 2, top = 10)
fviz_contrib(pcactl_F9, choice = "var", axes = 3, top = 10)
colnames(pcactl_F9$var$coord) <- c("ctl_dim1", "ctl_dim2", "ctl_dim3", "ctl_dim4", "ctl_dim5")  

# Tucker's congruence between PCs of the all data vs control data
factcor_F9<-factor.congruence(pca_F9$var$coord, pcactl_F9$var$coord,2)
corrplot(factcor_F9, method="number")
# Good congruence of PC1 and moderate congruence of PC3
# Note that PC2 has negative congruence of -0.28, but can be flipped so that it has a
# positive (though still low) congruence of 0.28


#
pca_F9$eig
pca_F9$var$coord
fviz_eig(pca_F9)


# Here are the individual coordinates for PC2 output with the qualitative vars
pca_F9_ind<-cbind(dataall_F9[,c(1:3)], pca_F9$ind$coord)
p_F9_v1<-plot(pca_F9,choix="var",axes=c(1,2),cex=1.4) + 
  theme(panel.grid = element_blank())
p_F9_v2<-plot(pca_F9,choix="var",axes=c(3,4),cex=1.4) + 
  theme(panel.grid = element_blank())
# simple ellipses with chosen colours

F9_c1<-fviz(pca_F9, title = "F9 PC1 vs PC2",
          element = "ind",
          habillage =  as.factor(dataall_F9$Country),
          geom = c("point","text"),
          label= "quali",
          pointsize = 2,
          pointshape = 18,
          palette = palette,
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          legend.title = "Treatment",invisible="quali")

p_F9_c1<-F9_c1+theme(text = element_text(size = 16), 
                axis.title = element_text(size = 14),
                axis.text = element_text(size = 14))+
  annotate("text", x = 1.6, y = 0.1, label = "AT", cex = 5, colour = "#E38800") + 
  annotate("text", x = 1.4, y = 0.4, label = "DE", cex = 5, colour = "#F6C200") +
  annotate("text", x = 1.9, y = 1.0, label = "RU", cex = 5, colour = "#095888") +
  annotate("text", x = -2.1, y = 1.0, label = "FI", cex = 5, colour = "#A00E00") +
  annotate("text", x = 0.9, y = -1.5, label = "UA", cex = 5, colour = "#0086A8") +
  annotate("text", x = 0.1, y = 1.8, label = "DK", cex = 5, colour = "#D04E00") +
  annotate("text", x = 2.8, y = -1.2, label = "PT", cex = 5, colour = "#7BA354") +
  annotate("text", x = -1.8, y = -1.6, label = "TR", cex = 5, colour = "#132B69") +
  annotate("text", x = -1.9, y = -0.2, label = "ES", cex = 5, colour = "#B82E00") +
  theme(legend.position = "none")
# labels can be added later

F9_c2<-fviz(pca_F9, title = "F9 PC3 vs PC4",
          element = "ind", axes = c(3, 4),
          habillage =  as.factor(dataall_F9$Country),
          geom = c("point","text"),
          label= "quali",
          pointsize = 2,
          pointshape = 18,
          palette = palette,
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          legend.title = "Treatment",invisible="quali")

p_F9_c2<-F9_c2+theme(text = element_text(size = 16), 
                axis.title = element_text(size = 14),
                axis.text = element_text(size = 14))+
  annotate("text", x = 0.2, y = -1.65, label = "AT", cex = 5, colour = "#E38800") + 
  annotate("text", x = 0.2, y = 1.2, label = "DE", cex = 5, colour = "#F6C200") +
  annotate("text", x = 0.45, y = -0.25, label = "RU", cex = 5, colour = "#095888") +
  annotate("text", x = 1.1, y = 1.4, label = "FI", cex = 5, colour = "#A00E00") +
  annotate("text", x = 2.3, y = -0.7, label = "UA", cex = 5, colour = "#0086A8") +
  annotate("text", x = -1.1, y = -0.9, label = "DK", cex = 5, colour = "#D04E00") +
  annotate("text", x = -0.5, y = 0.95, label = "PT", cex = 5, colour = "#7BA354") +
  annotate("text", x = -1.3, y = -0.1, label = "TR", cex = 5, colour = "#132B69") +
  annotate("text", x = -1.7, y = 0.6, label = "ES", cex = 5, colour = "#B82E00") +
  theme(legend.position = "none")
# labels can be added later

plot_grid(p_F9_v1,p_F9_c1)
plot_grid(p_F9_v2,p_F9_c2)

###########################
##
## Male 8 or 9 traits PCA
##
###########################

# Compare effect of diet control on PCA
pca_M8<-PCA(dataall_M8,scale.unit=TRUE,graph=FALSE,quali.sup=c(1:3),axes=c(1,2))
pcactl_M8<-PCA(datactl_M8,scale.unit=TRUE,graph=FALSE,quali.sup=c(1:3),axes=c(1,2))
colnames(pcactl_M8$var$coord) <- c("ctl_dim1", "ctl_dim2", "ctl_dim3", "ctl_dim4", "ctl_dim5")  

# Tucker's congruence between PCs of the all data vs control data
factcor_M8<-factor.congruence(pca_M8$var$coord, pcactl_M8$var$coord,2)
# Good congruence of PC1 and PC2, moderate congruence of PC3
corrplot(factcor_M8, method="number")

# PCA 1 = M9 (columns 1-3 are qualitative)
pca_M9<-PCA(dataall_M9,scale.unit=TRUE,graph=FALSE,quali.sup=c(1:3),axes=c(1,2))
pca_M9$eig
pca_M9$var$coord
fviz_eig(pca_M9)
fviz_contrib(pca_M9, choice = "var", axes = c(1,3), top = 10)

# Here are the individual coordinates for PC1 output with the qualitative vars
pca_M9_ind<-cbind(dataall_M9[,c(1:3)], pca_M9$ind$coord)

p_M9_v1<-plot(pca_M9,choix="var",axes=c(1,2),cex=1.4) + 
  theme(panel.grid = element_blank())
p_M9_v2<-plot(pca_M9,choix="var",axes=c(3,4),cex=1.4) + 
  theme(panel.grid = element_blank())

#############

M9_c1<-fviz(pca_M9, title = "M9 PC1 vs PC2",
          element = "ind",
          habillage = as.factor(dataall_M9$Country),
          geom = c("point","text"),
          label= "quali",
          pointsize = 2,
          pointshape = 18,
          palette = palette,
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          legend.title = "Treatment",invisible="quali")

# labels for countries have to be added manually. I can do this later for this figure if required
p_M9_c1<-M9_c1+theme(text = element_text(size = 16), 
                axis.title = element_text(size = 14),
                axis.text = element_text(size = 14))+
  ylim(c(-3.3,3.3)) +
  annotate("text", x = 2.2, y = 0.7, label = "AT", cex = 5, colour = "#E38800") + 
  annotate("text", x = -0.8, y = -0.2, label = "DE", cex = 5, colour = "#F6C200") +
  annotate("text", x = 2.15, y = 0.3, label = "RU", cex = 5, colour = "#095888") +
  annotate("text", x = -1.7, y = -1.2, label = "FI", cex = 5, colour = "#A00E00") +
  annotate("text", x = 2.3, y = -0.55, label = "UA", cex = 5, colour = "#0086A8") +
  annotate("text", x = 1.5, y = -1.8, label = "DK", cex = 5, colour = "#D04E00") +
  annotate("text", x = 2.4, y = 2.6, label = "PT", cex = 5, colour = "#7BA354") +
  annotate("text", x = -2.2, y = 1.8, label = "TR", cex = 5, colour = "#132B69") +
  annotate("text", x = -0.5, y = 1.4, label = "ES", cex = 5, colour = "#B82E00") +
  theme(legend.position = "none")

M9_c2<-fviz(pca_M9, title = "M9 PC3 vs PC4",
          element = "ind", axes = c(3, 4),
          habillage = as.factor(dataall_M9$Country),
          geom = c("point","text"),
          label= "quali",
          pointsize = 2,
          pointshape = 18,
          palette = palette,
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          legend.title = "Treatment",invisible="quali")

# labels for countries have to be added manually. I can do this later for this figure if required
p_M9_c2<-M9_c2+theme(text = element_text(size = 16),
                axis.title = element_text(size = 14),
                axis.text = element_text(size = 14))+
  annotate("text", x = -1.85, y = -1.2, label = "AT", cex = 5, colour = "#E38800") +
  annotate("text", x = 1.45, y = 1.25, label = "DE", cex = 5, colour = "#F6C200") +
  annotate("text", x = -1.2, y = 0.85, label = "RU", cex = 5, colour = "#095888") +
  annotate("text", x = -0.4, y = 1.9, label = "FI", cex = 5, colour = "#A00E00") +
  annotate("text", x = -1.5, y = -0.1, label = "UA", cex = 5, colour = "#0086A8") +
  annotate("text", x = 0.6, y = -1.6, label = "DK", cex = 5, colour = "#D04E00") +
  annotate("text", x = 0.5, y = 1.3, label = "PT", cex = 5, colour = "#7BA354") +
  annotate("text", x = -0.5, y = -1.0, label = "TR", cex = 5, colour = "#132B69") +
  annotate("text", x = 2.0, y = -0.3, label = "ES", cex = 5, colour = "#B82E00") +
  theme(legend.position = "none")

plot_grid(p_M9_v1,p_M9_c1)
plot_grid(p_M9_v2,p_M9_c2)


##########################################################
###
### Effect of diet plots
###
##########################################################

# F13 diet plot

F13_dietall<-fviz(pca_F13, title = "F13 all data",
             element = "ind",
             habillage =  as.factor(dataall_F13$Country),
             geom = c("point","text"),
             label= "quali",
             pointsize = 2,
             pointshape = 18,
             palette = palette,
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type="confidence",
             legend.title = "Treatment",invisible="quali")

p_F13_dietall<-F13_dietall+theme(text = element_text(size = 16), 
                       axis.title = element_text(size = 14),
                       axis.text = element_text(size = 14))+
  annotate("text", x = 0.6, y = -0.7, label = "AT", cex = 5, colour = "#E38800") + 
  annotate("text", x = -0.8, y = 0.2, label = "DE", cex = 5, colour = "#F6C200") +
  annotate("text", x = -0.2, y = -0.2, label = "RU", cex = 5, colour = "#095888") +
  annotate("text", x = -2.1, y = 0.8, label = "FI", cex = 5, colour = "#A00E00") +
  annotate("text", x = 1.4, y = 1.4, label = "UA", cex = 5, colour = "#0086A8") +
  annotate("text", x = -0.5, y = 1.9, label = "DK", cex = 5, colour = "#D04E00") +
  annotate("text", x = 2.8, y = -0.5, label = "PT", cex = 5, colour = "#7BA354") +
  annotate("text", x = -1.8, y = -1.8, label = "TR", cex = 5, colour = "#132B69") +
  annotate("text", x = -1.55, y = -0.4, label = "ES", cex = 5, colour = "#B82E00") +
  theme(legend.position = "none")

F13_dietctl<-fviz(pcactl_F13, title = "F13 ctrl diet",
             element = "ind", axes = c(1,2),
             habillage =  as.factor(datactl_F13$Country),
             geom = c("point","text"),
             label= "quali",
             pointsize = 2,
             pointshape = 18,
             palette = palette,
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type="confidence",
             legend.title = "Treatment",invisible="quali")

p_F13_dietctl<-F13_dietctl+theme(text = element_text(size = 16), 
                       axis.title = element_text(size = 14),
                       axis.text = element_text(size = 14))+
  annotate("text", x = 2.0, y = 1.4, label = "AT", cex = 5, colour = "#E38800") + 
  annotate("text", x = -0.5, y = 1.9, label = "DE", cex = 5, colour = "#F6C200") +
  annotate("text", x = 0.3, y = 2.0, label = "RU", cex = 5, colour = "#095888") +
  annotate("text", x = -2.9, y = 0.8, label = "FI", cex = 5, colour = "#A00E00") +
  annotate("text", x = 2.0, y = 0.3, label = "UA", cex = 5, colour = "#0086A8") +
  annotate("text", x = -1.15, y = 2.1, label = "DK", cex = 5, colour = "#D04E00") +
  annotate("text", x = 2.9, y = -0.65, label = "PT", cex = 5, colour = "#7BA354") +
  annotate("text", x = -2.0, y = -1.7, label = "TR", cex = 5, colour = "#132B69") +
  annotate("text", x = -1.1, y = -0.6, label = "ES", cex = 5, colour = "#B82E00") +
  theme(legend.position = "none")

plot_grid(p_F13_dietall,p_F13_dietctl)

#############################################

# F9 diet plot

F9_dietall<-fviz(pca_F9, title = "F9 all data",
            element = "ind",
            habillage =  as.factor(dataall_F9$Country),
            geom = c("point","text"),
            label= "quali",
            pointsize = 2,
            pointshape = 18,
            palette = palette,
            addEllipses = TRUE, # Concentration ellipses
            ellipse.type="confidence",
            legend.title = "Treatment",invisible="quali")

p_F9_dietall<-F9_dietall+theme(text = element_text(size = 16), 
                     axis.title = element_text(size = 14),
                     axis.text = element_text(size = 14))+
  annotate("text", x = 1.6, y = 0.1, label = "AT", cex = 5, colour = "#E38800") + 
  annotate("text", x = 1.4, y = 0.4, label = "DE", cex = 5, colour = "#F6C200") +
  annotate("text", x = 1.9, y = 1.0, label = "RU", cex = 5, colour = "#095888") +
  annotate("text", x = -2.1, y = 1.0, label = "FI", cex = 5, colour = "#A00E00") +
  annotate("text", x = 0.9, y = -1.5, label = "UA", cex = 5, colour = "#0086A8") +
  annotate("text", x = 0.1, y = 1.8, label = "DK", cex = 5, colour = "#D04E00") +
  annotate("text", x = 2.8, y = -1.2, label = "PT", cex = 5, colour = "#7BA354") +
  annotate("text", x = -1.8, y = -1.6, label = "TR", cex = 5, colour = "#132B69") +
  annotate("text", x = -1.9, y = -0.2, label = "ES", cex = 5, colour = "#B82E00") +
  theme(legend.position = "none")
# labels can be added later

F9_dietctl<-fviz(pcactl_F9, title = "F9 ctrl diet",
            element = "ind", axes = c(1, 2),
            habillage =  as.factor(datactl_F9$Country),
            geom = c("point","text"),
            label= "quali",
            pointsize = 2,
            pointshape = 18,
            palette = palette,
            addEllipses = TRUE, # Concentration ellipses
            ellipse.type="confidence",
            legend.title = "Treatment",invisible="quali")

p_F9_dietctl<-F9_dietctl+theme(text = element_text(size = 16), 
                     axis.title = element_text(size = 14),
                     axis.text = element_text(size = 14))+
  annotate("text", x = 1.3, y = -2.0, label = "AT", cex = 5, colour = "#E38800") + 
  annotate("text", x = -0.6, y = 1.1, label = "DE", cex = 5, colour = "#F6C200") +
  annotate("text", x = 1.2, y = 0.75, label = "RU", cex = 5, colour = "#095888") +
  annotate("text", x = -3.0, y = 1.0, label = "FI", cex = 5, colour = "#A00E00") +
  annotate("text", x = 1.1, y = 2.2, label = "UA", cex = 5, colour = "#0086A8") +
  annotate("text", x = -1.4, y = -1.4, label = "DK", cex = 5, colour = "#D04E00") +
  annotate("text", x = 2.8, y = -0.9, label = "PT", cex = 5, colour = "#7BA354") +
  annotate("text", x = -1.2, y = 1.2, label = "TR", cex = 5, colour = "#132B69") +
  annotate("text", x = 1.4, y = 1.0, label = "ES", cex = 5, colour = "#B82E00") +
  theme(legend.position = "none")
# labels can be added later

plot_grid(p_F9_dietall,p_F9_dietctl)

########################################

# M9 diet plot

M8_dietall<-fviz(pca_M8, title = "M8 all data",
            element = "ind",
            habillage = as.factor(dataall_M8$Country),
            geom = c("point","text"),
            label= "quali",
            pointsize = 2,
            pointshape = 18,
            palette = palette,
            addEllipses = TRUE, # Concentration ellipses
            ellipse.type="confidence",
            legend.title = "Treatment",invisible="quali")

# labels for countries have to be added manually. I can do this later for this figure if required
p_M8_dietall<-M8_dietall+theme(text = element_text(size = 16), 
                     axis.title = element_text(size = 14),
                     axis.text = element_text(size = 14))+
  annotate("text", x = 2.2, y = 0.6, label = "AT", cex = 5, colour = "#E38800") + 
  annotate("text", x = -0.8, y = -0.2, label = "DE", cex = 5, colour = "#F6C200") +
  annotate("text", x = 1.9, y = -0.4, label = "RU", cex = 5, colour = "#095888") +
  annotate("text", x = -1.3, y = -1.2, label = "FI", cex = 5, colour = "#A00E00") +
  annotate("text", x = 2.0, y = 1.0, label = "UA", cex = 5, colour = "#0086A8") +
  annotate("text", x = 1.5, y = -1.0, label = "DK", cex = 5, colour = "#D04E00") +
  annotate("text", x = 1.6, y = 2.0, label = "PT", cex = 5, colour = "#7BA354") +
  annotate("text", x = -3.5, y = 0.2, label = "TR", cex = 5, colour = "#132B69") +
  annotate("text", x = -1.0, y = 1.4, label = "ES", cex = 5, colour = "#B82E00") +
  theme(legend.position = "none")

M8_dietctl<-fviz(pcactl_M8, title = "M8 ctrl diet",
            element = "ind", axes = c(1, 2),
            habillage = as.factor(datactl_M8$Country),
            geom = c("point","text"),
            label= "quali",
            pointsize = 2,
            pointshape = 18,
            palette = palette,
            addEllipses = TRUE, # Concentration ellipses
            ellipse.type="confidence",
            legend.title = "Treatment",invisible="quali")

# labels for countries have to be added manually. I can do this later for this figure if required
p_M8_dietctl<-M8_dietctl+theme(text = element_text(size = 16),
                     axis.title = element_text(size = 14),
                     axis.text = element_text(size = 14))+
  annotate("text", x = 1.4, y = 0.7, label = "AT", cex = 5, colour = "#E38800") + 
  annotate("text", x = -0.8, y = -0.4, label = "DE", cex = 5, colour = "#F6C200") +
  annotate("text", x = 1.9, y = -0.4, label = "RU", cex = 5, colour = "#095888") +
  annotate("text", x = -1.3, y = -1.2, label = "FI", cex = 5, colour = "#A00E00") +
  annotate("text", x = 1.7, y = 1.3, label = "UA", cex = 5, colour = "#0086A8") +
  annotate("text", x = 1.5, y = -1.5, label = "DK", cex = 5, colour = "#D04E00") +
  annotate("text", x = 0.5, y = 2.2, label = "PT", cex = 5, colour = "#7BA354") +
  annotate("text", x = -3.5, y = 0.2, label = "TR", cex = 5, colour = "#132B69") +
  annotate("text", x = -0.8, y = 1.6, label = "ES", cex = 5, colour = "#B82E00") +
  theme(legend.position = "none")

plot_grid(p_M8_dietall,p_M8_dietctl)

##########################################################

saveRDS(pca_F13, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/pca_F13.rds")
saveRDS(pca_F9, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/pca_F9.rds")
saveRDS(pca_M9, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/pca_M9.rds")
saveRDS(pca_M8, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/pca_M8.rds")

saveRDS(pcactl_F13, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/pcactl_F13.rds")
saveRDS(pcactl_F9, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/pcactl_F9.rds")
saveRDS(pcactl_M8, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/pcactl_M8.rds")

##########################################################

write.csv(pca_F13_ind, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/pca_F13_PCcoords.csv", row.names = F)
write.csv(pca_F9_ind, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/pca_F9_PCcoords.csv", row.names = F)
write.csv(pca_M9_ind, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/pca_M9_PCcoords.csv", row.names = F)

save.image(file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/PCA_drosEU_240603.RData")
