rm(list=ls(all=TRUE))
library(dplyr)
library(FactoMineR)
library(factoextra)
library(ggforce)

master <- read.csv("/Users/ewanharney/Dropbox/Barcelona_IBE/DrosEU/all_models_line_meta_compound_random_coefs_wide_220916.csv")

# There are several different ways to do the PCA, depending on how we consider sex effects:

# 1). use all the traits (M and F as separate variables - will be skewed towards variables with M + F data, as these are represented more times), 
# 2). average across sexes (if a trait has M + F values, take the mean), 
# 3). combine sexes into a single variable (separate value for each sex), 
# 4). analyse the sexes separately: A) Females variables, B) Male variables C) Female variables plus DevTime_pupae and viability (which have no sex) as these are interesting and it seems a shame to lose them!

# Make variables for Country and the sexes
Country <- c(rep("Turkey", 20),rep("Portugal", 17),rep("Spain", 15),rep("Germany", 20),rep("Austria", 20),rep("Ukraine", 19),rep("Denmark", 20),rep("Russia", 20),rep("Finland", 22))
SexM <-c(rep("M", 173))
SexF <-c(rep("F", 173))

####
#### Formatting the data for the different versions of the PCA
####
#################################
# 1). All the traits = 23 quantitative variables.
# Note we drop locomotor activity, as there are many missing values
# Traits with mane and female values are over-represented - I don't think this is the best option!
# I have also made a decision to only include one pigmentation trait (Pgm_Total_F) and one set of wing areas (WA_L_F and WA_L_M).
master1 <- master[,c("Population","Line","CCRT_F","CCRT_M","CSM_F","CSM_M","DT_A_F","DT_A_M","DT_P_NA","Dia_F","DW_F","DW_M",
                      "Fec_F","HSM_F","HSM_M","LS_F","LS_M","Pgm_Total_F","SR_F","SR_M","TL_F","TL_M","Via_NA","WA_L_F","WA_L_M")]
data1 <- na.omit(cbind(Country, master1))

#################################
# 2). average across sexes (if there are data for both sexes) = 14 quant variables
# Again, drop locomotor acrivity, just focus on WA_L
CCRT_FM<-(master$CCRT_F+master$CCRT_M)/2
CSM_FM<-(master$CSM_F+master$CSM_M)/2
DT_A_FM<-(master$DT_A_F+master$DT_A_M)/2
DW_FM<-(master$DW_F+master$DW_M)/2
HSM_FM<-(master$HSM_F+master$HSM_M)/2
LS_FM<-(master$LS_F+master$LS_M)/2
SR_FM<-(master$SR_F+master$SR_M)/2
TL_FM<-(master$TL_F+master$TL_M)/2
WA_L_FM<-(master$WA_L_F+master$WA_L_M)/2
master2 <- cbind(Country,master[,c("Population","Line")],CCRT_FM,CSM_FM,DT_A_FM,DW_FM,HSM_FM,LS_FM,SR_FM,TL_FM,WA_L_FM)
# we also include traits that are female only, or without sex
data2 <- na.omit(cbind(master2,master[,c("DT_P_NA","Dia_F","Fec_F","Pgm_Total_F","Via_NA")]))

#################################
# 3). Combine sexes into a single variable = 9 quant variables
# Again, drop locomotor acrivity, just focus on WA_L
mastersubF <- master[,c("Population","Line","CCRT_F","CSM_F","DT_A_F","DW_F","HSM_F","LS_F","SR_F","TL_F","WA_L_F")]
mastersubM <- master[,c("Population","Line","CCRT_M","CSM_M","DT_A_M","DW_M","HSM_M","LS_M","SR_M","TL_M","WA_L_M")]
datasubF <- na.omit(cbind(SexF,Country,mastersubF))
datasubM <- na.omit(cbind(SexM,Country,mastersubM))
# the column names need to be renamed and match so that we can do the row bind
colnames(datasubF) <- c("Sex", "Country", "Population", "Line", "CCRT", "CSM", "DT_A", "DW", "HSM", "LS", "SR", "TL", "WA_L")
colnames(datasubM) <- c("Sex", "Country", "Population", "Line", "CCRT", "CSM", "DT_A", "DW", "HSM", "LS", "SR", "TL", "WA_L" )
data3 <- rbind(datasubF, datasubM)
# make another label that combines sex and country
data3$Country_Sex <-paste(data3$Country,data3$Sex,sep = "_")

#################################
# 4). analyse the sexes separately: 
# A) Females variables - 12 quant traits
masterallF <- master[,c("Population","Line","CCRT_F","CSM_F","DT_A_F","Dia_F","DW_F","Fec_F","HSM_F","LS_F","Pgm_Total_F","SR_F","TL_F","WA_L_F")]
data4A <- na.omit(cbind(Country,masterallF))
# B) Male variables - 9 quant traits
data4B <- na.omit(cbind(Country,mastersubM))
# C) Female variables plus DevTime_pupae and viability - 14 quant traits
masterallFplus <- master[,c("Population","Line","CCRT_F","CSM_F","DT_A_F","Dia_F","DW_F","Fec_F","HSM_F","LS_F","Pgm_Total_F","SR_F","TL_F","WA_L_F", "DT_P_NA", "Via_NA")]
data4C <- na.omit(cbind(Country,masterallFplus))

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
# PCA 1 (columns 1-3 are qualitative)
pca1<-PCA(data1,scale.unit=TRUE,graph=FALSE,quali.sup=c(1:3),axes=c(1,2))
pca1$eig
pca1$var$coord
# Here are the individual coordinates for PC1 output with the qualitative vars
PCA1_ind<-cbind(data1[,c(1:3)], pca1$ind$coord)

plot(pca1,choix="var",axes=c(1,2),cex=1.4)
plot(pca1,choix="var",axes=c(2,3),cex=1.4)
# simple ellipses with chosen colours
plotellipses(pca1, keepvar = c(1), axes = c(1, 2),label = "quali", level = 0.95, palette=palette)
p1a<-fviz(pca1, title = "PCA - all 23 traits",
            element = "ind",
                     habillage =  data1$Country,
                     geom = c("point","text"),
                     label= "quali",
                     pointsize = 2,
                     pointshape = 18,
                     palette = palette,
                     addEllipses = TRUE, # Concentration ellipses
                     ellipse.type="confidence",
                     legend.title = "Treatment",invisible="quali"
)
# labels for countries have to be added manually. I can do this later for this figure if required
p1a+theme(text = element_text(size = 16), 
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14))



#####################
# PCA 2 (columns 1-3 are qualitative)
pca2<-PCA(data2,scale.unit=TRUE,graph=FALSE,quali.sup=c(1:3),axes=c(1,2))
pca2$eig
pca2$var$coord
# Here are the individual coordinates for PC2 output with the qualitative vars
PCA2_ind<-cbind(data2[,c(1:3)], pca2$ind$coord)

plot(pca2,choix="var",axes=c(1,2),cex=1.4)
plot(pca2,choix="var",axes=c(2,3),cex=1.4)
# simple ellipses with chosen colours
plotellipses(pca2, keepvar = c(1), axes = c(1, 2),label = "quali", level = 0.95, palette=palette)
p2a<-fviz(pca2, title = "PCA - 14 traits (ignoring Sex)",
          element = "ind",
          habillage =  data2$Country,
          geom = c("point","text"),
          label= "quali",
          pointsize = 2,
          pointshape = 18,
          palette = palette,
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          legend.title = "Treatment",invisible="quali"
)
p2a+theme(text = element_text(size = 16), 
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14))+
  annotate("text", x = 1.8, y = 1.6, label = "Austria", cex = 5, colour = "#E38800") + 
  annotate("text", x = -1, y = -0.5, label = "Germany", cex = 5, colour = "#F6C200") +
  annotate("text", x = 0.5, y = -0.4, label = "Russia", cex = 5, colour = "#095888") +
  annotate("text", x = -2.7, y = 0.3, label = "Finland", cex = 5, colour = "#A00E00") +
  annotate("text", x = 0, y = 2.2, label = "Ukraine", cex = 5, colour = "#0086A8") +
  annotate("text", x = -2.7, y = 1.2, label = "Denmark", cex = 5, colour = "#D04E00") +
  annotate("text", x = 2.5, y = -1.2, label = "Portugal", cex = 5, colour = "#7BA354") +
  annotate("text", x = -0.6, y = -3.6, label = "Turkey", cex = 5, colour = "#132B69") +
  annotate("text", x = 1.4, y = -2.1, label = "Spain", cex = 5, colour = "#B82E00") 
# labels have been added

#####################
# PCA 3 (columns 1-4 and 14 are qualitative)
pca3<-PCA(data3,scale.unit=TRUE,graph=FALSE,quali.sup=c(1:4,14),axes=c(1,2))
pca3$eig
pca3$var$coord
# Here are the individual coordinates for PC3 output with the qualitative vars
PCA3_ind<-cbind(data3[,c(1:4,14)], pca3$ind$coord)
plot(pca3,choix="var",axes=c(1,2),cex=1.4)
plot(pca3,choix="var",axes=c(2,3),cex=1.4)
# simple ellipses with chosen colours
plotellipses(pca3, keepvar = c(2), axes = c(1, 2),label = "quali", level = 0.95, palette=palette)
plotellipses(pca3, keepvar = c(2), axes = c(3, 2),label = "quali", level = 0.95, palette=palette)
# note that axes 2 and 3 are not in the usual order (it makes comparison with other output easier)
p3a<-fviz(pca3, title = "PCA - 9 Female and Male traits",
          element = "ind", axes = c(3, 2),
          habillage =  data3$Country,
          geom = c("point","text"),
          label= "quali",
          pointsize = 2,
          pointshape = 18,
          alpha =0.5,
          palette = palette,
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          legend.title = "Treatment",invisible="quali"
)
p3a+theme(text = element_text(size = 16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14))+
  annotate("text", x = -1.0, y = -0.2, label = "Austria", cex = 5, colour = "#E38800") + 
  annotate("text", x = 1.1, y = 0.9, label = "Germany", cex = 5, colour = "#F6C200") +
  annotate("text", x = 1.1, y = -0.7, label = "Russia", cex = 5, colour = "#095888") +
  annotate("text", x = 1.7, y = -0.3, label = "Finland", cex = 5, colour = "#A00E00") +
  annotate("text", x = 1.6, y = 0.6, label = "Ukraine", cex = 5, colour = "#0086A8") +
  annotate("text", x = -1.5, y = -1, label = "Denmark", cex = 5, colour = "#D04E00") +
  annotate("text", x = 0.4, y = 1.2, label = "Portugal", cex = 5, colour = "#7BA354") +
  annotate("text", x = -1.7, y = 1.4, label = "Turkey", cex = 5, colour = "#132B69") +
  annotate("text", x = -1.2, y = 0.4, label = "Spain", cex = 5, colour = "#B82E00") 


p3b<-fviz(pca3, title = "PCA - 9 Female and Male traits",
          element = "ind", axes = c(1, 2),
          habillage = as.factor(data3$Country_Sex),
          geom = c("point","text"),
          label= "quali",
          pointsize = 2,
          pointshape = 18,
          alpha =0.5,
          palette = c("#E38800","#FFBC58","#D04E00","#FF904F","#A00E00","#E04737","#F6C200","#FFDD60","#7BA354",
                      "#CFE7B8","#095888","#3E81A8","#B82E00","#F26E41","#132B69","#3C5491","#0086A8","#06CCFB"),
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          legend.title = "Treatment",invisible="quali"
)
p3b+theme(text = element_text(size = 16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14))

p3c<-fviz(pca3, title = "PCA - 9 Female and Male traits",
          element = "ind", axes = c(3, 2),
          habillage = as.factor(data3$Country_Sex),
          geom = c("point","text"),
          label= "quali",
          pointsize = 2,
          pointshape = 18,
          alpha =0.5,
          palette = c("#E38800","#FFBC58","#D04E00","#FF904F","#A00E00","#E04737","#F6C200","#FFDD60","#7BA354",
                      "#CFE7B8","#095888","#3E81A8","#B82E00","#F26E41","#132B69","#3C5491","#0086A8","#06CCFB"),
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          legend.title = "Treatment",invisible="quali"
)
p3c+theme(text = element_text(size = 16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14))

p3d<-fviz(pca3, title = "PCA - 9 Female and Male traits",
          element = "ind", axes = c(3, 2),
          habillage = as.factor(data3$Country_Sex),
          geom = c("point","text"),
          label= "quali",
          pointsize = 3,
          pointshape = 18,
          alpha =0.2,
          palette = c("#E38800","#FFBC58","#D04E00","#FF904F",
                      "#A00E00","#E04737","#F6C200","#FFDD60","#7BA354","#CFE7B8",
                      "#095888","#3E81A8","#B82E00","#F26E41","#132B69","#3C5491",
                      "#0086A8","#06CCFB"),
          #addEllipses = TRUE, # Concentration ellipses
          #ellipse.type="confidence",
          legend.title = "Treatment"#,invisible="quali"
)

# Looking at how sexes differ in a bit more detail - have to redo pca with just 1 qualitative variable
pca3b<-PCA(data3[,c(5:14)],scale.unit=TRUE,graph=FALSE,quali.sup=c(10),axes=c(1,2))
QS=cbind.data.frame(rownames(pca3b$quali.sup$coord),pca3b$quali.sup$coord)

p3d+theme(text = element_text(size = 16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14))+
  geom_segment(aes(x = QS[2,4], y = QS[2,3], xend = QS[1,4], yend = QS[1,3]),arrow=arrow(angle = 30, length = unit(0.3, "cm"),ends = "last", type = "open"))+
  geom_segment(aes(x = QS[4,4], y = QS[4,3], xend = QS[3,4], yend = QS[3,3]),arrow=arrow(angle = 30, length = unit(0.3, "cm"),ends = "last", type = "open"))+
  geom_segment(aes(x = QS[6,4], y = QS[6,3], xend = QS[5,4], yend = QS[5,3]),arrow=arrow(angle = 30, length = unit(0.3, "cm"),ends = "last", type = "open"))+
  geom_segment(aes(x = QS[8,4], y = QS[8,3], xend = QS[7,4], yend = QS[7,3]),arrow=arrow(angle = 30, length = unit(0.3, "cm"),ends = "last", type = "open"))+
  geom_segment(aes(x = QS[10,4], y = QS[10,3], xend = QS[9,4], yend = QS[9,3]),arrow=arrow(angle = 30, length = unit(0.3, "cm"),ends = "last", type = "open"))+
  geom_segment(aes(x = QS[12,4], y = QS[12,3], xend = QS[11,4], yend = QS[11,3]),arrow=arrow(angle = 30, length = unit(0.3, "cm"),ends = "last", type = "open"))+
  geom_segment(aes(x = QS[14,4], y = QS[14,3], xend = QS[13,4], yend = QS[13,3]),arrow=arrow(angle = 30, length = unit(0.3, "cm"),ends = "last", type = "open"))+
  geom_segment(aes(x = QS[16,4], y = QS[16,3], xend = QS[15,4], yend = QS[15,3]),arrow=arrow(angle = 30, length = unit(0.3, "cm"),ends = "last", type = "open"))+
  geom_segment(aes(x = QS[18,4], y = QS[18,3], xend = QS[17,4], yend = QS[17,3]),arrow=arrow(angle = 30, length = unit(0.3, "cm"),ends = "last", type = "open"))+
  annotate("text", x = 1.7, y = -0.3, label = "Finland", cex = 5, colour = "#A00E00") +
  annotate("text", x = 0.9, y = -0.8, label = "Ukraine", cex = 5, colour = "#0086A8") +
  annotate("text", x = -1.5, y = 1.8, label = "Turkey", cex = 5, colour = "#132B69") +
  annotate("text", x = -0.1, y = 2.0, label = "Spain", cex = 5, colour = "#B82E00") 


#####################
# PCA 4A (columns 1-3 qualitative)
pca4A<-PCA(data4A,scale.unit=TRUE,graph=FALSE,quali.sup=c(1:3),axes=c(1,2))
pca4A$eig
pca4A$var$coord
# Here are the individual coordinates for PC4A output with the qualitative vars
PCA4A_ind<-cbind(data4A[,c(1:3)], pca4A$ind$coord)

plot(pca4A,choix="var",axes=c(1,2),cex=1.4)
# simple ellipses with chosen colours
plotellipses(pca4A, keepvar = c(1), axes = c(1, 2),label = "quali", level = 0.95, palette=palette)
p4Aa<-fviz(pca4A, title = "PCA - 12 Female traits",
          element = "ind", axes = c(1, 2),
          habillage =  data4A$Country,
          geom = c("point","text"),
          label= "quali",
          pointsize = 2,
          pointshape = 18,
          alpha =0.5,
          palette = palette,
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          legend.title = "Treatment",invisible="quali"
)
p4Aa+theme(text = element_text(size = 16),
           axis.title = element_text(size = 14),
           axis.text = element_text(size = 14))
# labels for countries have to be added manually. I can do this later for this figure if required

#####################
# PCA 4B (columns 1-3 are qualitative)
pca4B<-PCA(data4B,scale.unit=TRUE,graph=FALSE,quali.sup=c(1:3),axes=c(1,2))
pca4B$eig
pca4B$var$coord
# Here are the individual coordinates for PC4B output with the qualitative vars
PCA4B_ind<-cbind(data4B[,c(1:3)], pca4B$ind$coord)

plot(pca4B,choix="var",axes=c(1,2),cex=1.4)
# simple ellipses with chosen colours
plotellipses(pca4B, keepvar = c(1), axes = c(1, 2),label = "quali", level = 0.95, palette=palette)
p4Ba<-fviz(pca4B, title = "PCA - 9 Male traits",
          element = "ind", axes = c(3, 2),
          habillage =  data4B$Country,
          geom = c("point","text"),
          label= "quali",
          pointsize = 2,
          pointshape = 18,
          alpha =0.5,
          palette = palette,
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          legend.title = "Treatment",invisible="quali"
)
p4Ba+theme(text = element_text(size = 16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14))
# labels for countries have to be added manually. I can do this later for this figure if required

#####################
# PCA 4C (columns 1-3 are qualitative)
pca4C<-PCA(data4C,scale.unit=TRUE,graph=FALSE,quali.sup=c(1:3),axes=c(1,2))
pca4C$eig
pca4C$var$coord
# Here are the individual coordinates for PC4C output with the qualitative vars
PCA4C_ind<-cbind(data4C[,c(1:3)], pca4C$ind$coord)

plot(pca4C,choix="var",axes=c(1,2),cex=1.4)
# simple ellipses with chosen colours
plotellipses(pca4C, keepvar = c(1), axes = c(1, 2),label = "quali", level = 0.95, palette=palette)
p4Ca<-fviz(pca4C, title = "PCA - 12 Female traits + 2 extra non-sexed traits",
          element = "ind", axes = c(3, 2),
          habillage =  data4C$Country,
          geom = c("point","text"),
          label= "quali",
          pointsize = 2,
          pointshape = 18,
          alpha =0.5,
          palette = palette,
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          legend.title = "Treatment",invisible="quali"
)
p4Ca+theme(text = element_text(size = 16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14))

#######
save(pca1, file = "/Users/ewanharney/Dropbox/Barcelona_IBE/DrosEU/PCA1_drosEU.RData")
save(pca2, file = "/Users/ewanharney/Dropbox/Barcelona_IBE/DrosEU/PCA2_drosEU.RData")
save(pca3, file = "/Users/ewanharney/Dropbox/Barcelona_IBE/DrosEU/PCA3_drosEU.RData")
save(pca4A, file = "/Users/ewanharney/Dropbox/Barcelona_IBE/DrosEU/PCA4A_drosEU.RData")
save(pca4B, file = "/Users/ewanharney/Dropbox/Barcelona_IBE/DrosEU/PCA4B_drosEU.RData")
save(pca4C, file = "/Users/ewanharney/Dropbox/Barcelona_IBE/DrosEU/PCA4C_drosEU.RData")

