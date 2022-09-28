setwd("C:/Users/Венера/Dropbox/UoL/other/ClimateComponents/NASA/results/traits_combined/")

library(nasapower)
#Recarei - 05/10/2018
#Gimenells (Lleida) - 25/08/2018
#Karensminde - 05/10/2018
#Munich - 21-30/06/2018
#Mauternbach - 08/09/2018
#Akaa - 20/07/2018
#Uman - 18/08/2018
#Yesiloz - 27/09/2018
#Valday - 20-30/08/2018

#preparing a dataframe
all_climatology <- data.frame(matrix(NA, nrow = 9, ncol = 19))
col_names <- c("Longitude", "Latitude", "Altitude", "Population", "Country", 
               "TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
               "TS_MIN", "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS",
               "PRECTOTCORR", "ALLSKY_SFC_LW_DWN")

colnames(all_climatology) <-c(col_names)

all_climatology$Longitude <- c(-8.41, 0.62, 10.213, 11.61, 15.56, 23.52, 30.206, 32.26, 33.244)
all_climatology$Latitude <- c(41.15, 41.618, 55.945, 48.18, 48.375, 61.1, 48.753, 40.231, 57.979)
all_climatology$Altitude <- c(175, 173, 15, 520, 572, 88, 214, 680, 217)
all_climatology$Population <- c("Recarei", "Gimenells", "Karensminde", "Munich", "Mauternbach", "Akaa", "Uman", "Yesiloz", "Valday")
all_climatology$Country <- c("Portugal", "Spain", "Denmark", "Germany", "Austria", "Finland", "Ukaraine", "Turkey", "Russia")

#Recarei

Recarei_ag <- get_power(
  community = "ag",
  lonlat = c(-8.410, 41.150),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "ALLSKY_SFC_LW_DWN"),
  temporal_api = "climatology"
)
all_climatology[1,6:19] <- Recarei_ag$ANN 

#Gimenells

Gimenells_ag <- get_power(
  community = "ag",
  lonlat = c(0.62, 41.618),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "ALLSKY_SFC_LW_DWN"),
  temporal_api = "climatology"
)
all_climatology[2,6:19] <- Gimenells_ag$ANN 

#Karensminde

Karensminde_ag <- get_power(
  community = "ag",
  lonlat = c(10.213, 55.945),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "ALLSKY_SFC_LW_DWN"),
  temporal_api = "climatology"
)
all_climatology[3,6:19] <- Karensminde_ag$ANN 

#Munich

Munich_ag <- get_power(
  community = "ag",
  lonlat = c(11.61, 48.18),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "ALLSKY_SFC_LW_DWN"),
  temporal_api = "climatology"
)
all_climatology[4,6:19] <- Munich_ag$ANN 

#Mauternbach

Mauternbach_ag <- get_power(
  community = "ag",
  lonlat = c(15.560, 48.375),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "ALLSKY_SFC_LW_DWN"),
  temporal_api = "climatology"
)
all_climatology[5,6:19] <- Mauternbach_ag$ANN 

#Akaa

Akaa_ag <- get_power(
  community = "ag",
  lonlat = c(23.520, 61.100),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "ALLSKY_SFC_LW_DWN"),
  temporal_api = "climatology"
)
all_climatology[6,6:19] <- Akaa_ag$ANN 

#Uman

Uman_ag <- get_power(
  community = "ag",
  lonlat = c(30.206, 48.753),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "ALLSKY_SFC_LW_DWN"),
  temporal_api = "climatology"
)
all_climatology[7,6:19] <- Uman_ag$ANN 

#Yesiloz
Yesiloz_ag <- get_power(
  community = "ag",
  lonlat = c(32.26, 40.231),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "ALLSKY_SFC_LW_DWN"),
  temporal_api = "climatology"
)
all_climatology[8,6:19] <- Yesiloz_ag$ANN 

#Valday
Valday_ag <- get_power(
  community = "ag",
  lonlat = c(33.244, 57.979),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR","ALLSKY_SFC_LW_DWN"),
  temporal_api = "climatology"
)
all_climatology[9,6:19] <- Valday_ag$ANN 

write.csv(all_climatology,"all_param_30y.csv", row.names = F)

#PCA
setwd("C:/Users/Венера/Dropbox/UoL/other/ClimateComponents/NASA/results/traits_combined/")
all_climatology <- read.csv("all_param_30y.csv")
all_climatology <- all_climatology[,6:19]
rownames(all_climatology) <- c("Recarei", "Gimenells", "Karensminde", "Munich", "Mauternbach", "Akaa", "Uman", "Yesiloz", "Valday")

library("FactoMineR")
library("factoextra")
options(ggrepel.max.overlaps = Inf)
all_ann_clim_pca <- PCA(all_climatology, scale.unit = TRUE, graph = TRUE)
print(all_ann_clim_pca)

# matrix with eigenvalues
all_ann_clim_pca$eig
eig.val <- get_eigenvalue(all_ann_clim_pca)
eig.val #first 2 PCs >1

#Scree plot
fviz_eig(all_ann_clim_pca, addlabels = TRUE, xlab = "principal components")

fviz_pca_ind(all_ann_clim_pca, axes = c(1, 2), col.ind = "coord", 
             title = "Principal Component Analysis",
             subtitle = "Bioclimatic variables",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)

#PCs (scores)
all_ann_clim_pca$ind$coord

#Contribution to dimentions
# Change the gradient color
fviz_pca_var(all_ann_clim_pca, col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# Contributions of variables to PC1
fviz_contrib(all_ann_clim_pca, choice = "var", axes = 1, top = 19, title="Contribution of variables to PC1")

# Contributions of variables to PC2
fviz_contrib(all_ann_clim_pca, choice = "var", axes = 2, top = 19, title="Contribution of variables to PC2")

#pca result dataframe
all_climatology <- read.csv("all_param_30y.csv")
pops_data <- all_climatology[,1:5]
PC1_2 <- all_ann_clim_pca$ind$coord[1:9, 1:2] #there are 2 PCs with eigvalues >1
all_ann_clim_results <- cbind(pops_data, PC1_2)
write.csv(all_ann_clim_results,"all_30y_PCA.csv", row.names =F)

library(corrplot)
cor_df <- all_ann_clim_results[ -c(4:5) ]
corrplot(cor(cor_df), method = 'number') 
corrplot(cor(cor_df), type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)


#phenotype analyses
#response as output from the Ewan's PCA for:
#females 9 traits
#males 9 traits
#females max
#females max+

#
#females 9 traits
#
setwd("/Users/Венера/Dropbox/UoL/other/ClimateComponents/NASA/results/traits_combined/")
F9_all_PC <- read.csv("/Users/Венера/Dropbox/UoL/other/ClimateComponents/NASA/results/Ewan/F9_drosEU_PCcoords.csv")
F9_data <- F9_all_PC[,1:5]
colnames(F9_data) <- c("Country", "Population", "Line", "PC1_F9", "PC2_F9")
climate_data <- read.csv("all_30y_PCA.csv")
clim <- climate_data[ -c(5) ]

new_cols <- c("Longitude", "Latitude", "Altitude","Location", "PC1_clim", "PC2_clim")
F9_data[ , new_cols] <- NA
for (i in 1:nrow(F9_data)){
  if (F9_data[i,"Population"] == "RE"){ #Portugal
    F9_data[i, 6:11] <- clim[1,]
  }
  if (F9_data[i,"Population"] == "GI"){ #Spain
    F9_data[i, 6:11] <- clim[2,]
  }
  if (F9_data[i,"Population"] == "KA"){ #Denmark
    F9_data[i, 6:11] <- clim[3,]
  }
  if (F9_data[i,"Population"] == "MU"){ #Germany
    F9_data[i, 6:11] <- clim[4,]
  }
  if (F9_data[i,"Population"] == "MA"){ #Austria
    F9_data[i, 6:11] <- clim[5,]
  }
  if (F9_data[i,"Population"] == "AK"){ #Finland
    F9_data[i, 6:11] <- clim[6,]
  }
  if (F9_data[i,"Population"] == "UM"){ #Ukraine
    F9_data[i, 6:11] <- clim[7,]
  }
  if (F9_data[i,"Population"] == "YE"){ #Turkey
    F9_data[i, 6:11] <- clim[8,]
  }
  if (F9_data[i,"Population"] == "VA"){ #Russia
    F9_data[i, 6:11] <- clim[9,]
  }
} 
write.csv(F9_data,"F9_30y_data.csv", row.names =F)

#phenotypes
#PC1
F9_PC1_1 <- lmer(F9_data$PC1_F9 ~ F9_data$PC1_clim + F9_data$PC2_clim +(1|(as.factor(F9_data$Population))),data = F9_data)
summary(F9_PC1_1)
F9_PC1_2 <- lmer(F9_data$PC1_F9 ~ F9_data$PC2_clim +(1|(as.factor(F9_data$Population))),data = F9_data)
summary(F9_PC1_2)

F9_PC1_3 <- glm(F9_data$PC1_F9 ~ F9_data$PC1_clim + F9_data$PC2_clim,data = F9_data)
summary(F9_PC1_3)
F9_PC1_4 <- glm(F9_data$PC1_F9 ~ F9_data$PC1_clim + F9_data$PC2_clim,data = F9_data)
summary(F9_PC1_4)
anova(F9_PC1_1,F9_PC1_2,F9_PC1_3,F9_PC1_4, test="Chisq")


#PC2
F9_PC2_1 <- lmer(F9_data$PC2_F9 ~ F9_data$PC1_clim + F9_data$PC2_clim +(1|(as.factor(F9_data$Population))),data = F9_data)
summary(F9_PC2_1)
F9_PC2_2 <- lmer(F9_data$PC2_F9 ~ F9_data$PC2_clim +(1|(as.factor(F9_data$Population))),data = F9_data)
summary(F9_PC2_2)

F9_PC2_3 <- glm(F9_data$PC2_F9 ~ F9_data$PC1_clim + F9_data$PC2_clim,data = F9_data)
summary(F9_PC2_3)
anova(F9_PC2_1,F9_PC2_2,F9_PC2_3, test="Chisq")


#
#males 9 traits
#
setwd("/Users/Венера/Dropbox/UoL/other/ClimateComponents/NASA/results/traits_combined/")
M9_all_PC <- read.csv("/Users/Венера/Dropbox/UoL/other/ClimateComponents/NASA/results/Ewan/M9_drosEU_PCcoords.csv")
M9_data <- M9_all_PC[,1:5]
colnames(M9_data) <- c("Country", "Population", "Line", "PC1_M9", "PC2_M9")
climate_data <- read.csv("all_30y_PCA.csv")
clim <- climate_data[ -c(5) ]

new_cols <- c("Longitude", "Latitude", "Altitude","Location", "PC1_clim", "PC2_clim")
M9_data[ , new_cols] <- NA
for (i in 1:nrow(M9_data)){
  if (M9_data[i,"Population"] == "RE"){ #Portugal
    M9_data[i, 6:11] <- clim[1,]
  }
  if (M9_data[i,"Population"] == "GI"){ #Spain
    M9_data[i, 6:11] <- clim[2,]
  }
  if (M9_data[i,"Population"] == "KA"){ #Denmark
    M9_data[i, 6:11] <- clim[3,]
  }
  if (M9_data[i,"Population"] == "MU"){ #Germany
    M9_data[i, 6:11] <- clim[4,]
  }
  if (M9_data[i,"Population"] == "MA"){ #Austria
    M9_data[i, 6:11] <- clim[5,]
  }
  if (M9_data[i,"Population"] == "AK"){ #Finland
    M9_data[i, 6:11] <- clim[6,]
  }
  if (M9_data[i,"Population"] == "UM"){ #Ukraine
    M9_data[i, 6:11] <- clim[7,]
  }
  if (M9_data[i,"Population"] == "YE"){ #Turkey
    M9_data[i, 6:11] <- clim[8,]
  }
  if (M9_data[i,"Population"] == "VA"){ #Russia
    M9_data[i, 6:11] <- clim[9,]
  }
} 

write.csv(M9_data,"M9_30y_data.csv", row.names =F)

#phenotypes
#PC1
M9_PC1_1 <- lmer(M9_data$PC1_M9 ~ M9_data$PC1_clim + M9_data$PC2_clim +(1|(as.factor(M9_data$Population))),data = M9_data)
summary(M9_PC1_1)
M9_PC1_2 <- lmer(M9_data$PC1_M9 ~ M9_data$PC2_clim +(1|(as.factor(M9_data$Population))),data = M9_data)
summary(M9_PC1_2)

M9_PC1_3 <- glm(M9_data$PC1_M9 ~ M9_data$PC1_clim + M9_data$PC2_clim,data = M9_data)
summary(M9_PC1_3)
M9_PC1_4 <- glm(M9_data$PC1_M9 ~ M9_data$PC2_clim,data = M9_data)
summary(M9_PC1_4)

anova(M9_PC1_1,M9_PC1_2,M9_PC1_3,M9_PC1_4, test="Chisq")

#PC2
M9_PC2_1 <- lmer(M9_data$PC2_M9 ~ M9_data$PC1_clim + M9_data$PC2_clim +(1|(as.factor(M9_data$Population))),data = M9_data)
summary(M9_PC2_1)
M9_PC2_2 <- lmer(M9_data$PC2_M9 ~ M9_data$PC1_clim +(1|(as.factor(M9_data$Population))),data = M9_data)
summary(M9_PC2_2)

M9_PC2_3 <- glm(M9_data$PC2_M9 ~ M9_data$PC1_clim + M9_data$PC2_clim,data = M9_data)
summary(M9_PC2_3)

anova(M9_PC2_1,M9_PC2_2,M9_PC2_3, test="Chisq")

#females max
setwd("/Users/Венера/Dropbox/UoL/other/ClimateComponents/NASA/results/traits_combined/")
F9max_all_PC <- read.csv("/Users/Венера/Dropbox/UoL/other/ClimateComponents/NASA/results/Ewan/Fmax_drosEU_PCcoords.csv")
F9max_data <- F9max_all_PC[,1:5]
colnames(F9max_data) <- c("Country", "Population", "Line", "PC1_F9max", "PC2_F9max")
climate_data <- read.csv("all_30y_PCA.csv")
clim <- climate_data[ -c(5) ]

new_cols <- c("Longitude", "Latitude", "Altitude","Location", "PC1_clim", "PC2_clim")
F9max_data[ , new_cols] <- NA
for (i in 1:nrow(F9max_data)){
  if (F9max_data[i,"Population"] == "RE"){ #Portugal
    F9max_data[i, 6:11] <- clim[1,]
  }
  if (F9max_data[i,"Population"] == "GI"){ #Spain
    F9max_data[i, 6:11] <- clim[2,]
  }
  if (F9max_data[i,"Population"] == "KA"){ #Denmark
    F9max_data[i, 6:11] <- clim[3,]
  }
  if (F9max_data[i,"Population"] == "MU"){ #Germany
    F9max_data[i, 6:11] <- clim[4,]
  }
  if (F9max_data[i,"Population"] == "MA"){ #Austria
    F9max_data[i, 6:11] <- clim[5,]
  }
  if (F9max_data[i,"Population"] == "AK"){ #Finland
    F9max_data[i, 6:11] <- clim[6,]
  }
  if (F9max_data[i,"Population"] == "UM"){ #Ukraine
    F9max_data[i, 6:11] <- clim[7,]
  }
  if (F9max_data[i,"Population"] == "YE"){ #Turkey
    F9max_data[i, 6:11] <- clim[8,]
  }
  if (F9max_data[i,"Population"] == "VA"){ #Russia
    F9max_data[i, 6:11] <- clim[9,]
  }
} 
write.csv(F9max_data,"F9max_30y_data.csv", row.names =F)


#phenotypes
library(lme4)
library(lmerTest)
#PC1
F9max_PC1_1 <- lmer(F9max_data$PC1_F9max ~ F9max_data$PC1_clim + F9max_data$PC2_clim +(1|(as.factor(F9max_data$Population))),data = F9max_data)
summary(F9max_PC1_1)
F9max_PC1_2 <- lmer(F9max_data$PC1_F9max ~ F9max_data$PC1_clim +(1|(as.factor(F9max_data$Population))),data = F9max_data)
summary(F9max_PC1_2)

F9max_PC1_3 <- glm(F9max_data$PC1_F9max ~ F9max_data$PC1_clim + F9max_data$PC2_clim,data = F9max_data)
summary(F9max_PC1_3)
F9max_PC1_4 <- glm(F9max_data$PC1_F9max ~ F9max_data$PC1_clim,data = F9max_data)
summary(F9max_PC1_4)

anova(F9max_PC1_1,F9max_PC1_2,F9max_PC1_3, F9max_PC1_4, test="Chisq")


#PC2
F9max_PC2_1 <- lmer(F9max_data$PC2_F9max ~ F9max_data$PC1_clim + F9max_data$PC2_clim +(1|(as.factor(F9max_data$Population))),data = F9max_data)
summary(F9max_PC2_1)

F9max_PC2_2 <- glm(F9max_data$PC2_F9max ~ F9max_data$PC1_clim + F9max_data$PC2_clim,data = F9max_data)
summary(F9max_PC2_2)

anova(F9max_PC2_1,F9max_PC2_2,test="Chisq")

#females max+
setwd("/Users/Венера/Dropbox/UoL/other/ClimateComponents/NASA/results/traits_combined/")
F9maxP_all_PC <- read.csv("/Users/Венера/Dropbox/UoL/other/ClimateComponents/NASA/results/Ewan/FmaxP_drosEU_PCcoords.csv")
F9maxP_data <- F9maxP_all_PC[,1:5]
colnames(F9maxP_data) <- c("Country", "Population", "Line", "PC1_F9maxP", "PC2_F9maxP")
climate_data <- read.csv("all_30y_PCA.csv")
clim <- climate_data[ -c(5) ]

new_cols <- c("Longitude", "Latitude", "Altitude","Location", "PC1_clim", "PC2_clim")
F9maxP_data[ , new_cols] <- NA
for (i in 1:nrow(F9maxP_data)){
  if (F9maxP_data[i,"Population"] == "RE"){ #Portugal
    F9maxP_data[i, 6:11] <- clim[1,]
  }
  if (F9maxP_data[i,"Population"] == "GI"){ #Spain
    F9maxP_data[i, 6:11] <- clim[2,]
  }
  if (F9maxP_data[i,"Population"] == "KA"){ #Denmark
    F9maxP_data[i, 6:11] <- clim[3,]
  }
  if (F9maxP_data[i,"Population"] == "MU"){ #Germany
    F9maxP_data[i, 6:11] <- clim[4,]
  }
  if (F9maxP_data[i,"Population"] == "MA"){ #Austria
    F9maxP_data[i, 6:11] <- clim[5,]
  }
  if (F9maxP_data[i,"Population"] == "AK"){ #Finland
    F9maxP_data[i, 6:11] <- clim[6,]
  }
  if (F9maxP_data[i,"Population"] == "UM"){ #Ukraine
    F9maxP_data[i, 6:11] <- clim[7,]
  }
  if (F9maxP_data[i,"Population"] == "YE"){ #Turkey
    F9maxP_data[i, 6:11] <- clim[8,]
  }
  if (F9maxP_data[i,"Population"] == "VA"){ #Russia
    F9maxP_data[i, 6:11] <- clim[9,]
  }
} 
write.csv(F9maxP_data,"F9maxP_30y_data.csv", row.names =F)


#phenotypes
#PC1
F9maxP_PC1_1 <- lmer(F9maxP_data$PC1_F9maxP ~F9maxP_data$PC1_clim+ F9maxP_data$PC1_clim + F9maxP_data$PC2_clim +(1|(as.factor(F9maxP_data$Population))),data = F9maxP_data)
summary(F9maxP_PC1_1)
F9maxP_PC1_2 <- lmer(F9maxP_data$PC1_F9maxP ~ F9maxP_data$PC1_clim +(1|(as.factor(F9maxP_data$Population))),data = F9maxP_data)
summary(F9maxP_PC1_2)

F9maxP_PC1_3 <- glm(F9maxP_data$PC1_F9maxP ~ F9maxP_data$PC1_clim + F9maxP_data$PC2_clim,data = F9maxP_data)
summary(F9maxP_PC1_3)
F9maxP_PC1_4 <- glm(F9maxP_data$PC1_F9maxP ~ F9maxP_data$PC1_clim,data = F9maxP_data)
summary(F9maxP_PC1_4)

anova(F9maxP_PC1_1,F9maxP_PC1_2,F9maxP_PC1_3,F9maxP_PC1_4, test="Chisq")

#PC2
F9maxP_PC2_1 <- lmer(F9maxP_data$PC2_F9maxP ~ F9maxP_data$PC1_clim + F9maxP_data$PC2_clim +(1|(as.factor(F9maxP_data$Population))),data = F9maxP_data)
summary(F9maxP_PC2_1)

F9maxP_PC2_2 <- glm(F9maxP_data$PC2_F9maxP ~ F9maxP_data$PC1_clim + F9maxP_data$PC2_clim,data = F9maxP_data)
summary(F9maxP_PC2_2)

anova(F9maxP_PC2_1,F9maxP_PC2_2,test="Chisq")
