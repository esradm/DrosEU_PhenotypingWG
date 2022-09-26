setwd("C:/Users/Венера/Dropbox/UoL/other/ClimateComponents/NASA/all_param/")

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
all_climatology <- data.frame(matrix(NA, nrow = 9, ncol = 20))
col_names <- c("Longitude", "Latitude", "Altitude", "Population", "Country", 
               "TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
               "TS_MIN", "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS",
               "PRECTOTCORR", "PRECTOTCORR_SUM","ALLSKY_SFC_LW_DWN")

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
           "PRECTOTCORR", "PRECTOTCORR_SUM", "ALLSKY_SFC_LW_DWN"),
  temporal_api = "climatology"
)
all_climatology[1,6:20] <- Recarei_ag$ANN 

#Gimenells

Gimenells_ag <- get_power(
  community = "ag",
  lonlat = c(0.62, 41.618),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "PRECTOTCORR_SUM","ALLSKY_SFC_LW_DWN"),
  temporal_api = "climatology"
)
all_climatology[2,6:20] <- Gimenells_ag$ANN 

#Karensminde

Karensminde_ag <- get_power(
  community = "ag",
  lonlat = c(10.213, 55.945),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "PRECTOTCORR_SUM","ALLSKY_SFC_LW_DWN"),
  temporal_api = "climatology"
)
all_climatology[3,6:20] <- Karensminde_ag$ANN 

#Munich

Munich_ag <- get_power(
  community = "ag",
  lonlat = c(11.61, 48.18),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "PRECTOTCORR_SUM","ALLSKY_SFC_LW_DWN"),
  temporal_api = "climatology"
)
all_climatology[4,6:20] <- Munich_ag$ANN 

#Mauternbach

Mauternbach_ag <- get_power(
  community = "ag",
  lonlat = c(15.560, 48.375),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "PRECTOTCORR_SUM","ALLSKY_SFC_LW_DWN"),
  temporal_api = "climatology"
)
all_climatology[5,6:20] <- Mauternbach_ag$ANN 

#Akaa

Akaa_ag <- get_power(
  community = "ag",
  lonlat = c(23.520, 61.100),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "PRECTOTCORR_SUM","ALLSKY_SFC_LW_DWN"),
  temporal_api = "climatology"
)
all_climatology[6,6:20] <- Akaa_ag$ANN 

#Uman

Uman_ag <- get_power(
  community = "ag",
  lonlat = c(30.206, 48.753),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "PRECTOTCORR_SUM","ALLSKY_SFC_LW_DWN"),
  temporal_api = "climatology"
)
all_climatology[7,6:20] <- Uman_ag$ANN 

#Yesiloz
Yesiloz_ag <- get_power(
  community = "ag",
  lonlat = c(32.26, 40.231),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "PRECTOTCORR_SUM","ALLSKY_SFC_LW_DWN"),
  temporal_api = "climatology"
)
all_climatology[8,6:20] <- Yesiloz_ag$ANN 

#Valday
Valday_ag <- get_power(
  community = "ag",
  lonlat = c(33.244, 57.979),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "PRECTOTCORR_SUM","ALLSKY_SFC_LW_DWN"),
  temporal_api = "climatology"
)
all_climatology[9,6:20] <- Valday_ag$ANN 

write.csv(all_climatology,"all_param_climatology.csv", row.names = F)

#PCA
setwd("C:/Users/Венера/Dropbox/UoL/other/ClimateComponents/NASA/all_param/")
all_climatology <- read.csv("all_param_climatology.csv")
all_climatology <- all_climatology[,6:20]
rownames(all_climatology) <- c("Recarei", "Gimenells", "Karensminde", "Munich", "Mauternbach", "Akaa", "Uman", "Yesiloz", "Valday")

library("FactoMineR")
library("factoextra")
options(ggrepel.max.overlaps = Inf)
all_ann_clim_pca <- PCA(all_climatology, scale.unit = TRUE, graph = TRUE)
print(all_ann_clim_pca)

# matrix with eigenvalues
all_ann_clim_pca$eig
eig.val <- get_eigenvalue(all_ann_clim_pca)
eig.val #first 3 PCs >1

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

# Contributions of variables to PC3
fviz_contrib(all_ann_clim_pca, choice = "var", axes = 3, top = 19, title="Contribution of variables to PC3")

#pca result dataframe
all_climatology <- read.csv("all_param_climatology.csv")
pops_data <- all_climatology[,1:5]
PC1_3 <- all_ann_clim_pca$ind$coord[1:9, 1:3] #there are 3 PCs with eigvalues >1
all_ann_clim_results <- cbind(pops_data, PC1_3)
write.csv(all_ann_clim_results,"all_ann_clim_PCA.csv", row.names =T)

library(corrplot)
cor_df <- all_ann_clim_results[ -c(4:5) ]
corrplot(cor(cor_df), method = 'number') 
corrplot(cor(cor_df), type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)

#Phenotype analysis
setwd("C:/Users/Венера/Dropbox/UoL/other/ClimateComponents/NASA/all_param/phenotypes/")
all_clim_data <- read.csv("all_ann_clim_PCA.csv", row.names = 1)
#this file is from DrosEU drive
pheno <-  read.csv("all_models_line_meta_compound_random_coefs_wide.csv", header=T)
str(pheno)
nrow(pheno)

new_cols <- c("Country", "PC1", "PC2", "PC3")
pheno[ , new_cols] <- NA
for (i in 1:nrow(pheno)){
  if (pheno[i,"Population"] == "RE"){ #Portugal
    pheno[i, 35:38] <- all_clim_data[1,5:8]
  }
  if (pheno[i,"Population"] == "GI"){ #Spain
    pheno[i, 35:38] <- all_clim_data[2,5:8]
  }
  if (pheno[i,"Population"] == "KA"){ #Denmark
    pheno[i, 35:38] <- all_clim_data[3,5:8]
  }
  if (pheno[i,"Population"] == "MU"){ #Germany
    pheno[i, 35:38] <- all_clim_data[4,5:8]
  }
  if (pheno[i,"Population"] == "MA"){ #Austria
    pheno[i, 35:38] <- all_clim_data[5,5:8]
  }
  if (pheno[i,"Population"] == "AK"){ #Finland
    pheno[i, 35:38] <- all_clim_data[6,5:8]
  }
  if (pheno[i,"Population"] == "UM"){ #Ukraine
    pheno[i, 35:38] <- all_clim_data[7,5:8]
  }
  if (pheno[i,"Population"] == "YE"){ #Turkey
    pheno[i, 35:38] <- all_clim_data[8,5:8]
  }
  if (pheno[i,"Population"] == "VA"){ #Russia
    pheno[i, 35:38] <- all_clim_data[9,5:8]
  }
} 

#mixed-model with random effects
library(lme4)
library(lmerTest) 
pheno <- na.omit(pheno)

#CCRT_F
library(car)
CCRT_F_model1 <- lmer(pheno$CCRT_F ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
CCRT_F_model2 <- lmer(pheno$CCRT_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
CCRT_F_model3 <- lmer(pheno$CCRT_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
CCRT_F_model4 <- glm(pheno$CCRT_F ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
CCRT_F_model5 <- glm(pheno$CCRT_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
CCRT_F_model6 <- glm(pheno$CCRT_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(CCRT_F_model1, CCRT_F_model2,CCRT_F_model3, CCRT_F_model4,CCRT_F_model5,CCRT_F_model6,test="Chisq")
summary(CCRT_F_model2)
Anova(CCRT_F_model4)

#CCRT_M
CCRT_M_model1 <- lmer(pheno$CCRT_M ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
CCRT_M_model2 <- lmer(pheno$CCRT_M ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
CCRT_M_model3 <- lmer(pheno$CCRT_M ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
CCRT_M_model4 <- glm(pheno$CCRT_M ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
CCRT_M_model5 <- glm(pheno$CCRT_M ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
CCRT_M_model6 <- glm(pheno$CCRT_M ~ pheno$PC1 + pheno$Population,data = pheno)
anova(CCRT_M_model1, CCRT_M_model2,CCRT_M_model3, CCRT_M_model4,CCRT_M_model5,CCRT_M_model6,test="Chisq")
summary(CCRT_M_model3)
Anova(CCRT_M_model4)

#CSM_F
CSM_F_model1 <- lmer(pheno$CSM_F ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
CSM_F_model2 <- lmer(pheno$CSM_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
CSM_F_model3 <- lmer(pheno$CSM_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
CSM_F_model4 <- glm(pheno$CSM_F ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
CSM_F_model5 <- glm(pheno$CSM_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
CSM_F_model6 <- glm(pheno$CSM_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(CSM_F_model1, CSM_F_model2,CSM_F_model3, CSM_F_model4,CSM_F_model5,CSM_F_model6,test="Chisq")
summary(CSM_F_model3)
Anova(CSM_F_model4)

#CSM_M
CSM_M_model1 <- lmer(pheno$CSM_M ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
CSM_M_model2 <- lmer(pheno$CSM_M ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
CSM_M_model3 <- lmer(pheno$CSM_M ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
CSM_M_model4 <- glm(pheno$CSM_M ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
CSM_M_model5 <- glm(pheno$CSM_M ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
CSM_M_model6 <- glm(pheno$CSM_M ~ pheno$PC1 + pheno$Population,data = pheno)
anova(CSM_M_model1, CSM_M_model2,CSM_M_model3, CSM_M_model4,CSM_M_model5,CSM_M_model6,test="Chisq")
summary(CSM_M_model3)
Anova(CSM_M_model4)

#DT_A_F
DT_A_F_model1 <- lmer(pheno$DT_A_F ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
DT_A_F_model2 <- lmer(pheno$DT_A_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
DT_A_F_model3 <- lmer(pheno$DT_A_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
DT_A_F_model4 <- glm(pheno$DT_A_F ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
DT_A_F_model5 <- glm(pheno$DT_A_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
DT_A_F_model6 <- glm(pheno$DT_A_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(DT_A_F_model1, DT_A_F_model2,DT_A_F_model3, DT_A_F_model4,DT_A_F_model5,DT_A_F_model6,test="Chisq")
summary(DT_A_F_model1)
Anova(DT_A_F_model4)

#DT_A_M
DT_A_M_model1 <- lmer(pheno$DT_A_M ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
DT_A_M_model2 <- lmer(pheno$DT_A_M ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
DT_A_M_model3 <- lmer(pheno$DT_A_M ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
DT_A_M_model4 <- glm(pheno$DT_A_M ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
DT_A_M_model5 <- glm(pheno$DT_A_M ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
DT_A_M_model6 <- glm(pheno$DT_A_M ~ pheno$PC1 + pheno$Population,data = pheno)
anova(DT_A_M_model1, DT_A_M_model2,DT_A_M_model3, DT_A_M_model4,DT_A_M_model5,DT_A_M_model6,test="Chisq")
summary(DT_A_M_model1)
Anova(DT_A_M_model1)

#DT_P_NA
DT_P_NA_model1 <- lmer(pheno$DT_P_NA ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
DT_P_NA_model2 <- lmer(pheno$DT_P_NA ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
DT_P_NA_model3 <- lmer(pheno$DT_P_NA ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
DT_P_NA_model4 <- glm(pheno$DT_P_NA ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
DT_P_NA_model5 <- glm(pheno$DT_P_NA ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
DT_P_NA_model6 <- glm(pheno$DT_P_NA ~ pheno$PC1 + pheno$Population,data = pheno)
anova(DT_P_NA_model1, DT_P_NA_model2,DT_P_NA_model3, DT_P_NA_model4,DT_P_NA_model5,DT_P_NA_model6,test="Chisq")
summary(DT_P_NA_model1)
Anova(DT_P_NA_model1)

#Dia_F
Dia_F_model1 <- lmer(pheno$Dia_F ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
Dia_F_model2 <- lmer(pheno$Dia_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
Dia_F_model3 <- lmer(pheno$Dia_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
Dia_F_model4 <- glm(pheno$Dia_F ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
Dia_F_model5 <- glm(pheno$Dia_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
Dia_F_model6 <- glm(pheno$Dia_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(Dia_F_model1, Dia_F_model2,Dia_F_model3, Dia_F_model4,Dia_F_model5,Dia_F_model6,test="Chisq")
summary(Dia_F_model3)
Anova(Dia_F_model4)

#DW_F
DW_F_model1 <- lmer(pheno$DW_F ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
DW_F_model2 <- lmer(pheno$DW_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
DW_F_model3 <- lmer(pheno$DW_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
DW_F_model4 <- glm(pheno$DW_F ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
DW_F_model5 <- glm(pheno$DW_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
DW_F_model6 <- glm(pheno$DW_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(DW_F_model1, DW_F_model2,DW_F_model3, DW_F_model4,DW_F_model5,DW_F_model6,test="Chisq")
summary(DW_F_model3)
Anova(DW_F_model4)

#DW_M
DW_M_model1 <- lmer(pheno$DW_M ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
DW_M_model2 <- lmer(pheno$DW_M ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
DW_M_model3 <- lmer(pheno$DW_M ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
DW_M_model4 <- glm(pheno$DW_M ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
DW_M_model5 <- glm(pheno$DW_M ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
DW_M_model6 <- glm(pheno$DW_M ~ pheno$PC1 + pheno$Population,data = pheno)
anova(DW_M_model1, DW_M_model2,DW_M_model3, DW_M_model4,DW_M_model5,DW_M_model6,test="Chisq")
summary(DW_M_model2)
Anova(DW_M_model4)

#Fec_F
Fec_F_model1 <- lmer(pheno$Fec_F ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
Fec_F_model2 <- lmer(pheno$Fec_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
Fec_F_model3 <- lmer(pheno$Fec_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
Fec_F_model4 <- glm(pheno$Fec_F ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
Fec_F_model5 <- glm(pheno$Fec_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
Fec_F_model6 <- glm(pheno$Fec_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(Fec_F_model1, Fec_F_model2,Fec_F_model3, Fec_F_model4,Fec_F_model5,Fec_F_model6,test="Chisq")
summary(Fec_F_model2)
Anova(Fec_F_model4)

#HSM_F
HSM_F_model1 <- lmer(pheno$HSM_F ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
HSM_F_model2 <- lmer(pheno$HSM_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
HSM_F_model3 <- lmer(pheno$HSM_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
HSM_F_model4 <- glm(pheno$HSM_F ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
HSM_F_model5 <- glm(pheno$HSM_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
HSM_F_model6 <- glm(pheno$HSM_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(HSM_F_model1, HSM_F_model2,HSM_F_model3, HSM_F_model4,HSM_F_model5,HSM_F_model6,test="Chisq")
summary(HSM_F_model3)
Anova(HSM_F_model4)

#HSM_M
HSM_M_model1 <- lmer(pheno$HSM_M ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
HSM_M_model2 <- lmer(pheno$HSM_M ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
HSM_M_model3 <- lmer(pheno$HSM_M ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
HSM_M_model4 <- glm(pheno$HSM_M ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
HSM_M_model5 <- glm(pheno$HSM_M ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
HSM_M_model6 <- glm(pheno$HSM_M ~ pheno$PC1 + pheno$Population,data = pheno)
anova(HSM_M_model1, HSM_M_model2,HSM_M_model3, HSM_M_model4,HSM_M_model5,HSM_M_model6,test="Chisq")
summary(HSM_M_model3)
Anova(HSM_M_model4)

#LS_F
LS_F_model1 <- lmer(pheno$LS_F ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
LS_F_model2 <- lmer(pheno$LS_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
LS_F_model3 <- lmer(pheno$LS_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
LS_F_model4 <- glm(pheno$LS_F ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
LS_F_model5 <- glm(pheno$LS_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
LS_F_model6 <- glm(pheno$LS_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(LS_F_model1, LS_F_model2,LS_F_model3, LS_F_model4,LS_F_model5,LS_F_model6,test="Chisq")
summary(LS_F_model1)
Anova(LS_F_model4)

#LS_M
LS_M_model1 <- lmer(pheno$LS_M ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
LS_M_model2 <- lmer(pheno$LS_M ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
LS_M_model3 <- lmer(pheno$LS_M ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
LS_M_model4 <- glm(pheno$LS_M ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
LS_M_model5 <- glm(pheno$LS_M ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
LS_M_model6 <- glm(pheno$LS_M ~ pheno$PC1 + pheno$Population,data = pheno)
anova(LS_M_model1, LS_M_model2,LS_M_model3, LS_M_model4,LS_M_model5,LS_M_model6,test="Chisq")
summary(LS_M_model1)
Anova(LS_M_model4)

#LA_CircPhase_B
LA_CircPhase_B_model1 <- lmer(pheno$LA_CircPhase_B ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
LA_CircPhase_B_model2 <- lmer(pheno$LA_CircPhase_B ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
LA_CircPhase_B_model3 <- lmer(pheno$LA_CircPhase_B ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
LA_CircPhase_B_model4 <- glm(pheno$LA_CircPhase_B ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
LA_CircPhase_B_model5 <- glm(pheno$LA_CircPhase_B ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
LA_CircPhase_B_model6 <- glm(pheno$LA_CircPhase_B ~ pheno$PC1 + pheno$Population,data = pheno)
anova(LA_CircPhase_B_model1, LA_CircPhase_B_model2,LA_CircPhase_B_model3, LA_CircPhase_B_model4,LA_CircPhase_B_model5,LA_CircPhase_B_model6,test="Chisq")
summary(LA_CircPhase_B_model3)
Anova(LA_CircPhase_B_model4)


#LA_Activity_B
LA_Activity_B_model1 <- lmer(pheno$LA_Activity_B ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
LA_Activity_B_model2 <- lmer(pheno$LA_Activity_B ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
LA_Activity_B_model3 <- lmer(pheno$LA_Activity_B ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
LA_Activity_B_model4 <- glm(pheno$LA_Activity_B ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
LA_Activity_B_model5 <- glm(pheno$LA_Activity_B ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
LA_Activity_B_model6 <- glm(pheno$LA_Activity_B ~ pheno$PC1 + pheno$Population,data = pheno)
anova(LA_Activity_B_model1, LA_Activity_B_model2,LA_Activity_B_model3, LA_Activity_B_model4,LA_Activity_B_model5,LA_Activity_B_model6,test="Chisq")
summary(LA_Activity_B_model3)
Anova(LA_Activity_B_model4)

#LA_NDlog2_B
LA_NDlog2_B_model1 <- lmer(pheno$LA_NDlog2_B ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
LA_NDlog2_B_model2 <- lmer(pheno$LA_NDlog2_B ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
LA_NDlog2_B_model3 <- lmer(pheno$LA_NDlog2_B ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
LA_NDlog2_B_model4 <- glm(pheno$LA_NDlog2_B ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
LA_NDlog2_B_model5 <- glm(pheno$LA_NDlog2_B ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
LA_NDlog2_B_model6 <- glm(pheno$LA_NDlog2_B ~ pheno$PC1 + pheno$Population,data = pheno)
anova(LA_NDlog2_B_model1, LA_NDlog2_B_model2,LA_NDlog2_B_model3, LA_NDlog2_B_model4,LA_NDlog2_B_model5,LA_NDlog2_B_model6,test="Chisq")
summary(LA_NDlog2_B_model3)
Anova(LA_NDlog2_B_model4)

#LA_Period_B
LA_Period_B_model1 <- lmer(pheno$LA_Period_B ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
LA_Period_B_model2 <- lmer(pheno$LA_Period_B ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
LA_Period_B_model3 <- lmer(pheno$LA_Period_B ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
LA_Period_B_model4 <- glm(pheno$LA_Period_B ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
LA_Period_B_model5 <- glm(pheno$LA_Period_B ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
LA_Period_B_model6 <- glm(pheno$LA_Period_B ~ pheno$PC1 + pheno$Population,data = pheno)
anova(LA_Period_B_model1, LA_Period_B_model2,LA_Period_B_model3, LA_Period_B_model4,LA_Period_B_model5,LA_Period_B_model6,test="Chisq")
summary(LA_Period_B_model3)
Anova(LA_Period_B_model4)

#Pgm_T4_F
Pgm_T4_F_model1 <- lmer(pheno$Pgm_T4_F ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
Pgm_T4_F_model2 <- lmer(pheno$Pgm_T4_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
Pgm_T4_F_model3 <- lmer(pheno$Pgm_T4_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
Pgm_T4_F_model4 <- glm(pheno$Pgm_T4_F ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
Pgm_T4_F_model5 <- glm(pheno$Pgm_T4_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
Pgm_T4_F_model6 <- glm(pheno$Pgm_T4_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(Pgm_T4_F_model1, Pgm_T4_F_model2,Pgm_T4_F_model3, Pgm_T4_F_model4,Pgm_T4_F_model5,Pgm_T4_F_model6,test="Chisq")
summary(Pgm_T4_F_model3)
Anova(Pgm_T4_F_model4)

#Pgm_T5_F
Pgm_T5_F_model1 <- lmer(pheno$Pgm_T5_F ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
Pgm_T5_F_model2 <- lmer(pheno$Pgm_T5_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
Pgm_T5_F_model3 <- lmer(pheno$Pgm_T5_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
Pgm_T5_F_model4 <- glm(pheno$Pgm_T5_F ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
Pgm_T5_F_model5 <- glm(pheno$Pgm_T5_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
Pgm_T5_F_model6 <- glm(pheno$Pgm_T5_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(Pgm_T5_F_model1, Pgm_T5_F_model2,Pgm_T5_F_model3, Pgm_T5_F_model4,Pgm_T5_F_model5,Pgm_T5_F_model6,test="Chisq")
summary(Pgm_T5_F_model3)
Anova(Pgm_T5_F_model4)

#Pgm_T6_F
Pgm_T6_F_model1 <- lmer(pheno$Pgm_T6_F ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
Pgm_T6_F_model2 <- lmer(pheno$Pgm_T6_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
Pgm_T6_F_model3 <- lmer(pheno$Pgm_T6_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
Pgm_T6_F_model4 <- glm(pheno$Pgm_T6_F ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
Pgm_T6_F_model5 <- glm(pheno$Pgm_T6_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
Pgm_T6_F_model6 <- glm(pheno$Pgm_T6_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(Pgm_T6_F_model1, Pgm_T6_F_model2,Pgm_T6_F_model3, Pgm_T6_F_model4,Pgm_T6_F_model5,Pgm_T6_F_model6,test="Chisq")
summary(Pgm_T6_F_model3)
Anova(Pgm_T6_F_model4)

#Pgm_Total_F
Pgm_Total_F_model1 <- lmer(pheno$Pgm_Total_F ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
Pgm_Total_F_model2 <- lmer(pheno$Pgm_Total_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
Pgm_Total_F_model3 <- lmer(pheno$Pgm_Total_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
Pgm_Total_F_model4 <- glm(pheno$Pgm_Total_F ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
Pgm_Total_F_model5 <- glm(pheno$Pgm_Total_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
Pgm_Total_F_model6 <- glm(pheno$Pgm_Total_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(Pgm_Total_F_model1, Pgm_Total_F_model2,Pgm_Total_F_model3, Pgm_Total_F_model4,Pgm_Total_F_model5,Pgm_Total_F_model6,test="Chisq")
summary(Pgm_Total_F_model3)
Anova(Pgm_Total_F_model4)

#SR_F
SR_F_model1 <- lmer(pheno$SR_F ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
SR_F_model2 <- lmer(pheno$SR_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
SR_F_model3 <- lmer(pheno$SR_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
SR_F_model4 <- glm(pheno$SR_F ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
SR_F_model5 <- glm(pheno$SR_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
SR_F_model6 <- glm(pheno$SR_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(SR_F_model1, SR_F_model2,SR_F_model3, SR_F_model4,SR_F_model5,SR_F_model6,test="Chisq")
summary(SR_F_model2)
Anova(SR_F_model4)

#SR_M
SR_M_model1 <- lmer(pheno$SR_M ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
SR_M_model2 <- lmer(pheno$SR_M ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
SR_M_model3 <- lmer(pheno$SR_M ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
SR_M_model4 <- glm(pheno$SR_M ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
SR_M_model5 <- glm(pheno$SR_M ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
SR_M_model6 <- glm(pheno$SR_M ~ pheno$PC1 + pheno$Population,data = pheno)
anova(SR_M_model1, SR_M_model2,SR_M_model3, SR_M_model4,SR_M_model5,SR_M_model6,test="Chisq")
summary(SR_M_model2)
Anova(SR_M_model4)

#TL_F
TL_F_model1 <- lmer(pheno$TL_F ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
TL_F_model2 <- lmer(pheno$TL_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
TL_F_model3 <- lmer(pheno$TL_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
TL_F_model4 <- glm(pheno$TL_F ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
TL_F_model5 <- glm(pheno$TL_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
TL_F_model6 <- glm(pheno$TL_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(TL_F_model1, TL_F_model2,TL_F_model3, TL_F_model4,TL_F_model5,TL_F_model6,test="Chisq")
summary(TL_F_model1)
Anova(TL_F_model4)

#TL_M
TL_M_model1 <- lmer(pheno$TL_M ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
TL_M_model2 <- lmer(pheno$TL_M ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
TL_M_model3 <- lmer(pheno$TL_M ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
TL_M_model4 <- glm(pheno$TL_M ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
TL_M_model5 <- glm(pheno$TL_M ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
TL_M_model6 <- glm(pheno$TL_M ~ pheno$PC1 + pheno$Population,data = pheno)
anova(TL_M_model1, TL_M_model2,TL_M_model3, TL_M_model4,TL_M_model5,TL_M_model6,test="Chisq")
summary(TL_M_model1)
Anova(TL_M_model4)

#Via_NA
Via_NA_model1 <- lmer(pheno$Via_NA ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
Via_NA_model2 <- lmer(pheno$Via_NA ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
Via_NA_model3 <- lmer(pheno$Via_NA ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
Via_NA_model4 <- glm(pheno$Via_NA ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
Via_NA_model5 <- glm(pheno$Via_NA ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
Via_NA_model6 <- glm(pheno$Via_NA ~ pheno$PC1 + pheno$Population,data = pheno)
anova(Via_NA_model1, Via_NA_model2,Via_NA_model3, Via_NA_model4,Via_NA_model5,Via_NA_model6,test="Chisq")
summary(Via_NA_model2)
Anova(Via_NA_model4)

#WA_L_F
WA_L_F_model1 <- lmer(pheno$WA_L_F ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
WA_L_F_model2 <- lmer(pheno$WA_L_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
WA_L_F_model3 <- lmer(pheno$WA_L_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
WA_L_F_model4 <- glm(pheno$WA_L_F ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
WA_L_F_model5 <- glm(pheno$WA_L_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
WA_L_F_model6 <- glm(pheno$WA_L_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(WA_L_F_model1, WA_L_F_model2,WA_L_F_model3, WA_L_F_model4,WA_L_F_model5,WA_L_F_model6,test="Chisq")
summary(WA_L_F_model1)
Anova(WA_L_F_model4)


#WA_L_M
WA_L_M_model1 <- lmer(pheno$WA_L_M ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
WA_L_M_model2 <- lmer(pheno$WA_L_M ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
WA_L_M_model3 <- lmer(pheno$WA_L_M ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
WA_L_M_model4 <- glm(pheno$WA_L_M ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
WA_L_M_model5 <- glm(pheno$WA_L_M ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
WA_L_M_model6 <- glm(pheno$WA_L_M ~ pheno$PC1 + pheno$Population,data = pheno)
anova(WA_L_M_model1, WA_L_M_model2,WA_L_M_model3, WA_L_M_model4,WA_L_M_model5,WA_L_M_model6,test="Chisq")
summary(WA_L_M_model1)
Anova(WA_L_M_model4)


#WA_R_F
WA_R_F_model1 <- lmer(pheno$WA_R_F ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
WA_R_F_model2 <- lmer(pheno$WA_R_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
WA_R_F_model3 <- lmer(pheno$WA_R_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
WA_R_F_model4 <- glm(pheno$WA_R_F ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
WA_R_F_model5 <- glm(pheno$WA_R_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
WA_R_F_model6 <- glm(pheno$WA_R_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(WA_R_F_model1, WA_R_F_model2,WA_R_F_model3, WA_R_F_model4,WA_R_F_model5,WA_R_F_model6,test="Chisq")
summary(WA_R_F_model1)
Anova(WA_R_F_model4)


#WA_R_M
WA_R_M_model1 <- lmer(pheno$WA_R_M ~ pheno$PC1 + pheno$PC2 +pheno$PC3 +(1|(as.factor(pheno$Population))),data = pheno)
WA_R_M_model2 <- lmer(pheno$WA_R_M ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
WA_R_M_model3 <- lmer(pheno$WA_R_M ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
WA_R_M_model4 <- glm(pheno$WA_R_M ~ pheno$PC1 + pheno$PC2 + pheno$PC3 + pheno$Population,data = pheno)
WA_R_M_model5 <- glm(pheno$WA_R_M ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
WA_R_M_model6 <- glm(pheno$WA_R_M ~ pheno$PC1 + pheno$Population,data = pheno)
anova(WA_R_M_model1, WA_R_M_model2,WA_R_M_model3, WA_R_M_model4,WA_R_M_model5,WA_R_M_model6,test="Chisq")
summary(WA_R_M_model1)
Anova(WA_R_M_model4)
