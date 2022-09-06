#GDM DrosoEu 2022

# Data preparation
################################################################################
# CHUNK 1: Data preparation
################################################################################
# loading packages
library(gdm)
library(raster)
library(ggplot2)
library(readr)
library(dplyr)
library(vegan)
library(forcats)

setwd("D:/DrosoEu_GDM/Results/WA_L_F/")
output<- "D:/DrosoEu_GDM/Results/WA_L_F/"

#Data input, change line 21 for the trait you want to model
gdm_df <- readRDS("D:/DrosoEu_GDM/gdm_df_meta_population_estimates.rds")
df_meta_estimates<-gdm_df[["WA_L_F_meta_pop_compound_estimates"]]
df_meta_estimates<-df_meta_estimates[,c(1,4,11,12)]
col_order <- c("Population", "Longitude", "Latitude", "Mstar")
df_meta_estimates<- df_meta_estimates[, col_order]

#load environmental variables, if you don't have you need to 
#make one according GDM documentaion
envTabvar <- read_csv("D:/DrosoEu_GDM/envTabvar_All.csv", show_col_types = FALSE)

# remove lat and long to make dissimilarity matrix. 

df_meta_estimates_no_xy<- df_meta_estimates[4]
dissimilarity_matrix<-as.matrix(vegdist(decostand(df_meta_estimates_no_xy,
                                           MARGIN =2, "range"), "euclidean"))
#Prepare environmental variables to make formatsitepair df. 
site2 <- (df_meta_estimates[,c(1:3)])
envTab <- cbind(site2,envTabvar)
site3 <- site2[c(1)]
gdmDissim <- cbind(site3, dissimilarity_matrix)

# formatsitepair df 
gdmTab_rast <- formatsitepair(gdmDissim, 3, XColumn="Longitude", 
                              YColumn="Latitude", predData=envTab,
                              siteColumn="Population")
# make sure there are no NA values
sum(is.na(gdmTab_rast))
gdmTab_rast <- na.omit(gdmTab_rast)



###############################################################################
# CHUNK 2: Variable selection routine
# follows Ferrier et al. (2007) 
################################################################################
####General model
GDM1 <- gdm(gdmTab_rast, geo=F)
cat("Deviance full model = ", GDM1$explained, fill=T)
summary(GDM1)
length(GDM1$predictors) 
plot(GDM1, plot.layout=c(3,3))
# Select significant predictors
# loop to select significant predictors
modTest <- gdm.varImp(gdmTab_rast, geo=T, predSelect = T, nPerm = 10000, pValue=0.15)

#If using second modtest use below code to make selected variables df
#making data frame from modTest result, check the modTest to find out if more 
#variable were removed
#third line number must be updated for each run
modtest_pred_import <- as.data.frame(modTest [["Predictor Importance"]])
modtest_pred_import$names <- rownames(modtest_pred_import)
#change col or row num according predictor importance
modtest_pred_import <- modtest_pred_import[c(4,5)]
colnames(modtest_pred_import) <- c('predictors_Importance','predictor')
rownames(modtest_pred_import) <- 1:nrow(modtest_pred_import)

modtest_pred_pValue <- as.data.frame(modTest [["Predictor p-values"]])
modtest_pred_pValue $names <- rownames(modtest_pred_pValue )
##change col or row num according predictor importance
modtest_pred_pValue <- modtest_pred_pValue[c(4,5)]
colnames(modtest_pred_pValue) <- c('predictors_pValues','predictor')
rownames(modtest_pred_pValue) <- 1:nrow(modtest_pred_pValue )

modtest_model_assessment <- as.data.frame(modTest [["Model assessment"]])
sink(file='Model_assessment_in_backward _elimination.txt')
modtest_model_assessment
sink()

selected_predictor <- merge(modtest_pred_import, modtest_pred_pValue, by = 'predictor')
sum(is.na(selected_predictor))
selected_predictor <- na.omit(selected_predictor)
sink(file='Selected_predictors.txt')
selected_predictor
sink()

write.csv(selected_predictor , file=paste(output, "selected_predictor.csv" , sep = ""), 
          row.names = FALSE)


#Plot predictors
selected_predictor %>% 
  select(predictor, predictors_Importance) %>% 
  mutate(predictor = fct_reorder(predictor, predictors_Importance)) %>% 
  ggplot(aes(predictor, predictors_Importance, fill=predictor)) +
  geom_col()+
  scale_fill_grey(start=0.7, end=0.1) + 
  xlab("Environmental variables") +
  ylab("Predictor Importance") +
  theme(legend.position="none", axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  coord_flip()

# new formatsitepair with selected environmental variables 
selected_variables <-selected_predictor[1]
envTab2<-envTab[c(1,2,3)]
for(i in 1:nrow(selected_variables)){
  output = c(envTab[selected_variables$predictor[i]])
  envTab2 <- cbind(envTab2, output)
}
SelectedgdmTab_rast <- formatsitepair(gdmDissim, 3, XColumn="Longitude", 
                                      YColumn="Latitude", predData=envTab2,
                                      siteColumn="Population")
# make sure there are no NA values
sum(is.na(SelectedgdmTab_rast))
SelectedgdmTab_rast <- na.omit(SelectedgdmTab_rast)

# Second GDM model with selected environmental variables
cat("Modeling using reduced env dataset...", fill=T)
flush.console()
GDM2 <- gdm(SelectedgdmTab_rast, geo=F)
devDiff <- 100*(abs(GDM1$explained - GDM2$explained)/mean(c(GDM1$explained,
                                                            GDM2$explained)))


cat("Deviance full model = ", GDM1$explained, fill=T)
cat("Deviance reduced model = ", GDM2$explained, fill=T)
cat("Percent difference in deviance= ", devDiff, fill=T)
# if the difference in deviance is large than >0.5%,  run
# nperm permutations to test for significance of the GDM2 
nperm = 10000
if(devDiff>0.5){
  permDev1 <- NULL
  permDev2 <- NULL
  
  for (i in 1:nperm){ # loop for permutations
    sppdatPerm <- df_meta_estimates[sample(seq(1, nrow(df_meta_estimates)), nrow(df_meta_estimates), replace=F),]
    permgdmTab_rast1 <- formatsitepair(sppdatPerm, 1, siteColumn="Population", XColumn="Longitude", YColumn="Latitude",
                                       predData=envTab)
    permSelectedgdmTab_rast <- formatsitepair(sppdatPerm, 1, siteColumn="Population", XColumn="Longitude", YColumn="Latitude",
                                              predData=envTab2)
    permGDM1 <- gdm(permgdmTab_rast1, geo=T)
    permGDM2 <- gdm(permSelectedgdmTab_rast, geo=T)
    permDev1[i] <- permGDM1$explained
    permDev2[i] <- permGDM2$explained
    permDiffDev <- 100*(abs(permGDM1$explained - permGDM2$explained)/
                          mean(c(permGDM1$explained,permGDM2$explained)))
                                                                            
    cat((i/nperm)*100, "% COMPLETE", fill=TRUE)
    flush.console()
  }
  pval <- sum(devDiff < permDiffDev)/nperm
}
if(pval > 0.05){permgdmTab_rast1[-c(1:7)]<- permSelectedgdmTab_rast[-c(1:7)]}
flush.console()
sink(file='Sig_variables.txt')
cat("pval=", pval, fill=TRUE)
cat("Vars= ", names(permSelectedgdmTab_rast[-c(1:7)]), fill=T)
sink()
################################################################################
################################################################################
# CHUNK 3: Fit GDMs
# now fit three GDMs using: (1) both env and geo vars, (2) env only, (3) geo only
#follows Thomassen et. al (2009) & Baldassarre et al. (2013)
################################################################################
gdm_full <- gdm(SelectedgdmTab_rast, geo=T) # Selected Variables sing. vars
summary(gdm_full)
sink(file='gdm_full.txt')
summary(gdm_full)
sink()

gdm_Env <- gdm(SelectedgdmTab_rast, geo=F) # env only
summary(gdm_Env)
sink(file='gdm_selected_env.txt')
summary(gdm_Env)
sink()

gdm_XY <- gdm(gdmTab_rast[,1:6], geo=T) # geo only
summary(gdm_XY)
sink(file='gdm_xy.txt')
summary(gdm_XY)
sink()
#Make random environmental variables
site2 <- (df_meta_estimates[,c(1:3)])
nperm = 100000
devince_list<- list()
for (i in 1:nperm){
  envTabvar_random <- as.data.frame(matrix(runif
                                           (n=(nrow(selected_predictor)*nrow(df_meta_estimates)),
                                             min=0, max=100), ncol=nrow(selected_predictor)))
  envTab_random <- cbind(site2,envTabvar_random)
  # formatsitepair df 
  gdmTab_random <- formatsitepair(gdmDissim, 3, XColumn="Longitude", 
                                  YColumn="Latitude", predData=envTab_random,
                                  siteColumn="Population")
  #random environmental variables gdm
  gdm_Random <- gdm(gdmTab_random, geo=T)
  devince_random <- gdm_Random$explained
  devince_list[[i]] = devince_random
}
devince_df <- as.data.frame(do.call(rbind, devince_list))
colnames(devince_df) <- c("devince_each_model")
# extract deviance explained 
cat("Deviance full model = ", gdm_full$explained, fill=T)
cat("Deviance env model = ", gdm_Env$explained, fill=T)
cat("Deviance geo_distance_only model = ", gdm_XY$explained, fill=T)
cat("Deviance random model = ", mean(devince_df$devince_each_model), fill=T)
#save to text file
sink(file='deviance_All_models.txt')
cat("Deviance full model = ", gdm_full$explained, fill=T)
cat("Deviance env model = ", gdm_Env$explained, fill=T)
cat("Deviance geo_distance_only model = ", gdm_XY$explained, fill=T)
cat("Deviance random model = ", mean(devince_df$devince_each_model), fill=T)
sink()
################################################################################
################################################################################
# CHUNK 4: final model
#Final model for making maps and graphs 
################################################################################
gdm_final<- gdm(SelectedgdmTab_rast, geo=F)
#summary(gdm.1)
summary(gdm_final)
#Short summary
str(gdm_final)
#gdm plots
length(gdm_final$predictors)
plot(gdm_final, plot.layout=c(1,1))

# I-spline.
gdm_splineDat <- isplineExtract(gdm_final)
str(gdm_splineDat)
gdm_splineDat_y <- gdm_splineDat[["y"]]
gdm_splineDat_x <- gdm_splineDat[["x"]]

##plot GDM uncertainty 
plotUncertainty(SelectedgdmTab_rast, sampleSites=0.90, 
                bsIters=1000, geo=F, plot.layout=c(1,1), save=TRUE,
                fileName="gdm.plotUncertainy.csv")


#Predicting biological distances between sites
gdm.pred <- predict(gdm_final, SelectedgdmTab_rast)
head(gdm.pred)

#Plot predicated against observed
lm_formula <- y ~ x

predicted_plot <- ggplot() + 
  aes(x = SelectedgdmTab_rast$distance, y = gdm.pred) +
  geom_smooth(method = "lm", se = TRUE, color = "black", formula = lm_formula)+
  geom_point() + 
  ylab("Predicted dissimilarity") +
  theme_bw() +
  xlab("Observed dissimilarity")
predicted_plot

#Transforming predictors and visualizing biological patterns
# reordering environmental data to create table for transformation
# transform climate rasters & plot pattern
SelectedLayer <- list.files("D:/DrosoEu_GDM/Selected_ASC/",
                            pattern="*.asc", full.names=TRUE)
Selectedenvgrids <- c(SelectedLayer)
SelectedEnvstack <- stack(Selectedenvgrids)

rastTrans <- gdm.transform(gdm_final, SelectedEnvstack)

#Visualizing multi-dimensional biological patterns

rastDat <- na.omit(getValues(rastTrans))
#rastDat <- sampleRandom(rastTrans, 50000) # can use if rasters are large
pcaSamp <- prcomp(rastDat)
# note the use of the 'index' argument
pcaRast <- predict(rastTrans, pcaSamp, index=1:3)
# scale rasters
pcaRast[[1]] <- (pcaRast[[1]]-pcaRast[[1]]@data@min) /
  (pcaRast[[1]]@data@max-pcaRast[[1]]@data@min)*255
pcaRast[[2]] <- (pcaRast[[2]]-pcaRast[[2]]@data@min) /
  (pcaRast[[2]]@data@max-pcaRast[[2]]@data@min)*255
pcaRast[[3]] <- (pcaRast[[3]]-pcaRast[[3]]@data@min) /
  (pcaRast[[3]]@data@max-pcaRast[[3]]@data@min)*255
plotRGB(pcaRast, r=1, g=2, b=3, interpolate=T, stretch='lin')


#Export as EPS
setEPS()
postscript("TL_F.eps")
plotRGB(pcaRast, r=1, g=2, b=3, interpolate=T, stretch='lin')
dev.off()

