setwd("C:/Users/Венера/Dropbox/UoL/other/ClimateComponents/NASA/all_param/")
#30days prior the collection date

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
d30_df <- data.frame(matrix(NA, nrow = 9, ncol = 19))
col_names <- c("Longitude", "Latitude", "Altitude", "Population", "Country", 
               "TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
               "TS_MIN", "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS",
               "PRECTOTCORR", "ALLSKY_SFC_LW_DWN")

colnames(d30_df) <-c(col_names)

d30_df$Longitude <- c(-8.41, 0.62, 10.213, 11.61, 15.56, 23.52, 30.206, 32.26, 33.244)
d30_df$Latitude <- c(41.15, 41.618, 55.945, 48.18, 48.375, 61.1, 48.753, 40.231, 57.979)
d30_df$Altitude <- c(175, 173, 15, 520, 572, 88, 214, 680, 217)
d30_df$Population <- c("Recarei", "Gimenells", "Karensminde", "Munich", "Mauternbach", "Akaa", "Uman", "Yesiloz", "Valday")
d30_df$Country <- c("Portugal", "Spain", "Denmark", "Germany", "Austria", "Finland", "Ukaraine", "Turkey", "Russia")

#Recarei

Recarei_ag <- get_power(
  community = "ag",
  lonlat = c(-8.410, 41.150),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "ALLSKY_SFC_LW_DWN"),
  dates = c("2018-09-06", "2018-10-05"),
  temporal_api = "daily"
)
d30_df$TS[1] <- mean(Recarei_ag$TS)
d30_df$T2M[1] <- mean(Recarei_ag$T2M)
d30_df$QV2M[1]<- mean(Recarei_ag$QV2M)
d30_df$RH2M[1]<- mean(Recarei_ag$RH2M)
d30_df$T2MDEW[1]<- mean(Recarei_ag$T2MDEW)
d30_df$T2MWET[1]<- mean(Recarei_ag$T2MWET)

a <- -50
for (i in 1:nrow(Recarei_ag)){
  if (Recarei_ag$TS_MAX[i] > a){
    a <- Recarei_ag$TS_MAX[i]
  }
}
d30_df$TS_MAX[1] <- a

b <- 50 
for (i in 1:nrow(Recarei_ag)){
  if (Recarei_ag$TS_MIN[i] < b){
    b <- Recarei_ag$TS_MIN[i]
  }
}
d30_df$TS_MIN[1] <- b

c <- -50
for (i in 1:nrow(Recarei_ag)){
  if (Recarei_ag$T2M_MAX[i] > c){
    c <- Recarei_ag$T2M_MAX[i]
  }
}
d30_df$T2M_MAX[1] <- c

d <- 50 
for (i in 1:nrow(Recarei_ag)){
  if (Recarei_ag$T2M_MIN[i] < d){
    d <- Recarei_ag$T2M_MIN[i]
  }
}
d30_df$T2M_MIN[1] <- d
d30_df$T2M_RANGE[1] <- mean(Recarei_ag$T2M_RANGE)
d30_df$FROST_DAYS[1] <- sum(Recarei_ag$FROST_DAYS)
d30_df$PRECTOTCORR[1] <- mean(Recarei_ag$PRECTOTCORR)
d30_df$ALLSKY_SFC_LW_DWN[1] <- mean(Recarei_ag$ALLSKY_SFC_LW_DWN)

rm(a, b, c, d, i)

#Gimenells

Gimenells_ag <- get_power(
  community = "ag",
  lonlat = c(0.62, 41.618),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "ALLSKY_SFC_LW_DWN"),
  dates = c("2018-07-27", "2018-08-25"),
  temporal_api = "daily"
)
d30_df$TS[2] <- mean(Gimenells_ag$TS)
d30_df$T2M[2] <- mean(Gimenells_ag$T2M)
d30_df$QV2M[2]<- mean(Gimenells_ag$QV2M)
d30_df$RH2M[2]<- mean(Gimenells_ag$RH2M)
d30_df$T2MDEW[2]<- mean(Gimenells_ag$T2MDEW)
d30_df$T2MWET[2]<- mean(Gimenells_ag$T2MWET)

a <- -50
for (i in 1:nrow(Gimenells_ag)){
  if (Gimenells_ag$TS_MAX[i] > a){
    a <- Gimenells_ag$TS_MAX[i]
  }
}
d30_df$TS_MAX[2] <- a

b <- 50 
for (i in 1:nrow(Gimenells_ag)){
  if (Gimenells_ag$TS_MIN[i] < b){
    b <- Gimenells_ag$TS_MIN[i]
  }
}
d30_df$TS_MIN[2] <- b

c <- -50
for (i in 1:nrow(Gimenells_ag)){
  if (Gimenells_ag$T2M_MAX[i] > c){
    c <- Gimenells_ag$T2M_MAX[i]
  }
}
d30_df$T2M_MAX[2] <- c

d <- 50 
for (i in 1:nrow(Gimenells_ag)){
  if (Gimenells_ag$T2M_MIN[i] < d){
    d <- Gimenells_ag$T2M_MIN[i]
  }
}
d30_df$T2M_MIN[2] <- d
d30_df$T2M_RANGE[2] <- mean(Gimenells_ag$T2M_RANGE)
d30_df$FROST_DAYS[2] <- sum(Gimenells_ag$FROST_DAYS)
d30_df$PRECTOTCORR[2] <- mean(Gimenells_ag$PRECTOTCORR)
d30_df$ALLSKY_SFC_LW_DWN[2] <- mean(Gimenells_ag$ALLSKY_SFC_LW_DWN)

rm(a, b, c, d, i)


#Karensminde

Karensminde_ag <- get_power(
  community = "ag",
  lonlat = c(10.213, 55.945),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "ALLSKY_SFC_LW_DWN"),
  dates = c("2018-09-06", "2018-10-05"),
  temporal_api = "daily"
)
d30_df$TS[3] <- mean(Karensminde_ag$TS)
d30_df$T2M[3] <- mean(Karensminde_ag$T2M)
d30_df$QV2M[3]<- mean(Karensminde_ag$QV2M)
d30_df$RH2M[3]<- mean(Karensminde_ag$RH2M)
d30_df$T2MDEW[3]<- mean(Karensminde_ag$T2MDEW)
d30_df$T2MWET[3]<- mean(Karensminde_ag$T2MWET)

a <- -50
for (i in 1:nrow(Karensminde_ag)){
  if (Karensminde_ag$TS_MAX[i] > a){
    a <- Karensminde_ag$TS_MAX[i]
  }
}
d30_df$TS_MAX[3] <- a

b <- 50 
for (i in 1:nrow(Karensminde_ag)){
  if (Karensminde_ag$TS_MIN[i] < b){
    b <- Karensminde_ag$TS_MIN[i]
  }
}
d30_df$TS_MIN[3] <- b

c <- -50
for (i in 1:nrow(Karensminde_ag)){
  if (Karensminde_ag$T2M_MAX[i] > c){
    c <- Karensminde_ag$T2M_MAX[i]
  }
}
d30_df$T2M_MAX[3] <- c

d <- 50 
for (i in 1:nrow(Karensminde_ag)){
  if (Karensminde_ag$T2M_MIN[i] < d){
    d <- Karensminde_ag$T2M_MIN[i]
  }
}
d30_df$T2M_MIN[3] <- d
d30_df$T2M_RANGE[3] <- mean(Karensminde_ag$T2M_RANGE)
d30_df$FROST_DAYS[3] <- sum(Karensminde_ag$FROST_DAYS)
d30_df$PRECTOTCORR[3] <- mean(Karensminde_ag$PRECTOTCORR)
d30_df$ALLSKY_SFC_LW_DWN[3] <- mean(Karensminde_ag$ALLSKY_SFC_LW_DWN)

rm(a, b, c, d, i)


#Munich

Munich_ag <- get_power(
  community = "ag",
  lonlat = c(11.61, 48.18),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "ALLSKY_SFC_LW_DWN"),
  dates = c("2018-09-06", "2018-10-05"),
  temporal_api = "daily"
)
d30_df$TS[4] <- mean(Munich_ag$TS)
d30_df$T2M[4] <- mean(Munich_ag$T2M)
d30_df$QV2M[4]<- mean(Munich_ag$QV2M)
d30_df$RH2M[4]<- mean(Munich_ag$RH2M)
d30_df$T2MDEW[4]<- mean(Munich_ag$T2MDEW)
d30_df$T2MWET[4]<- mean(Munich_ag$T2MWET)

a <- -50
for (i in 1:nrow(Munich_ag)){
  if (Munich_ag$TS_MAX[i] > a){
    a <- Munich_ag$TS_MAX[i]
  }
}
d30_df$TS_MAX[4] <- a

b <- 50 
for (i in 1:nrow(Munich_ag)){
  if (Munich_ag$TS_MIN[i] < b){
    b <- Munich_ag$TS_MIN[i]
  }
}
d30_df$TS_MIN[4] <- b

c <- -50
for (i in 1:nrow(Munich_ag)){
  if (Munich_ag$T2M_MAX[i] > c){
    c <- Munich_ag$T2M_MAX[i]
  }
}
d30_df$T2M_MAX[4] <- c

d <- 50 
for (i in 1:nrow(Munich_ag)){
  if (Munich_ag$T2M_MIN[i] < d){
    d <- Munich_ag$T2M_MIN[i]
  }
}
d30_df$T2M_MIN[4] <- d
d30_df$T2M_RANGE[4] <- mean(Munich_ag$T2M_RANGE)
d30_df$FROST_DAYS[4] <- sum(Munich_ag$FROST_DAYS)
d30_df$PRECTOTCORR[4] <- mean(Munich_ag$PRECTOTCORR)
d30_df$ALLSKY_SFC_LW_DWN[4] <- mean(Munich_ag$ALLSKY_SFC_LW_DWN)

rm(a, b, c, d, i)


#Mauternbach

Mauternbach_ag <- get_power(
  community = "ag",
  lonlat = c(15.560, 48.375),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "ALLSKY_SFC_LW_DWN"),
  dates = c("2018-09-06", "2018-10-05"),
  temporal_api = "daily"
)
d30_df$TS[5] <- mean(Mauternbach_ag$TS)
d30_df$T2M[5] <- mean(Mauternbach_ag$T2M)
d30_df$QV2M[5]<- mean(Mauternbach_ag$QV2M)
d30_df$RH2M[5]<- mean(Mauternbach_ag$RH2M)
d30_df$T2MDEW[5]<- mean(Mauternbach_ag$T2MDEW)
d30_df$T2MWET[5]<- mean(Mauternbach_ag$T2MWET)

a <- -50
for (i in 1:nrow(Mauternbach_ag)){
  if (Mauternbach_ag$TS_MAX[i] > a){
    a <- Mauternbach_ag$TS_MAX[i]
  }
}
d30_df$TS_MAX[5] <- a

b <- 50 
for (i in 1:nrow(Mauternbach_ag)){
  if (Mauternbach_ag$TS_MIN[i] < b){
    b <- Mauternbach_ag$TS_MIN[i]
  }
}
d30_df$TS_MIN[5] <- b

c <- -50
for (i in 1:nrow(Mauternbach_ag)){
  if (Mauternbach_ag$T2M_MAX[i] > c){
    c <- Mauternbach_ag$T2M_MAX[i]
  }
}
d30_df$T2M_MAX[5] <- c

d <- 50 
for (i in 1:nrow(Mauternbach_ag)){
  if (Mauternbach_ag$T2M_MIN[i] < d){
    d <- Mauternbach_ag$T2M_MIN[i]
  }
}
d30_df$T2M_MIN[5] <- d
d30_df$T2M_RANGE[5] <- mean(Mauternbach_ag$T2M_RANGE)
d30_df$FROST_DAYS[5] <- sum(Mauternbach_ag$FROST_DAYS)
d30_df$PRECTOTCORR[5] <- mean(Mauternbach_ag$PRECTOTCORR)
d30_df$ALLSKY_SFC_LW_DWN[5] <- mean(Mauternbach_ag$ALLSKY_SFC_LW_DWN)

rm(a, b, c, d, i)


#Akaa
Akaa_ag <- get_power(
  community = "ag",
  lonlat = c(23.520, 61.100),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "ALLSKY_SFC_LW_DWN"),
  dates = c("2018-06-21", "2018-07-20"),
  temporal_api = "daily"
)
d30_df$TS[6] <- mean(Akaa_ag$TS)
d30_df$T2M[6] <- mean(Akaa_ag$T2M)
d30_df$QV2M[6]<- mean(Akaa_ag$QV2M)
d30_df$RH2M[6]<- mean(Akaa_ag$RH2M)
d30_df$T2MDEW[6]<- mean(Akaa_ag$T2MDEW)
d30_df$T2MWET[6]<- mean(Akaa_ag$T2MWET)

a <- -50
for (i in 1:nrow(Akaa_ag)){
  if (Akaa_ag$TS_MAX[i] > a){
    a <- Akaa_ag$TS_MAX[i]
  }
}
d30_df$TS_MAX[6] <- a

b <- 50 
for (i in 1:nrow(Akaa_ag)){
  if (Akaa_ag$TS_MIN[i] < b){
    b <- Akaa_ag$TS_MIN[i]
  }
}
d30_df$TS_MIN[6] <- b

c <- -50
for (i in 1:nrow(Akaa_ag)){
  if (Akaa_ag$T2M_MAX[i] > c){
    c <- Akaa_ag$T2M_MAX[i]
  }
}
d30_df$T2M_MAX[6] <- c

d <- 50 
for (i in 1:nrow(Akaa_ag)){
  if (Akaa_ag$T2M_MIN[i] < d){
    d <- Akaa_ag$T2M_MIN[i]
  }
}
d30_df$T2M_MIN[6] <- d
d30_df$T2M_RANGE[6] <- mean(Akaa_ag$T2M_RANGE)
d30_df$FROST_DAYS[6] <- sum(Akaa_ag$FROST_DAYS)
d30_df$PRECTOTCORR[6] <- mean(Akaa_ag$PRECTOTCORR)
d30_df$ALLSKY_SFC_LW_DWN[6] <- mean(Akaa_ag$ALLSKY_SFC_LW_DWN)

rm(a,b,c,d,i)


#Uman
Uman_ag <- get_power(
  community = "ag",
  lonlat = c(30.206, 48.753),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "ALLSKY_SFC_LW_DWN"),
  dates = c("2018-07-20", "2018-08-18"),
  temporal_api = "daily"
)
d30_df$TS[7] <- mean(Uman_ag$TS)
d30_df$T2M[7] <- mean(Uman_ag$T2M)
d30_df$QV2M[7]<- mean(Uman_ag$QV2M)
d30_df$RH2M[7]<- mean(Uman_ag$RH2M)
d30_df$T2MDEW[7]<- mean(Uman_ag$T2MDEW)
d30_df$T2MWET[7]<- mean(Uman_ag$T2MWET)

a <- -50
for (i in 1:nrow(Uman_ag)){
  if (Uman_ag$TS_MAX[i] > a){
    a <- Uman_ag$TS_MAX[i]
  }
}
d30_df$TS_MAX[7] <- a

b <- 50 
for (i in 1:nrow(Uman_ag)){
  if (Uman_ag$TS_MIN[i] < b){
    b <- Uman_ag$TS_MIN[i]
  }
}
d30_df$TS_MIN[7] <- b

c <- -50
for (i in 1:nrow(Uman_ag)){
  if (Uman_ag$T2M_MAX[i] > c){
    c <- Uman_ag$T2M_MAX[i]
  }
}
d30_df$T2M_MAX[7] <- c

d <- 50 
for (i in 1:nrow(Uman_ag)){
  if (Uman_ag$T2M_MIN[i] < d){
    d <- Uman_ag$T2M_MIN[i]
  }
}
d30_df$T2M_MIN[7] <- d
d30_df$T2M_RANGE[7] <- mean(Uman_ag$T2M_RANGE)
d30_df$FROST_DAYS[7] <- sum(Uman_ag$FROST_DAYS)
d30_df$PRECTOTCORR[7] <- mean(Uman_ag$PRECTOTCORR)
d30_df$ALLSKY_SFC_LW_DWN[7] <- mean(Uman_ag$ALLSKY_SFC_LW_DWN)

rm(a,b,c,d,i)

#Yesiloz

Yesiloz_ag <- get_power(
  community = "ag",
  lonlat = c(32.26, 40.231),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "ALLSKY_SFC_LW_DWN"),
  dates = c("2018-08-29", "2018-09-27"),
  temporal_api = "daily"
)
d30_df$TS[8] <- mean(Yesiloz_ag$TS)
d30_df$T2M[8] <- mean(Yesiloz_ag$T2M)
d30_df$QV2M[8]<- mean(Yesiloz_ag$QV2M)
d30_df$RH2M[8]<- mean(Yesiloz_ag$RH2M)
d30_df$T2MDEW[8]<- mean(Yesiloz_ag$T2MDEW)
d30_df$T2MWET[8]<- mean(Yesiloz_ag$T2MWET)

a <- -50
for (i in 1:nrow(Yesiloz_ag)){
  if (Yesiloz_ag$TS_MAX[i] > a){
    a <- Yesiloz_ag$TS_MAX[i]
  }
}
d30_df$TS_MAX[8] <- a

b <- 50 
for (i in 1:nrow(Yesiloz_ag)){
  if (Yesiloz_ag$TS_MIN[i] < b){
    b <- Yesiloz_ag$TS_MIN[i]
  }
}
d30_df$TS_MIN[8] <- b

c <- -50
for (i in 1:nrow(Yesiloz_ag)){
  if (Yesiloz_ag$T2M_MAX[i] > c){
    c <- Yesiloz_ag$T2M_MAX[i]
  }
}
d30_df$T2M_MAX[8] <- c

d <- 50 
for (i in 1:nrow(Yesiloz_ag)){
  if (Yesiloz_ag$T2M_MIN[i] < d){
    d <- Yesiloz_ag$T2M_MIN[i]
  }
}
d30_df$T2M_MIN[8] <- d
d30_df$T2M_RANGE[8] <- mean(Yesiloz_ag$T2M_RANGE)
d30_df$FROST_DAYS[8] <- sum(Yesiloz_ag$FROST_DAYS)
d30_df$PRECTOTCORR[8] <- mean(Yesiloz_ag$PRECTOTCORR)
d30_df$ALLSKY_SFC_LW_DWN[8] <- mean(Yesiloz_ag$ALLSKY_SFC_LW_DWN)

rm(a,b,c,d,i)

#Valday

Valday_ag <- get_power(
  community = "ag",
  lonlat = c(33.244, 57.979),
  pars = c("TS", "T2M", "QV2M", "RH2M", "T2MDEW", "T2MWET", "TS_MAX", 
           "TS_MIN",   "T2M_MAX", "T2M_MIN", "T2M_RANGE", "FROST_DAYS", 
           "PRECTOTCORR", "ALLSKY_SFC_LW_DWN"),
  dates = c("2018-07-22", "2018-08-20"),
  temporal_api = "daily"
)
d30_df$TS[9] <- mean(Valday_ag$TS)
d30_df$T2M[9] <- mean(Valday_ag$T2M)
d30_df$QV2M[9]<- mean(Valday_ag$QV2M)
d30_df$RH2M[9]<- mean(Valday_ag$RH2M)
d30_df$T2MDEW[9]<- mean(Valday_ag$T2MDEW)
d30_df$T2MWET[9]<- mean(Valday_ag$T2MWET)

a <- -50
for (i in 1:nrow(Valday_ag)){
  if (Valday_ag$TS_MAX[i] > a){
    a <- Valday_ag$TS_MAX[i]
  }
}
d30_df$TS_MAX[9] <- a

b <- 50 
for (i in 1:nrow(Valday_ag)){
  if (Valday_ag$TS_MIN[i] < b){
    b <- Valday_ag$TS_MIN[i]
  }
}
d30_df$TS_MIN[9] <- b

c <- -50
for (i in 1:nrow(Valday_ag)){
  if (Valday_ag$T2M_MAX[i] > c){
    c <- Valday_ag$T2M_MAX[i]
  }
}
d30_df$T2M_MAX[9] <- c

d <- 50 
for (i in 1:nrow(Valday_ag)){
  if (Valday_ag$T2M_MIN[i] < d){
    d <- Valday_ag$T2M_MIN[i]
  }
}
d30_df$T2M_MIN[9] <- d
d30_df$T2M_RANGE[9] <- mean(Valday_ag$T2M_RANGE)
d30_df$FROST_DAYS[9] <- sum(Valday_ag$FROST_DAYS)
d30_df$PRECTOTCORR[9] <- mean(Valday_ag$PRECTOTCORR)
d30_df$ALLSKY_SFC_LW_DWN[9] <- mean(Valday_ag$ALLSKY_SFC_LW_DWN)

write.csv(d30_df,"all_30d.csv", row.names =F)

#PCA
setwd("C:/Users/Венера/Dropbox/UoL/other/ClimateComponents/NASA/all_param/")
d30_df <- read.csv("all_30d.csv")
d30_df <- d30_df[,6:19]
rownames(d30_df) <- c("Recarei", "Gimenells", "Karensminde", "Munich", "Mauternbach", "Akaa", "Uman", "Yesiloz", "Valday")

library("FactoMineR")
library("factoextra")
options(ggrepel.max.overlaps = Inf)
d30_pca <- PCA(d30_df, scale.unit = TRUE, graph = TRUE)
print(d30_pca)

# matrix with eigenvalues
d30_pca$eig
eig.val <- get_eigenvalue(d30_pca)
eig.val #first 2 PCs >1

#Scree plot
fviz_eig(d30_pca, addlabels = TRUE, xlab = "principal components")

fviz_pca_ind(d30_pca, axes = c(1, 2), col.ind = "coord", 
             title = "Principal Component Analysis",
             subtitle = "Bioclimatic variables",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)

#PCs (scores)
d30_pca$ind$coord

#Contribution to dimentions
# Change the gradient color
fviz_pca_var(d30_pca, col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# Contributions of variables to PC1
fviz_contrib(d30_pca, choice = "var", axes = 1, top = 19, title="Contribution of variables to PC1")

# Contributions of variables to PC2
fviz_contrib(d30_pca, choice = "var", axes = 2, top = 19, title="Contribution of variables to PC2")

#pca result dataframe
d30_df <- read.csv("all_30d.csv")
pops_data <- d30_df[,1:5]
PC1_2 <- d30_pca$ind$coord[1:9, 1:2] #there are 2 PCs with eigvalues >1
d30_pca_results <- cbind(pops_data, PC1_2)
write.csv(d30_pca_results,"all_d30_PCA.csv", row.names =T)

library(corrplot)
cor_df <- d30_pca_results[ -c(4:5) ]
corrplot(cor(cor_df), method = 'number') 
corrplot(cor(cor_df), type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)

#Phenotype analysis
setwd("C:/Users/Венера/Dropbox/UoL/other/ClimateComponents/NASA/all_param/phenotypes/")
d30_data <- read.csv("all_d30_PCA.csv", row.names = 1)
#this file is from DrosEU drive
pheno <-  read.csv("all_models_line_meta_compound_random_coefs_wide.csv", header=T)
str(pheno)
nrow(pheno)

new_cols <- c("Country", "PC1", "PC2")
pheno[ , new_cols] <- NA
for (i in 1:nrow(pheno)){
  if (pheno[i,"Population"] == "RE"){ #Portugal
    pheno[i, 35:37] <- d30_data[1,5:7]
  }
  if (pheno[i,"Population"] == "GI"){ #Spain
    pheno[i, 35:37] <- d30_data[2,5:7]
  }
  if (pheno[i,"Population"] == "KA"){ #Denmark
    pheno[i, 35:37] <- d30_data[3,5:7]
  }
  if (pheno[i,"Population"] == "MU"){ #Germany
    pheno[i, 35:37] <- d30_data[4,5:7]
  }
  if (pheno[i,"Population"] == "MA"){ #Austria
    pheno[i, 35:37] <- d30_data[5,5:7]
  }
  if (pheno[i,"Population"] == "AK"){ #Finland
    pheno[i, 35:37] <- d30_data[6,5:7]
  }
  if (pheno[i,"Population"] == "UM"){ #Ukraine
    pheno[i, 35:37] <- d30_data[7,5:7]
  }
  if (pheno[i,"Population"] == "YE"){ #Turkey
    pheno[i, 35:37] <- d30_data[8,5:7]
  }
  if (pheno[i,"Population"] == "VA"){ #Russia
    pheno[i, 35:37] <- d30_data[9,5:7]
  }
} 

#mixed-model with random effects
library(lme4)
library(lmerTest) 
pheno <- na.omit(pheno)

#CCRT_F
CCRT_F_model1 <- lmer(pheno$CCRT_F ~ pheno$PC1 + pheno$PC2+(1|(as.factor(pheno$Population))),data = pheno)
CCRT_F_model2 <- lmer(pheno$CCRT_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
CCRT_F_model3 <- glm(pheno$CCRT_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
CCRT_F_model4 <- glm(pheno$CCRT_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(CCRT_F_model1, CCRT_F_model2,CCRT_F_model3, CCRT_F_model4,test="Chisq")
summary(CCRT_F_model2)
Anova(CCRT_F_model2)

#CCRT_M
library(car)
CCRT_M_model1 <- lmer(pheno$CCRT_M ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
CCRT_M_model2 <- lmer(pheno$CCRT_M ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
CCRT_M_model3 <- glm(pheno$CCRT_M ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
CCRT_M_model4 <- glm(pheno$CCRT_M ~ pheno$PC1 + pheno$Population,data = pheno)
anova(CCRT_M_model1, CCRT_M_model2,CCRT_M_model3, CCRT_M_model4,test="Chisq")
summary(CCRT_M_model2)
Anova(CCRT_M_model3)

#CSM_F
CSM_F_model1 <- lmer(pheno$CSM_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
CSM_F_model2 <- lmer(pheno$CSM_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
CSM_F_model3 <- glm(pheno$CSM_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
CSM_F_model4 <- glm(pheno$CSM_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(CSM_F_model1, CSM_F_model2,CSM_F_model3, CSM_F_model4,test="Chisq")
summary(CSM_F_model1)
Anova(CSM_F_model3)

#CSM_M
CSM_M_model1 <- lmer(pheno$CSM_M ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
CSM_M_model2 <- lmer(pheno$CSM_M ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
CSM_M_model3 <- glm(pheno$CSM_M ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
CSM_M_model4 <- glm(pheno$CSM_M ~ pheno$PC1 + pheno$Population,data = pheno)
anova(CSM_M_model1, CSM_M_model2,CSM_M_model3, CSM_M_model4,test="Chisq")
summary(CSM_M_model1)
Anova(CSM_M_model3)

#DT_A_F
DT_A_F_model1 <- lmer(pheno$DT_A_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
DT_A_F_model2 <- lmer(pheno$DT_A_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
DT_A_F_model3 <- glm(pheno$DT_A_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
DT_A_F_model4 <- glm(pheno$DT_A_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(DT_A_F_model1, DT_A_F_model2,DT_A_F_model3, DT_A_F_model4,test="Chisq")
summary(DT_A_F_model2)
Anova(DT_A_F_model3)

#DT_A_M
DT_A_M_model1 <- lmer(pheno$DT_A_M ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
DT_A_M_model2 <- lmer(pheno$DT_A_M ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
DT_A_M_model3 <- glm(pheno$DT_A_M ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
DT_A_M_model4 <- glm(pheno$DT_A_M ~ pheno$PC1 + pheno$Population,data = pheno)
anova(DT_A_M_model1, DT_A_M_model2,DT_A_M_model3, DT_A_M_model4,test="Chisq")
summary(DT_A_M_model2)
Anova(DT_A_M_model3)

#DT_P_NA
DT_P_NA_model1 <- lmer(pheno$DT_P_NA ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
DT_P_NA_model2 <- lmer(pheno$DT_P_NA ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
DT_P_NA_model3 <- glm(pheno$DT_P_NA ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
DT_P_NA_model4 <- glm(pheno$DT_P_NA ~ pheno$PC1 + pheno$Population,data = pheno)
anova(DT_P_NA_model1, DT_P_NA_model2,DT_P_NA_model3, DT_P_NA_model4,test="Chisq")
summary(DT_P_NA_model2)
Anova(DT_P_NA_model3)

#Dia_F
Dia_F_model1 <- lmer(pheno$Dia_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
Dia_F_model2 <- lmer(pheno$Dia_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
Dia_F_model3 <- glm(pheno$Dia_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
Dia_F_model4 <- glm(pheno$Dia_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(Dia_F_model1, Dia_F_model2,Dia_F_model3, Dia_F_model4,test="Chisq")
summary(Dia_F_model2)
Anova(Dia_F_model3)

#DW_F
DW_F_model1 <- lmer(pheno$DW_F ~ pheno$PC1 + pheno$PC2 +(1|(as.factor(pheno$Population))),data = pheno)
DW_F_model2 <- lmer(pheno$DW_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
DW_F_model3 <- glm(pheno$DW_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
DW_F_model4 <- glm(pheno$DW_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(DW_F_model1, DW_F_model2,DW_F_model3, DW_F_model4,test="Chisq")
summary(DW_F_model2)
Anova(DW_F_model3)

#DW_M
DW_M_model1 <- lmer(pheno$DW_M ~ pheno$PC1 + pheno$PC2+(1|(as.factor(pheno$Population))),data = pheno)
DW_M_model2 <- lmer(pheno$DW_M ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
DW_M_model3 <- glm(pheno$DW_M ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
DW_M_model4 <- glm(pheno$DW_M ~ pheno$PC1 + pheno$Population,data = pheno)
anova(DW_M_model1, DW_M_model2,DW_M_model3, DW_M_model4,test="Chisq")
summary(DW_M_model2)
Anova(DW_M_model3)

#Fec_F
Fec_F_model1 <- lmer(pheno$Fec_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
Fec_F_model2 <- lmer(pheno$Fec_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
Fec_F_model3 <- glm(pheno$Fec_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
Fec_F_model4 <- glm(pheno$Fec_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(Fec_F_model1, Fec_F_model2,Fec_F_model3, Fec_F_model4,test="Chisq")
summary(Fec_F_model1)
Anova(Fec_F_model3)

#HSM_F
HSM_F_model1 <- lmer(pheno$HSM_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
HSM_F_model2 <- lmer(pheno$HSM_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
HSM_F_model3 <- glm(pheno$HSM_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
HSM_F_model4 <- glm(pheno$HSM_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(HSM_F_model1, HSM_F_model2,HSM_F_model3, HSM_F_model4,test="Chisq")
summary(HSM_F_model2)
Anova(HSM_F_model3)

#HSM_M
HSM_M_model1 <- lmer(pheno$HSM_M ~ pheno$PC1 + pheno$PC2+(1|(as.factor(pheno$Population))),data = pheno)
HSM_M_model2 <- lmer(pheno$HSM_M ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
HSM_M_model3 <- glm(pheno$HSM_M ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
HSM_M_model4 <- glm(pheno$HSM_M ~ pheno$PC1 + pheno$Population,data = pheno)
anova(HSM_M_model1, HSM_M_model2,HSM_M_model3, HSM_M_model4,test="Chisq")
summary(HSM_M_model1)
Anova(HSM_M_model3)

#LS_F
LS_F_model1 <- lmer(pheno$LS_F ~ pheno$PC1 + pheno$PC2 +(1|(as.factor(pheno$Population))),data = pheno)
LS_F_model2 <- lmer(pheno$LS_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
LS_F_model3 <- glm(pheno$LS_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
LS_F_model4 <- glm(pheno$LS_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(LS_F_model1, LS_F_model2,LS_F_model3, LS_F_model4,test="Chisq")
summary(LS_F_model1)
Anova(LS_F_model3)

#LS_M
LS_M_model1 <- lmer(pheno$LS_M ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
LS_M_model2 <- lmer(pheno$LS_M ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
LS_M_model3 <- glm(pheno$LS_M ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
LS_M_model4 <- glm(pheno$LS_M ~ pheno$PC1 + pheno$Population,data = pheno)
anova(LS_M_model1, LS_M_model2,LS_M_model3, LS_M_model4,test="Chisq")
summary(LS_M_model1)
Anova(LS_M_model3)


#LA_Activity_B
LA_Activity_B_model1 <- lmer(pheno$LA_Activity_B ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
LA_Activity_B_model2 <- lmer(pheno$LA_Activity_B ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
LA_Activity_B_model3 <- glm(pheno$LA_Activity_B ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
LA_Activity_B_model4 <- glm(pheno$LA_Activity_B ~ pheno$PC1 + pheno$Population,data = pheno)
anova(LA_Activity_B_model1, LA_Activity_B_model2,LA_Activity_B_model3, LA_Activity_B_model4,test="Chisq")
summary(LA_Activity_B_model1)
Anova(LA_Activity_B_model3)

#LA_CircPhase_B
LA_CircPhase_B_model1 <- lmer(pheno$LA_CircPhase_B ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
LA_CircPhase_B_model2 <- lmer(pheno$LA_CircPhase_B ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
LA_CircPhase_B_model3 <- glm(pheno$LA_CircPhase_B ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
LA_CircPhase_B_model4 <- glm(pheno$LA_CircPhase_B ~ pheno$PC1 + pheno$Population,data = pheno)
anova(LA_CircPhase_B_model1, LA_CircPhase_B_model2,LA_CircPhase_B_model3, LA_CircPhase_B_model4,test="Chisq")
summary(LA_CircPhase_B_model2)
Anova(LA_CircPhase_B_model3)

#LA_NDlog2_B
LA_NDlog2_B_model1 <- lmer(pheno$LA_NDlog2_B ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
LA_NDlog2_B_model2 <- lmer(pheno$LA_NDlog2_B ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
LA_NDlog2_B_model3 <- glm(pheno$LA_NDlog2_B ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
LA_NDlog2_B_model4 <- glm(pheno$LA_NDlog2_B ~ pheno$PC1 + pheno$Population,data = pheno)
anova(LA_NDlog2_B_model1, LA_NDlog2_B_model2,LA_NDlog2_B_model3, LA_NDlog2_B_model4,test="Chisq")
summary(LA_NDlog2_B_model1)
Anova(LA_NDlog2_B_model3)

#LA_Period_B
LA_Period_B_model1 <- lmer(pheno$LA_Period_B ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
LA_Period_B_model2 <- lmer(pheno$LA_Period_B ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
LA_Period_B_model3 <- glm(pheno$LA_Period_B ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
LA_Period_B_model4 <- glm(pheno$LA_Period_B ~ pheno$PC1 + pheno$Population,data = pheno)
anova(LA_Period_B_model1, LA_Period_B_model2,LA_Period_B_model3, LA_Period_B_model4,test="Chisq")
summary(LA_Period_B_model2)
Anova(LA_Period_B_model3)

#Pgm_T4_F
Pgm_T4_F_model1 <- lmer(pheno$Pgm_T4_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
Pgm_T4_F_model2 <- lmer(pheno$Pgm_T4_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
Pgm_T4_F_model3 <- glm(pheno$Pgm_T4_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
Pgm_T4_F_model4 <- glm(pheno$Pgm_T4_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(Pgm_T4_F_model1, Pgm_T4_F_model2,Pgm_T4_F_model3, Pgm_T4_F_model4,test="Chisq")
summary(Pgm_T4_F_model1)
Anova(Pgm_T4_F_model3)

#Pgm_T5_F
Pgm_T5_F_model1 <- lmer(pheno$Pgm_T5_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
Pgm_T5_F_model2 <- lmer(pheno$Pgm_T5_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
Pgm_T5_F_model3 <- glm(pheno$Pgm_T5_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
Pgm_T5_F_model4 <- glm(pheno$Pgm_T5_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(Pgm_T5_F_model1, Pgm_T5_F_model2,Pgm_T5_F_model3, Pgm_T5_F_model4,test="Chisq")
summary(Pgm_T5_F_model2)
Anova(Pgm_T5_F_model3)

#Pgm_T6_F
Pgm_T6_F_model1 <- lmer(pheno$Pgm_T6_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
Pgm_T6_F_model2 <- lmer(pheno$Pgm_T6_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
Pgm_T6_F_model3 <- glm(pheno$Pgm_T6_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
Pgm_T6_F_model4 <- glm(pheno$Pgm_T6_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(Pgm_T6_F_model1, Pgm_T6_F_model2,Pgm_T6_F_model3, Pgm_T6_F_model4,test="Chisq")
summary(Pgm_T6_F_model2)
Anova(Pgm_T6_F_model3)

#Pgm_Total_F
Pgm_Total_F_model1 <- lmer(pheno$Pgm_Total_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
Pgm_Total_F_model2 <- lmer(pheno$Pgm_Total_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
Pgm_Total_F_model3 <- glm(pheno$Pgm_Total_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
Pgm_Total_F_model4 <- glm(pheno$Pgm_Total_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(Pgm_Total_F_model1, Pgm_Total_F_model2,Pgm_Total_F_model3, Pgm_Total_F_model4,test="Chisq")
summary(Pgm_Total_F_model2)
Anova(Pgm_Total_F_model3)

#SR_F
SR_F_model1 <- lmer(pheno$SR_F ~ pheno$PC1 + pheno$PC2 +(1|(as.factor(pheno$Population))),data = pheno)
SR_F_model2 <- lmer(pheno$SR_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
SR_F_model3 <- glm(pheno$SR_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
SR_F_model4 <- glm(pheno$SR_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(SR_F_model1, SR_F_model2,SR_F_model3, SR_F_model4,test="Chisq")
summary(SR_F_model1)
Anova(SR_F_model3)

#SR_M
SR_M_model1 <- lmer(pheno$SR_M ~ pheno$PC1 + pheno$PC2 +(1|(as.factor(pheno$Population))),data = pheno)
SR_M_model2 <- lmer(pheno$SR_M ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
SR_M_model3 <- glm(pheno$SR_M ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
SR_M_model4 <- glm(pheno$SR_M ~ pheno$PC1 + pheno$Population,data = pheno)
anova(SR_M_model1, SR_M_model2,SR_M_model3, SR_M_model4,test="Chisq")
summary(SR_M_model1)
Anova(SR_M_model3)

#TL_F
TL_F_model1 <- lmer(pheno$TL_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
TL_F_model2 <- lmer(pheno$TL_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
TL_F_model3 <- glm(pheno$TL_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
TL_F_model4 <- glm(pheno$TL_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(TL_F_model1, TL_F_model2,TL_F_model3, TL_F_model4,test="Chisq")
summary(TL_F_model1)
Anova(TL_F_model3)


#TL_M
TL_M_model1 <- lmer(pheno$TL_M ~ pheno$PC1 + pheno$PC2+(1|(as.factor(pheno$Population))),data = pheno)
TL_M_model2 <- lmer(pheno$TL_M ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
TL_M_model3 <- glm(pheno$TL_M ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
TL_M_model4 <- glm(pheno$TL_M ~ pheno$PC1 + pheno$Population,data = pheno)
anova(TL_M_model1, TL_M_model2,TL_M_model3, TL_M_model4,test="Chisq")
summary(TL_M_model1)
Anova(TL_M_model3)


#Via_NA
Via_NA_model1 <- lmer(pheno$Via_NA ~ pheno$PC1 + pheno$PC2 +(1|(as.factor(pheno$Population))),data = pheno)
Via_NA_model2 <- lmer(pheno$Via_NA ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
Via_NA_model3 <- glm(pheno$Via_NA ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
Via_NA_model4 <- glm(pheno$Via_NA ~ pheno$PC1 + pheno$Population,data = pheno)
anova(Via_NA_model1, Via_NA_model2,Via_NA_model3, Via_NA_model4,test="Chisq")
summary(Via_NA_model1)
Anova(Via_NA_model3)

#WA_L_F
WA_L_F_model1 <- lmer(pheno$WA_L_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
WA_L_F_model2 <- lmer(pheno$WA_L_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
WA_L_F_model3 <- glm(pheno$WA_L_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
WA_L_F_model4 <- glm(pheno$WA_L_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(WA_L_F_model1, WA_L_F_model2,WA_L_F_model3, WA_L_F_model4,test="Chisq")
summary(WA_L_F_model1)
Anova(WA_L_F_model3)


#WA_L_M
WA_L_M_model1 <- lmer(pheno$WA_L_M ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
WA_L_M_model2 <- lmer(pheno$WA_L_M ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
WA_L_M_model3 <- glm(pheno$WA_L_M ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
WA_L_M_model4 <- glm(pheno$WA_L_M ~ pheno$PC1 + pheno$Population,data = pheno)
anova(WA_L_M_model1, WA_L_M_model2,WA_L_M_model3, WA_L_M_model4,test="Chisq")
summary(WA_L_M_model2)
Anova(WA_L_M_model3)

#WA_R_F
WA_R_F_model1 <- lmer(pheno$WA_R_F ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
WA_R_F_model2 <- lmer(pheno$WA_R_F ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
WA_R_F_model3 <- glm(pheno$WA_R_F ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
WA_R_F_model4 <- glm(pheno$WA_R_F ~ pheno$PC1 + pheno$Population,data = pheno)
anova(WA_R_F_model1, WA_R_F_model2,WA_R_F_model3, WA_R_F_model4,test="Chisq")
summary(WA_R_F_model1)
Anova(WA_R_F_model3)

#WA_R_M
WA_R_M_model1 <- lmer(pheno$WA_R_M ~ pheno$PC1 + pheno$PC2  +(1|(as.factor(pheno$Population))),data = pheno)
WA_R_M_model2 <- lmer(pheno$WA_R_M ~ pheno$PC1 + (1|(as.factor(pheno$Population))),data = pheno)
WA_R_M_model3 <- glm(pheno$WA_R_M ~ pheno$PC1 + pheno$PC2 + pheno$Population,data = pheno)
WA_R_M_model4 <- glm(pheno$WA_R_M ~ pheno$PC1 + pheno$Population,data = pheno)
anova(WA_R_M_model1, WA_R_M_model2,WA_R_M_model3, WA_R_M_model4,test="Chisq")
summary(WA_R_M_model2)
Anova(WA_R_M_model3)
