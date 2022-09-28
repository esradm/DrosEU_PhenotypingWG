#setwd("/Users/venera/Dropbox/UoL/other/ClimateComponents/NASA/results/traits_combined/")
setwd("C:/Users/Венера/Dropbox/UoL/other/ClimateComponents/NASA/results/traits_combined/")

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
setwd("/Users/Венера/Dropbox/UoL/other/ClimateComponents/NASA/results/traits_combined/")
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
colnames(d30_pca_results) <- c("Longitude", "Latitude", "Altitude", "Population",
                               "Country", "PC1_clim", "PC2_clim")
write.csv(d30_pca_results,"all_30d_PCA.csv", row.names =F)

library(corrplot)
cor_df <- d30_pca_results[ -c(4:5) ]
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
climate_data <- read.csv("all_30d_PCA.csv")
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
write.csv(F9_data,"F9_30d_data.csv", row.names =F)

#phenotypes
#PC1
F9_PC1_1 <- lmer(F9_data$PC1_F9 ~ F9_data$PC1_clim + F9_data$PC2_clim +(1|(as.factor(F9_data$Population))),data = F9_data)
summary(F9_PC1_1)
F9_PC1_2 <- lmer(F9_data$PC1_F9 ~ F9_data$PC1_clim +(1|(as.factor(F9_data$Population))),data = F9_data)
summary(F9_PC1_2)

F9_PC1_3 <- glm(F9_data$PC1_F9 ~ F9_data$PC1_clim + F9_data$PC2_clim,data = F9_data)
summary(F9_PC1_3)
F9_PC1_4 <- glm(F9_data$PC1_F9 ~ F9_data$PC2_clim,data = F9_data)
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
climate_data <- read.csv("all_30d_PCA.csv")
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

write.csv(M9_data,"M9_30d_data.csv", row.names =F)

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
M9_PC2_2 <- lmer(M9_data$PC2_M9 ~ M9_data$PC2_clim +(1|(as.factor(M9_data$Population))),data = M9_data)
summary(M9_PC2_2)

M9_PC2_3 <- glm(M9_data$PC2_M9 ~ M9_data$PC1_clim + M9_data$PC2_clim,data = M9_data)
summary(M9_PC2_3)

anova(M9_PC2_1,M9_PC2_2,M9_PC2_3, test="Chisq")

#females max
setwd("/Users/Венера/Dropbox/UoL/other/ClimateComponents/NASA/results/traits_combined/")
F9max_all_PC <- read.csv("/Users/Венера/Dropbox/UoL/other/ClimateComponents/NASA/results/Ewan/Fmax_drosEU_PCcoords.csv")
F9max_data <- F9max_all_PC[,1:5]
colnames(F9max_data) <- c("Country", "Population", "Line", "PC1_F9max", "PC2_F9max")
climate_data <- read.csv("all_30d_PCA.csv")
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
write.csv(F9max_data,"F9max_30d_data.csv", row.names =F)


#phenotypes
#PC1
F9max_PC1_1 <- lmer(F9max_data$PC1_F9max ~ F9max_data$PC1_clim + F9max_data$PC2_clim +(1|(as.factor(F9max_data$Population))),data = F9max_data)
summary(F9max_PC1_1)
F9max_PC1_2 <- lmer(F9max_data$PC1_F9max ~ F9max_data$PC2_clim +(1|(as.factor(F9max_data$Population))),data = F9max_data)
summary(F9max_PC1_2)

F9max_PC1_3 <- glm(F9max_data$PC1_F9max ~ F9max_data$PC1_clim + F9max_data$PC2_clim,data = F9max_data)
summary(F9max_PC1_3)
F9max_PC1_4 <- glm(F9max_data$PC1_F9max ~ F9max_data$PC2_clim,data = F9max_data)
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
climate_data <- read.csv("all_30d_PCA.csv")
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
write.csv(F9maxP_data,"F9maxP_30d_data.csv", row.names =F)


#phenotypes
#PC1
F9maxP_PC1_1 <- lmer(F9maxP_data$PC1_F9maxP ~F9maxP_data$PC1_clim+ F9maxP_data$PC1_clim + F9maxP_data$PC2_clim +(1|(as.factor(F9maxP_data$Population))),data = F9maxP_data)
summary(F9maxP_PC1_1)
F9maxP_PC1_2 <- lmer(F9maxP_data$PC1_F9maxP ~ F9maxP_data$PC2_clim +(1|(as.factor(F9maxP_data$Population))),data = F9maxP_data)
summary(F9maxP_PC1_2)

F9maxP_PC1_3 <- glm(F9maxP_data$PC1_F9maxP ~ F9maxP_data$PC1_clim + F9maxP_data$PC2_clim,data = F9maxP_data)
summary(F9maxP_PC1_3)
F9maxP_PC1_4 <- glm(F9maxP_data$PC1_F9maxP ~ F9maxP_data$PC2_clim,data = F9maxP_data)
summary(F9maxP_PC1_4)

anova(F9maxP_PC1_1,F9maxP_PC1_2,F9maxP_PC1_3,F9maxP_PC1_4, test="Chisq")

#PC2
F9maxP_PC2_1 <- lmer(F9maxP_data$PC2_F9maxP ~ F9maxP_data$PC1_clim + F9maxP_data$PC2_clim +(1|(as.factor(F9maxP_data$Population))),data = F9maxP_data)
summary(F9maxP_PC2_1)

F9maxP_PC2_2 <- glm(F9maxP_data$PC2_F9maxP ~ F9maxP_data$PC1_clim + F9maxP_data$PC2_clim,data = F9maxP_data)
summary(F9maxP_PC2_2)
