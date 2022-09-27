### libraries
  library(data.table)
  library(gdata)
  library(nasapower)
  library(lubridate)
  library(foreach)
  library(doMC)
  registerDoMC(4)
  library(readxl)

### load sample data
  samps <- read_excel("/Users/alanbergland/Documents/GitHub/misc/DrosEU_phenotype_2022/DrosEU_PhenotypingWG_SamplingDates.xlsx")
  samps <- as.data.table(samps)
  samps[,month:=as.numeric(tstrsplit(midDate, "/")[[2]])]
  samps[,day:=as.numeric(tstrsplit(midDate, "/")[[1]])]
  samps[,year:=as.numeric(tstrsplit(midDate, "/")[[3]])]
  samps[,Date:=date(paste(year, month, day, sep="-"))]

### function to use nasapower
  getPower <- function(i) {
      daily_single_ag <- get_power(
        community = "ag",
        lonlat = c(samps[i]$Longitude, samps[i]$Latitude),
        pars = c("RH2M", "T2M", "PRECTOTCORR"),
        dates = c(paste(samps[i]$year, "-01-01", sep=""), paste(samps[i]$year, "-12-31", sep="")),
        temporal_api = "hourly",
        time_standard="UTC"
      )
      daily_single_ag <- as.data.table(daily_single_ag)
      daily_single_ag
  }

### iterate through
  power.dt <- list()

  for(i in 1:dim(samps)[1]) {
    power.dt[[i]] <- getPower(i)
    Sys.sleep(5)
    power.dt[[i]][,sampleId:=samps[i]$LocationAbbr]
    power.dt[[i]][,date:=ymd_hms(paste(paste(YEAR, MO, DY, sep="-"), paste(HR, ":00:00", sep=""), sep=" "))]
    power.dt[[i]][,queryDate:=Sys.time()]

  }
  power.dt <- rbindlist(power.dt)

### merge
  power.dt <- merge(power.dt, samps[,c("LocationAbbr", "Date"),with=F], by.x="sampleId", by.y="LocationAbbr")
  setnames(power.dt, "Date", "collectionDate")
  setnames(power.dt, c("T2M", "RH2M", "PRECTOTCORR"), c("temp", "humidity", "precip"))

### summarize
  sets <- data.table(mod=c(1:7),
                     start=c(0, 0,  0,  0,  0,   0, 0),
                     end=	 c(1, 7, 30, 60, 90, 180, 360))

  setkey(power.dt, sampleId)
  weather.ave <- foreach(i=unique(power.dt$sampleId), .combine="rbind")%dopar%{
    # i <- unique(power.dt$sampleId)[1]
    power.tmp <- power.dt[J(i)]

    power.tmp[,delta:=julian(date)-julian(collectionDate)]

    foreach(k=1:dim(sets)[1], .combine="rbind")%do%{
      message(paste(i, k, sep=" / "))
      power.mod <- power.tmp[delta>= -1*(sets[k]$end) & delta<=-1*(sets[k]$start)]

      power.mod.ag <- power.mod[,list(dailyMax=max((temp), na.rm=T),
                                  dailyMin=min((temp), na.rm=T)),
                             list(delta=round(delta))]

      data.table(sampleId=i, mod=k,
                temp.ave=mean(power.mod$temp),
                temp.var=var(power.mod$temp),
                temp.max=max(power.mod$temp),
                temp.min=min(power.mod$temp),
                temp.propMax=mean(power.mod.ag$dailyMax>32),
                temp.propmin=mean(power.mod.ag$dailyMin<5),

                humidity.ave=mean(power.mod$humidity),
                humidity.var=var(power.mod$humidity),

                precip.ave=mean(power.mod$precip),
                precip.var=var(power.mod$precip))
      }

    }
  wl <- melt(weather.ave, id.vars=c("sampleId", "mod"))
  summary(weather.ave)

### save
  save(power.dt, weather.ave, file="/Users/alanbergland/Documents/GitHub/misc/DrosEU_phenotype_2022/weather_ave.Rdata")
  load(file="/Users/alanbergland/Documents/GitHub/misc/DrosEU_phenotype_2022/weather_ave.Rdata")

### load in phenotype data
  #pheno <- fread("/Users/alanbergland/Documents/GitHub/DrosEU_PhenotypingWG/LinearModelsPop/all_models_line_compound_random_coefs.csv")
  #pheno <- fread("/Users/alanbergland/Documents/GitHub/DrosEU_PhenotypingWG/MetaAnalyses/all_models_pop_meta_compound_estimates.csv")


  pheno <- fread("/Users/alanbergland/Documents/GitHub/DrosEU_PhenotypingWG/MetaAnalyses/all_models_line_meta_compound_random_coefs_wide.csv")
  pheno <- melt(pheno, id.vars=c("Population", "Line"), variable.name="Trait")
  setnames(pheno, "value", "Estimate")
  pheno[,Sex:=tstrsplit(Trait, "_")[[2]]]
  pheno[!Sex%in%c("F", "M"), Sex:=NA]
  pheno[,Trait:=gsub("_F", "", Trait)]
  pheno[,Trait:=gsub("_M", "", Trait)]

  m <- merge(wl, pheno, by.x="sampleId", by.y="Population", allow.cartesian=T)

  o <-
  foreach(mod.i=unique(m$mod), .combine="rbind", .errorhandling="remove")%do%{
    foreach(var.i=unique(m$variable), .combine="rbind", .errorhandling="remove")%do%{
      foreach(trait.i=unique(m$Trait), .combine="rbind", .errorhandling="remove")%do%{
        foreach(sex.i=c("F", "M"), .combine="rbind", .errorhandling="remove")%dopar%{
          # mod.i=1; var.i="temp.min"; trait.i="CSM"; sex.i <- "M"
          message(paste(mod.i, var.i, trait.i, sep=" / "))
          tmp <- m[mod==mod.i][variable==var.i][Trait==trait.i][Sex==sex.i]
          tmp[,pheno_norm:=(Estimate-mean(Estimate, na.rm=T))/sd(Estimate, na.rm=T)]
          tmp[,env_norm:=(value-mean(value))/sd(value)]

          tmp.ag <- tmp[,list(pheno_norm_mu=mean(pheno_norm), env_norm=mean(env_norm)), list(sampleId)]

          t1 <- lm(pheno_norm~env_norm, data=tmp)
          t0 <- lm(pheno_norm~1, data=tmp)
          t1a <- lm(pheno_norm_mu~env_norm, data=tmp.ag)
          t0a <- lm(pheno_norm_mu~1, data=tmp.ag)

        #  t2 <- lm(Value~value, data=tmp)
          #t3 <- lmer(pheno_norm~env_norm + (1|sampleId), data=tmp)
          #t3a <- lmer(pheno_norm~1 + (1|sampleId), data=tmp)

          es <- eta_squared(t1 , partial=F)
          esa <- eta_squared(t1a , partial=F)

          rbind(
          data.table(mod=mod.i, var=var.i, trait=trait.i, sex=sex.i,
                      Eta2=es$Eta2, eta2_lci=es$CI_low, eta2_uci=es$CI_high,
                      Eta2a=esa$Eta2, eta2_lci_a=esa$CI_low, eta2_uci_a=esa$CI_high,
                      beta_norm=summary(t1)$coef[2,1], p_norm=summary(t1)$coef[2,4],
                      beta_norm_a=summary(t1a)$coef[2,1], p_norm_a=summary(t1a)$coef[2,4],
                      AIC=AIC(t1),
                      AICa=AIC(t1a)),
          data.table(mod=NA, var="NULL", trait=trait.i, sex=sex.i,
                      AIC=AIC(t0), AICa=AIC(t0a)),
          fill=T)

                    #  beta_raw=summary(t2)$coef[2,1], p_raw=summary(t2)$coef[2,4],
                    #  beta_lmm=summary(t3)$coef[2,1], lrt_lmm=anova(t3, t3a, test="Chisq")[2,8])
        }
      }
    }
  }

  o.aic <- o[,list(mod=mod[which.min(AIC)], var=var[which.min(AIC)], minAIC=min(AIC)), list(trait, sex)]
  setkey(o,     trait, sex)
  setkey(o.aic, trait, sex)

  o <- merge(o, o.aic[,-"mod"][,-"var"])
  o[,deltaAIC:=AIC - minAIC]

### save
  save(m, o, file="/Users/alanbergland/Documents/GitHub/misc/DrosEU_phenotype_2022/nasa_power_correlation_corrected.Rdata")

### here is another version
  o[, pa:=p.adjust(p_norm_a)]
  sets[,set:=paste("Days ", start, "-", end, "\nprior to collection", sep="")]

  o <- merge(o, sets, by.x="mod", by.y="mod")

  o[,set:=factor(mod, levels=sets$set)]

  tilep <-
  ggplot(data=o, aes(x=var, y=trait, fill=Eta2a)) +
  geom_tile() +
  geom_tile(data=o[deltaAIC<2], aes(fill=Eta2a), color="green", size=1) +
  geom_tile(data=o[deltaAIC==0], aes(fill=Eta2a), color="orange", size=1) +
  geom_point(data=o[eta2_lci_a>0], aes(color=as.factor(beta_norm_a>0))) +
  facet_grid(set.y~sex, scales="free_x") + coord_flip() +
  theme_bw() +
  scale_fill_viridis(option="E") +
  theme(panel.border = element_blank(),
       panel.background = element_blank(),
       panel.grid = element_blank(),
       panel.spacing = unit(.05,"line")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.background = element_rect(size = 0.5))

  ggsave(tilep, file="~/tileplot4.png",h=11, w=11)



o[trait=="CSM"][sex=="F"][mod==3]












  ### plot
    p1 <-
    ggplot(data=o) +
    geom_point(aes(x=mod, y=-log10(p))) +
    facet_grid(trait~var) +
    geom_hline(aes(yintercept=-log10(.05)))

    ggsave(p1, file="~/drosEU_nasapowe.pdf", h=14, w=12)


    o.ag <- o[,list(nSig=sum(pa<.05), uniqP=length(unique(trait[pa<.05]))), list(mod, var)]

    ggplot(data=o.ag, aes(x=mod, y=nSig)) + geom_point() + facet_grid(~var)

    m[,pheno:=Value]
    m[,env:=value]



#### here is one version
  modmod <- 3
  o[,pa:=NA]
  o[, pa:=p.adjust(p_norm)]


  examp <-
  ggplot(data=m[variable=="temp.max"][mod==modmod][Trait=="HSM"],
        aes(y=pheno, x=env, color=sampleId)) +
  geom_point() +
  facet_grid(~Sex) +
  geom_smooth(method = "lm", se=F) +
  geom_text(data=m[variable=="temp.max"][mod==modmod][Trait=="HSM"][,list(pheno=max(pheno), env=mean(env)), list(sampleId, Sex)],
            aes(x=env, y=pheno*1.05, label=sampleId)) +
  xlab("Temp.max 0-30") + ylab("HSM") +
  ggtitle("Heat Shock Mortality ~ Temp.max 0-30")

  o[var=="temp.max"][mod==modmod][trait=="HSM"]



  male <-
  ggplot(data=o[mod==modmod][sex=="M"]) +
  geom_line(aes(x=beta_norm, y=-log10(p_norm), group=interaction(sex, mod), linetype=sex)) +
  geom_point(aes(x=beta_norm, y=-log10(p_norm), color=as.factor(pa<.05), shape=sex)) +
  facet_wrap(~trait, ncol=5) +
  geom_hline(aes(yintercept=-log10(.05))) +
  ggtitle("Males, days 0-30 prior to collection") +
  geom_text_repel(data=o[mod==modmod][pa<.05][sex=="M"], aes(label=var, x=beta_norm, y=-log10(p_norm)), box.padding = 0.5, max.overlaps = Inf)

  female <-
  ggplot(data=o[mod==modmod][sex=="F"]) +
  geom_line(aes(x=beta_norm, y=-log10(p_norm), group=interaction(sex, mod), linetype=sex)) +
  geom_point(aes(x=beta_norm, y=-log10(p_norm), color=as.factor(pa<.05), shape=sex)) +
  facet_wrap(~trait, ncol=5) +
  geom_hline(aes(yintercept=-log10(.05))) +
  ggtitle("Females, days 0-30 prior to collection") +
  geom_text_repel(data=o[mod==modmod][pa<.05][sex=="F"], aes(label=var, x=beta_norm, y=-log10(p_norm)), box.padding = 0.5, max.overlaps = Inf)

  layout <- "
  AABB
  AABB
  CCCC"

  mega <- male + female + examp + plot_layout(design=layout)
  ggsave(mega, file="~/nasa_power_pheno.pdf", h=15, w=20)

### noddle
  o[,pa:=p.adjust(p, "fdr")]

  o.ag <- o[mod==1 & var=="temp.ave",list(nSig=mean(p<.05), n=length(p)), list(trait, var)]
  o.ag


  o[trait=="LS"][var=="humidity.var"][mod==3]




  ### load in phenotype data
    pheno <- fread("/Users/alanbergland/Documents/GitHub/DrosEU_PhenotypingWG/MetaAnalyses/all_metas_pop_compound_estimates.csv")

    m <- merge(wl, pheno, by.x="sampleId", by.y="Population", allow.cartesian=T)

    o2 <-
    foreach(mod.i=unique(m$mod), .combine="rbind", .errorhandling="remove")%do%{
      foreach(var.i=unique(m$variable), .combine="rbind", .errorhandling="remove")%do%{
        foreach(trait.i=unique(m$Trait), .combine="rbind", .errorhandling="remove")%dopar%{
          foreach(sex.i=unique(m$Sex), .combine="rbind", .errorhandling="remove")%do%{
            # mod.i=3; var.i="temp.ave"; trait.i="CCRT"; sex.i="M"
            message(paste(mod.i, var.i, trait.i, sex.i, sep=" / "))
            t1 <- lm(Mstar~value, m[mod==mod.i][variable==var.i][Trait==trait.i][Sex==sex.i])
            data.table(mod=mod.i, var=var.i, trait=trait.i, sex=sex.i,
                        beta=summary(t1)$coef[2,1], p=summary(t1)$coef[2,4])
          }
        }
      }
    }
    o2[,pa:=p.adjust(p)]
    table(o2$var, o2$p<.05)
55/264
  o2[var=="temp.ave"]
