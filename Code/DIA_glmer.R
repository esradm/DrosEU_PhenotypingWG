### library
  library(data.table)
  library(lme4)
  library(dplyr)
  library(foreach)

### data
  #dia <- fread("/Users/alanbergland/Documents/GitHub/DrosEU_PhenotypingWG/Data/MasterSheets_May22_git/DIA_MasterSheet_Feb22.csv")
  dia <- fread("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/Data/MasterSheets_May22_git/DIA_MasterSheet_Feb22.csv")

###
  source("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/Functions/lab_correlations_functions.R")

### my version
  dia <- dia[,list(Prop_Max_Stage9=mean(MostAdvancedStage<=9 & NumberOfEggs==0), .N),
              list(Supervisor.PI, Diet, Batch, Population, Line)]

### model estiamtes per population, linear model
  PI <- unique(dia$Supervisor.PI)
  dia.est <- foreach(pii = PI, .combine="rbind", .errorhandling="remove")%do%{


    dia.glmer <- glmer(Prop_Max_Stage9 ~ Population + (1|Line:Population),
                        data=dia[Supervisor.PI==pii],
                        weights=dia[Supervisor.PI==pii]$N,
                        family=binomial())

    line <- getEstSE3(dia.glmer)
  #  line[,fitted:=fitted(dia.glmer)]
    pop <- as.data.table(getEstSE(dia.glmer))
    line[,pii:=pii]
    pop[,pii:=pii]

    list(line, pop)
  }

  line <- rbindlist(dia.est[c(1,2,3)])
  pop <- rbindlist(dia.est[c(4,5,6)])

  pop.w <- dcast(pop, Population~ pii, value.var="Estimate")

  cor.test(pop.w$Bergland, pop.w$Flatt)
  cor.test(pop.w$Bergland, pop.w$Schlotterer)
  cor.test(pop.w$Flatt, pop.w$Schlotterer)

  cor.test(rank(pop.w$Bergland), rank(pop.w$Schlotterer))


  line.w <- dcast(line, gr~ pii, value.var="Estimate")

  cor.test(line.w$Bergland, line.w$Flatt)
  cor.test(line.w$Bergland, line.w$Schlotterer)
  cor.test(line.w$Flatt, line.w$Schlotterer)

  cor.test(rank(line.w$Bergland), rank(line.w$Schlotterer))

  library(GGally)
  ggpairs(line.w[,-1])

  line[,pop:=tstrsplit(gr, ":")[[1]]]
  line[,line:=tstrsplit(gr, ":")[[2]]]

  ggplot(data=line, aes(x=pop, y=Estimate, color=pii, group=line)) +
  geom_point(position=position_dodge(width = 1)) +
  geom_line(position=position_dodge(width = 1))




tmp1 <- filter(droseu$dia, Supervisor.PI == "Bergland")
tmp2 <- dia[Supervisor.PI==pii]

tmp1 <- as.data.table(tmp1)
tmp2 <- as.data.table(tmp2)
tmp1[,Batch:=as.factor(Batch)]
tmp2[,Batch:=as.factor(Batch)]

setkey(tmp1, Supervisor.PI, Batch, Population, Line)
setkey(tmp2, Supervisor.PI, Batch, Population, Line)



m <- merge(tmp1[,c("Supervisor.PI", "Batch", "Population", "Line", "Prop_Max_Stage9"), with=F],
          tmp2[,c("Supervisor.PI", "Batch", "Population", "Line", "Prop_Max_Stage9"), with=F])
)])

tmp1$Line==tmp2$Line



### glmer
  dia.Bergland <- glmer(Prop_Max_Stage9 ~ Population + (1|Population:Line),
                        weights=dia[Supervisor.PI=="Bergland"]$N,
                        data=dia[Supervisor.PI=="Bergland"],
                        family=binomial())

  
  b <- ranef(dia.Bergland)[[1]]
  b$Line <- rownames(b)
  colnames(b)[1] <- "Bergland"
  
  dia.Flatt <- glmer(Prop_Max_Stage9 ~ Population + (1|Population:Line),
                        weights=dia[Supervisor.PI=="Flatt"]$N,
                        data=dia[Supervisor.PI=="Flatt"],
                        family=binomial())
  
  
  f <- ranef(dia.Flatt)[[1]]
  f$Line <- rownames(f)
  colnames(f)[1] <- "Flatt"
  
  
  
  dia.Schlot <- glmer(Prop_Max_Stage9 ~ Population + (1|Population:Line),
                     weights=dia[Supervisor.PI=="Schlotterer"]$N,
                     data=dia[Supervisor.PI=="Schlotterer"],
                     family=binomial())
  
  s <- ranef(dia.Schlot)[[1]]
  s$Line <- rownames(s)
  colnames(s)[1] <- "Schlotterer"
  
  
  
  
  
  m <- full_join(b,f)
  m <- full_join(m,s)
  

  cor.test(m$Bergland, m$Flatt)
  cor.test(m$Bergland, m$Schlotterer)
  cor.test(m$Flatt, m$Schlotterer)
  
  
  
  library(GGally)
  ggpairs(m[,-2])
  
  
  
  
  
  
  
  
  
  m <- merge(m, s, by = 0)
  
  cor.test(m[,2], m[,3])
  
  
  dia.Bergland <- lmer(asin(sqrt(Prop_Max_Stage9)) ~ 1 + (1|Population),
                        weights=dia[Supervisor.PI=="Bergland"]$N,
                        data=dia[Supervisor.PI=="Bergland"],family=binomial())


  dia.Flatt <- glmer(Prop_Max_Stage9 ~ 1 + (1|Population) + (1|Population:Line),
                        weights=dia[Supervisor.PI=="Flatt"]$N,
                        data=dia[Supervisor.PI=="Flatt"],
                        family=binomial())
