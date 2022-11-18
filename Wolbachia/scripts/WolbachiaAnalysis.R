
# remotes::install_github("bbolker/broom.mixed")
# install.packages("dfoptim")
# install.packages("optimx")
# install.packages("sjPlot")
# install.packages("sjmisc")
# install.packages("see")
library(broom.mixed)
library(afex)
library(tidyverse)
library(sjPlot)
library(sjmisc)
library(car)
library(see)

## adjust contrasts to fit Type - III ANOVAs in R, see here https://rcompanion.org/rcompanion/d_04.html and https://www.r-bloggers.com/2011/03/anova-%E2%80%93-type-iiiiii-ss-explained/
## essentially this uses sum contrasts to compare each group against grand mean.
options(contrasts = c("contr.sum","contr.poly"))

## make list with all phenotypes
Phenotype=c("ChillComa",
  "ColdShockMortality",
  "DevelopmentTime_ETA",
  "Diapause",
  "DryWeight",
  "Fecundity",
  "HeatShock",
  "Lifespan",
  "Pigmentation",
  "Starvation",
  "ThoraxLength",
  "Viability",
  "DevelopmentTime_ETP",
  "WingArea")

## make list with all traits per phenotype
Trait=list(
  list("CCRT_seconds"),
  list("CSM_PropDead"),
  list("DT_EggAdult"),
  #list("Prop_Max_Stage7","Prop_Max_Stage8","Prop_Max_Stage9","Prop_Max_Stage7_asin","Prop_Max_Stage8_asin","Prop_Max_Stage9_asin"),
  list("Prop_Max_Stage8_asin"),
  list("DW_micrograms"),
  list("NumberOfAdultsEclosed"),
  list("TimeDeath_min"),
  list("LSL_AgeAtDeath_days"),
  #list("PercT4","PercT5","PercT6","TotalPerc","PercT4_asin","PercT5_asin","PercT6_asin","TotalPerc_asin"),
  list("TotalPerc_asin"),
  list("AgeAtDeath_hours"),
  list("TL_micrometers"),
  list("ProportionEggtoAdultSurvival"),
  list("DT_EggPupa"),
  list("DevAssym","AvArea"))

## set working directory
setwd("/media/inter/mkapun/projects/DrosEU_PhenotypingWG/")

## Function to plot lineplots with Standard Deviations using ggplot
data_summary <- function(x) {
  m <- mean(x)
  se <- sd(x)
  # se <- sd(x)/sqrt(length(x))
  ymin <- m-se
  ymax <- m+se
  return(c(y=m,ymin=ymin,ymax=ymax))
}

## Function to plot lineplots with Standard Errors using ggplot
data_summary_se <- function(x) {
  m <- mean(x)
  #se <- sd(x)
  se <- sd(x)/sqrt(length(x))
  ymin <- m-se
  ymax <- m+se
  return(c(y=m,ymin=ymin,ymax=ymax))
}

## make output directory
dir.create("Wolbachia/results")

## loop through phenotypes
for (i in seq(1,length(Phenotype),1)){

  ## label Phenotype
  PH=Phenotype[i]
  TR=unlist(Trait[i])

  ## loop through traits
  for (j in seq(1,length(TR),1)){

    ## label Trait
    TRt=unlist(TR[j])

    ## create output directory
    out_dir <- "Wolbachia/results"
    dir.create(file.path(out_dir,PH), showWarnings = F)

    ## start output file for stats
    sink(paste0(out_dir,"/",PH,"/",TRt,".txt"))
    print(paste0("__________________", TRt,"_______________"))

    ## Read Raw Data Table
    DATA=read.table(paste0("Wolbachia/data/",PH,"/",PH,".txt"),
      header=T,
      comment.char = "",
      sep="\t")

    ## Define as factors
    DATA$Wolbachia<-as.factor(DATA$Wolbachia)
    DATA$Country<-as.factor(DATA$Country)
    DATA$Lab<-as.factor(DATA$Supervisor.PI)

    ## Test if Factor "Sex" in Table, if no, add "F", then replace "False" with "F"
    if (!("Sex" %in% colnames(DATA))){
      DATA$Sex <- rep("F",nrow(DATA))
    }
    DATA$Sex[DATA$Sex=="FALSE"]<-"F"

    ## Exclude Finland and Russia, since all lines are Wolbachia positive there.
    DATA<- DATA%>%
    filter(!DATA$Country %in% c("Finland","Russia"))

    ## Exclude lines that have not been phenotyped in all labs
    LABS=length(levels(DATA$Lab))

    tmp<-DATA %>%
    group_by(Line,Lab)%>%
    summarise(N=length(Lab))

    tmp2 <- tmp %>%
    group_by(Line) %>%
    summarise(N=n())

    ## only retain lines that were phenotyped in at least n-1 labs
    tmp3 <- tmp2$Line[tmp2$N>=(LABS-1)]

    DATA<- DATA%>%
      filter(DATA$Line %in% tmp3 )

    ## identify countries where less than 2 lines are either Wol+ or Wol- and remove them to avoid sampling bias.
    tmp<-DATA %>%
    group_by(Country,Wolbachia, Sex,Line)%>%
    summarise(N=length(Line))

    tmp2 <- na.omit(tmp %>%
    group_by(Country,Wolbachia,Sex)%>%
    summarise(N=length(Line)) %>%
    spread(Wolbachia,N))

    ## remove countries with less than 2 lines per infection type.
    tmp3 <- tmp2$Country[tmp2[["-"]]>=2 & tmp2[["+"]]>=2]

    DATA<- DATA %>%
      filter(DATA$Country %in% tmp3 )

    ## only continue if any data left after filtering;-)
    if (nrow(DATA)>0){

      tmp<-DATA %>%
      group_by(Country,Wolbachia, Sex,Line)%>%
      summarise(N=length(Line))

      tmp2 <- na.omit(tmp %>%
      group_by(Country,Wolbachia,Sex)%>%
      summarise(N=length(Line)) %>%
      spread(Wolbachia,N))

      print("______ Table with line counts per country_________")
      print(tmp2)

      ## reorder countries by longitude
      DATA$Country = with(DATA, reorder(Country, Longitude))

      ## Make formula for LMM and include factor Country if it has >1 level
      if (length(unique(as.factor(DATA$Country)))>1) {
        FIXED="Wolbachia * Country"
        TERMS=c("Wolbachia","Country")
        FACET=c("Country")
      } else {
        FIXED="Wolbachia"
        TERMS=c("Wolbachia")
        FACET=c()
      }

      ## Extend formula for LMM and include factor sex if it has >1 level
      if (length(levels(as.factor(DATA$Sex)))>1) {
        FIXED=paste0(FIXED," * Sex")
        TERMS=c(TERMS,"Sex")
      }
      FACET=c(FACET,"Sex")

      ## Set up Formula for Facets
      if (length(FACET) == 2){
        FACET = formula(paste0(FACET[2],"~",FACET[1]))
      } else {
        FACET = formula(paste0("~",FACET[1]))
      }

      ## make new column for plotting
      DATA$Trait <-DATA[[TRt]]

      ## Extend formula for LMM and include factor LAB if it has >1 level and generate plots split by Sex, Country and PI. The error bars are SD
      if (length(levels(as.factor(DATA$Lab)))>1){
        FIXED=paste0(FIXED," + PC.ratio+(1|Lab)+(1|Line:Country)+(1|Batch)")
        TERMS=c("Lab",TERMS)

        ## include line connecting the labs
        DATA.p=ggplot(DATA, aes(x=Lab, y=Trait,col=Wolbachia)) +
        stat_summary(fun.data=data_summary,aes(group=Wolbachia,color=Wolbachia),geom="line",linewidth=0.2,show.legend=F)+
        stat_summary(fun.data=data_summary,aes(colour=Wolbachia,pch=Wolbachia),size=1)+
          theme_bw()+
          facet_grid(FACET,scales="free_y")+
          theme(axis.title.y = element_text(size = 20, angle = 90)) +
          theme(axis.title.x = element_text(size = 20, angle = 00))+
          theme(axis.text=element_text(size=10))+
          theme(legend.text=element_text(size=20))+
          theme(legend.title=element_text(size=20))+
          theme(strip.text =element_text(size=20))+
          ylab(TRt)+
          xlab("Laboratory")

      } else {
        FIXED=paste0(FIXED," + (1|Line:Country)+(1|Batch)")

        DATA.p=ggplot(DATA, aes(x=Lab, y=Trait,col=Wolbachia)) +
        stat_summary(fun.data=data_summary,aes(colour=Wolbachia,pch=Wolbachia),size=1)+
          theme_bw()+
          facet_grid(FACET,scales="free_y")+
          theme(axis.title.y = element_text(size = 20, angle = 90)) +
          theme(axis.title.x = element_text(size = 20, angle = 00))+
          theme(axis.text=element_text(size=10))+
          theme(legend.text=element_text(size=20))+
          theme(legend.title=element_text(size=20))+
          theme(strip.text =element_text(size=20))+
          ylab(TRt)+
          xlab("Laboratory")

      }

      ggsave(paste0(out_dir,"/",PH,"/",TRt,"_raw.pdf"),
        DATA.p,
        width=15,
        height=10)

      ## Make plot averaging across PIs.Note that the Error bars are SE
      DATA.p2=ggplot(DATA, aes(x=Country, y=Trait,col=Wolbachia)) +
      stat_summary(fun.data=data_summary_se,aes(group=Wolbachia,color=Wolbachia),geom="line",linewidth=0.2,show.legend=F)+
      stat_summary(fun.data=data_summary_se,aes(colour=Wolbachia,pch=Wolbachia),size=1)+
        theme_bw()+
        facet_grid(~Sex, scales="free_y")+
        theme(axis.title.y = element_text(size = 20, angle = 90)) +
        theme(axis.title.x = element_text(size = 20, angle = 00))+
        theme(axis.text=element_text(size=10))+
        theme(legend.text=element_text(size=20))+
        theme(legend.title=element_text(size=20))+
        theme(strip.text =element_text(size=20))+
        ylab(TRt)+
        xlab("Country")

      ggsave(paste0(out_dir,"/",PH,"/",TRt,"_raw2.pdf"),
        DATA.p2,
        width=15,
        height=10)

      print("_______________ LMM __________________")
      print(FIXED)

      ## Test for effects of Country and Wolbachia and optionally include Sex as a fixed factor and PI as a random factor.
      f <- formula(paste0(TRt,"~",FIXED))
      Test<-lmer(f,data=DATA)

      ## optionally test how different optimizers affect the model fit
      #aa <- allFit(Test)
      #print(glance(aa) |> select(optimizer, AIC, NLL_rel) |> arrange(NLL_rel))

      print(Anova(Test,
      type=3))

      ## plot model estimates
      F.test<-plot_model(Test,
          type="pred",
          terms=TERMS)

      ggsave(paste0(out_dir,"/",PH,"/",TRt,"_model.pdf"),
        F.test,
        width=15,
        height=15)


      ## Test for the effect of Wolbachia on the Trait variance
      # if (length(levels(as.factor(DATA$Lab)))>1){
      #
      #   print("_______________ LMM Variance__________________")
      #
      #   tmp <- DATA %>%
      #   group_by(Country,Wolbachia,Sex,Lab)%>%
      #   summarise(Mean=mean(Trait),SD=sd(Trait))
      #
      #   if (length(levels(as.factor(DATA$Sex)))>1) {
      #     Test=lm(SD~Country*Wolbachia*Sex,
      #       data=tmp)
      #
      #     print(Anova(Test,
      #       type=3))
      #     } else {
      #       Test=lm(SD~Country*Wolbachia,
      #         data=tmp)
      #
      #       print(Anova(Test,
      #         type=3))
      #     }
      # }

      ## close output file
      sink()
      }
    }
  }
