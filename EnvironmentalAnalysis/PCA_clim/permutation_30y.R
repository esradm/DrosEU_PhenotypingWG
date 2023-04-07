### libraries
library(data.table)
library(ggplot2)
library(patchwork)
library(data.table)
library(patchwork)
library(ggpubr)
library(doMC)
registerDoMC(4)


### load 30d data
#load("/Users/alanbergland/Documents/GitHub/DrosEU_PhenotypingWG/EnvironmentalAnalysis/PCA_clim/TraitsCombined/data/all_traits_30d_WS.RData")
### load 30y data
load("/Users/alanbergland/Documents/GitHub/DrosEU_PhenotypingWG/EnvironmentalAnalysis/PCA_clim/TraitsCombined/data/all_traits_30y_WS.RData")

F9_all_PC <- as.data.table(F9_all_PC)
M9_all_PC <- as.data.table(M9_all_PC)

F9_all_PC[,pheno:="F9_all_PC"]
M9_all_PC[,pheno:="M9_all_PC"]

phenoPC <- rbind(F9_all_PC, M9_all_PC)

#for 30d
#m <- as.data.table(merge(d30_pca_results, phenoPC, by="Country"))
#for 30y
m <- as.data.table(merge(climate_data, phenoPC, by="Country"))

nPerm <- 1000

o <-
  foreach(pheno.i=c("F9_all_PC", "M9_all_PC"), .combine="rbind")%do%{
    foreach(pheno.pc.i=paste("Dim", c(1:5), sep="."), .combine="rbind")%do%{
      foreach(env.pc.i=c("PC1_clim", "PC2_clim"), .combine="rbind")%do%{
        foreach(perm=0:nPerm, .combine="rbind")%dopar%{
          message(paste(pheno.i, pheno.pc.i, env.pc.i, perm, sep=" / "))
          # pheno.i <- "F9_all_PC"; pheno.pc.i <- "Dim.1"; env.pc.i <- "PC1_clim"; perm=0
          
          tmp <- m[pheno==pheno.i][,c(pheno.pc.i, env.pc.i, "Population.x"), with=F]
          
          setnames(tmp, c(pheno.pc.i, env.pc.i, "Population.x"), c("pheno_pc", "env_pc", "pop"))
          
          tmp.ag <- tmp[,list(env_pc=mean(env_pc)), list(pop)]
          
          
          # if(perm==0) {
          #   tmp[,pheno_pc_perm:=pheno_pc]
          # } else {
          #   tmp[,pheno_pc_perm:=sample(pheno_pc, replace=F)]
          # }
          
          if(perm==0) {
            tmp.ag[,env_pc_perm:=env_pc]
          } else {
            tmp.ag[,env_pc_perm:=sample(env_pc, replace=F)]
          }
          
          tmp <- merge(tmp, tmp.ag, by="pop")
          
          t1 <- lm(pheno_pc~env_pc_perm, tmp)
          
          data.table(traitPC=pheno.pc.i, envPC=env.pc.i, perm=perm, pheno=pheno.i,
                     p=summary(t1)$coef[2,4], r2=summary(t1)$r.squared)
          
          
        }
      }
    }
  }



o.ag <- o[,list(p_thr=quantile(p, .05), r2_thr=quantile(r2, .95)), list(traitPC, envPC, pheno, perm=as.factor(perm!=0))]

permPlot <- ggplot() +
  geom_point(data=o.ag[perm==F], aes(x=-log10(p_thr), y=r2_thr, color=pheno)) +
  geom_hline(data=o.ag[perm==T], aes(yintercept=r2_thr,         color=pheno)) +
  geom_vline(data=o.ag[perm==T], aes(xintercept=-log10(p_thr),  color=pheno)) +
  ylab("R^2") + xlab("-log10(p)") + facet_grid(traitPC~envPC)

ggsave(permPlot, file="/Users/alanbergland/Documents/GitHub/DrosEU_PhenotypingWG/EnvironmentalAnalysis/PCA_clim/env_30y_trait_pc.pdf")