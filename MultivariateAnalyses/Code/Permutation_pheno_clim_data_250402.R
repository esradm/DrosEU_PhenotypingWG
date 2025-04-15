rm(list=ls(all=TRUE))

### libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(cowplot)
#library(doMC)
#registerDoMC(4)

library(doSNOW)
#registerDoSNOW(1)

#cl <- makeCluster(4, type="SOCK")
#registerDoSNOW(cl)
#m <- matrix(rnorm(9), 3, 3)
#foreach(i=1:nrow(m), .combine=rbind) %dopar% (m[i,] / mean(m[i,]))
#stopCluster(cl)

### load 30y data
load("C://Users//ewanh//Dropbox//Barcelona_IBE//DrosEU//all_traits_30y_WS.RData")

climate_data2<-rename(climate_data,
                      PC1_clim = Dim.1,
                      PC2_clim = Dim.2)

F13_trait_PC <- as.data.table(F9maxP_all_PC)
F9_trait_PC <- as.data.table(F9_all_PC)
M9_trait_PC <- as.data.table(M9_all_PC)

F13_trait_PC[,pheno:="F13"]
F9_trait_PC[,pheno:="F9"]
M9_trait_PC[,pheno:="M9"]

F13_trait_clim_y <- as.data.table(merge(climate_data2, F13_trait_PC, by="Country"))
F9_trait_clim_y <- as.data.table(merge(climate_data2, F9_trait_PC, by="Country"))
M9_trait_clim_y <- as.data.table(merge(climate_data2, M9_trait_PC, by="Country"))

saveRDS(F13_trait_clim_y, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/F13_trait_clim_y.rds")
saveRDS(F9_trait_clim_y, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/F9_trait_clim_y.rds")
saveRDS(M9_trait_clim_y, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/M9_trait_clim_y.rds")

##################################
###
### Remove all and load 30 day data so that you can create the 3 objects for day data
###

rm(list=ls(all=TRUE))

### load 30d data
load("C://Users//ewanh//Dropbox//Barcelona_IBE//DrosEU//all_traits_30d_WS.RData")

## No need to rename variables - already good names
# climate_data2<-rename(climate_data,
#                      CDim.1 = Dim.1,
#                      CDim.2 = Dim.2)

climate_data2<-climate_data

F13_trait_PC <- as.data.table(F9maxP_all_PC)
F9_trait_PC <- as.data.table(F9_all_PC)
M9_trait_PC <- as.data.table(M9_all_PC)

F13_trait_PC[,pheno:="F13"]
F9_trait_PC[,pheno:="F9"]
M9_trait_PC[,pheno:="M9"]

F13_trait_clim_d <- as.data.table(merge(climate_data2, F13_trait_PC, by="Country"))
F9_trait_clim_d <- as.data.table(merge(climate_data2, F9_trait_PC, by="Country"))
M9_trait_clim_d <- as.data.table(merge(climate_data2, M9_trait_PC, by="Country"))

saveRDS(F13_trait_clim_d, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/F13_trait_clim_d.rds")
saveRDS(F9_trait_clim_d, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/F9_trait_clim_d.rds")
saveRDS(M9_trait_clim_d, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/M9_trait_clim_d.rds")

##################################
rm(list=ls(all=TRUE))

pca_F13_y <- readRDS("C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/F13_trait_clim_y.rds")
pca_F9_y <- readRDS("C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/F9_trait_clim_y.rds")
pca_M9_y <- readRDS("C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/M9_trait_clim_y.rds")
pca_F13_d <- readRDS("C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/F13_trait_clim_d.rds")
pca_F9_d <- readRDS("C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/F9_trait_clim_d.rds")
pca_M9_d <- readRDS("C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/M9_trait_clim_d.rds")

##################################

m <- pca_M9_d
nPerm <- 1000

o <-
  foreach(pheno.i=c("M9"), .combine="rbind")%do%{
    foreach(pheno.pc.i=paste("Dim", c(1:3), sep="."), .combine="rbind")%do%{
      foreach(env.pc.i=c("PC1_clim", "PC2_clim"), .combine="rbind")%do%{
        foreach(perm=0:nPerm, .combine="rbind")%do%{
          message(paste(pheno.i, pheno.pc.i, env.pc.i, perm, sep=" / "))
          # pheno.i <- "F9_all_PC"; pheno.pc.i <- "Dim.1"; env.pc.i <- "PC1_clim"; perm=0
          
          tmp <- m[pheno==pheno.i][,c(pheno.pc.i, env.pc.i, "Population.x"), with=F]
          
          setnames(tmp, c(pheno.pc.i, env.pc.i, "Population.x"), c("pheno_pc", "env_pc", "pop"))
          
          tmp.ag <- tmp[,list(env_pc=mean(env_pc)), list(pop)]
          
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



o_F13y <- o
o_F9y <- o
o_M9y <- o

o_F13d <- o
o_F9d <- o
o_M9d <- o

o.ag_F13y <- o_F13y[,list(p_thr=quantile(p, .05), r2_thr=quantile(r2, .95)), list(traitPC, envPC, pheno, perm=as.factor(perm!=0))]
o.ag_F9y <- o_F9y[,list(p_thr=quantile(p, .05), r2_thr=quantile(r2, .95)), list(traitPC, envPC, pheno, perm=as.factor(perm!=0))]
o.ag_M9y <- o_M9y[,list(p_thr=quantile(p, .05), r2_thr=quantile(r2, .95)), list(traitPC, envPC, pheno, perm=as.factor(perm!=0))]

o.ag_F13d <- o_F13d[,list(p_thr=quantile(p, .05), r2_thr=quantile(r2, .95)), list(traitPC, envPC, pheno, perm=as.factor(perm!=0))]
o.ag_F9d <- o_F9d[,list(p_thr=quantile(p, .05), r2_thr=quantile(r2, .95)), list(traitPC, envPC, pheno, perm=as.factor(perm!=0))]
o.ag_M9d <- o_M9d[,list(p_thr=quantile(p, .05), r2_thr=quantile(r2, .95)), list(traitPC, envPC, pheno, perm=as.factor(perm!=0))]

###############

emp_pr_F13_y <- o_F13y[,list(pr_p=mean(p[perm==0]<p[perm!=0])), list(traitPC, envPC, pheno)]
emp_pr_F9_y <- o_F9y[,list(pr_p=mean(p[perm==0]<p[perm!=0])), list(traitPC, envPC, pheno)]
emp_pr_M9_y <- o_M9y[,list(pr_p=mean(p[perm==0]<p[perm!=0])), list(traitPC, envPC, pheno)]

emp_pr_F13_d <- o_F13d[,list(pr_p=mean(p[perm==0]<p[perm!=0])), list(traitPC, envPC, pheno)]
emp_pr_F9_d <- o_F9d[,list(pr_p=mean(p[perm==0]<p[perm!=0])), list(traitPC, envPC, pheno)]
emp_pr_M9_d <- o_M9d[,list(pr_p=mean(p[perm==0]<p[perm!=0])), list(traitPC, envPC, pheno)]

###############

emp_pr_F13_y[order(pr_p)]
emp_pr_F9_y[order(pr_p)]
emp_pr_M9_y[order(pr_p)]

emp_pr_F13_d[order(pr_p)]
emp_pr_F9_d[order(pr_p)]
emp_pr_M9_d[order(pr_p)]

#######################

"#cab2d6"
"#6a3d9a"
"#fdbf6f"
"#ff7f00"
"#fb9a99"
"#e31a1c"

###########################

permPlot1 <- ggplot(data=o_F13d[perm!=0], aes(x=-log10(p), y=r2)) +
  geom_hline(data=o.ag_F13d[perm==T], aes(yintercept=r2_thr, color=pheno), colour = "#EF4949") +
  geom_vline(data=o.ag_F13d[perm==T], aes(xintercept=-log10(p_thr), color=pheno), colour = "#EF4949") +
  geom_point(data=o_F13d[perm!=0], aes(colour= pheno)) + 
  scale_colour_manual(values = "grey85") +
  facet_grid(envPC~traitPC) +
  geom_point(data=o_F13d[perm==0], aes(x=-log10(p), y=r2, color=pheno), cex = 2.6, shape=21, fill="#EF4949") +
  ylab(expression(paste("R"^"2"))) + xlab(expression(paste("-log" [10]," p"))) + 
  ggtitle("F13 - 30 day") +  
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none")


permPlot2 <- ggplot(data=o_F9d[perm!=0], aes(x=-log10(p), y=r2)) +
  geom_hline(data=o.ag_F9d[perm==T], aes(yintercept=r2_thr, color=pheno), colour = "#FF9B39") +
  geom_vline(data=o.ag_F9d[perm==T], aes(xintercept=-log10(p_thr), color=pheno), colour = "#FF9B39") +
  geom_point(data=o_F9d[perm!=0], aes(colour= pheno)) + 
  scale_colour_manual(values = "grey85") +
  facet_grid(envPC~traitPC) +
  geom_point(data=o_F9d[perm==0], aes(x=-log10(p), y=r2, color=pheno), cex = 2.6, shape=21, fill="#FF9B39") +
  ylab(expression(paste("R"^"2"))) + xlab(expression(paste("-log" [10]," p"))) + 
  ggtitle("F9 - 30 day") +  
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none")


permPlot3 <- ggplot(data=o_M9d[perm!=0], aes(x=-log10(p), y=r2)) +
  geom_hline(data=o.ag_M9d[perm==T], aes(yintercept=r2_thr, color=pheno), colour = "#855CB1") +
  geom_vline(data=o.ag_M9d[perm==T], aes(xintercept=-log10(p_thr), color=pheno), colour = "#855CB1") +
  geom_point(data=o_M9d[perm!=0], aes(colour= pheno)) + 
  scale_colour_manual(values = "grey85") +
  facet_grid(envPC~traitPC) +
  geom_point(data=o_M9d[perm==0], aes(x=-log10(p), y=r2, color=pheno), cex = 2.6, shape=21, fill="#855CB1") +
  ylab(expression(paste("R"^"2"))) + xlab(expression(paste("-log" [10]," p"))) + 
  ggtitle("M9 - 30 day") +  
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none")


permPlot4 <- ggplot(data=o_F13y[perm!=0], aes(x=-log10(p), y=r2)) +
  geom_hline(data=o.ag_F13y[perm==T], aes(yintercept=r2_thr, color=pheno), colour = "#B20C0C") +
  geom_vline(data=o.ag_F13y[perm==T], aes(xintercept=-log10(p_thr), color=pheno), colour = "#B20C0C") +
  geom_point(data=o_F13y[perm!=0], aes(colour= pheno)) + 
  scale_colour_manual(values = "grey85") +
  facet_grid(envPC~traitPC) +
  geom_point(data=o_F13y[perm==0], aes(x=-log10(p), y=r2, color=pheno), cex = 2.6, shape=21, fill="#B20C0C") +
  ylab(expression(paste("R"^"2"))) + xlab(expression(paste("-log" [10]," p"))) + 
  ggtitle("F13 - 30 year") +  
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none")


permPlot5 <- ggplot(data=o_F9y[perm!=0], aes(x=-log10(p), y=r2)) +
  geom_hline(data=o.ag_F9y[perm==T], aes(yintercept=r2_thr, color=pheno), colour = "#C66200") +
  geom_vline(data=o.ag_F9y[perm==T], aes(xintercept=-log10(p_thr), color=pheno), colour = "#C66200") +
  geom_point(data=o_F9y[perm!=0], aes(colour= pheno)) + 
  scale_colour_manual(values = "grey85") +
  facet_grid(envPC~traitPC) +
  geom_point(data=o_F9y[perm==0], aes(x=-log10(p), y=r2, color=pheno), cex = 2.6, shape=21, fill="#C66200") +
  ylab(expression(paste("R"^"2"))) + xlab(expression(paste("-log" [10]," p"))) + 
  ggtitle("F9 - 30 year") +  
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none")


permPlot6 <- ggplot(data=o_M9y[perm!=0], aes(x=-log10(p), y=r2)) +
  geom_hline(data=o.ag_M9y[perm==T], aes(yintercept=r2_thr, color=pheno), colour = "#56238D") +
  geom_vline(data=o.ag_M9y[perm==T], aes(xintercept=-log10(p_thr), color=pheno), colour = "#56238D") +
  geom_point(data=o_M9y[perm!=0], aes(colour= pheno)) + 
  scale_colour_manual(values = "grey85") +
  facet_grid(envPC~traitPC) +
  geom_point(data=o_M9y[perm==0], aes(x=-log10(p), y=r2, color=pheno), cex = 2.6, shape=21, fill="#56238D") +
  ylab(expression(paste("R"^"2"))) + xlab(expression(paste("-log" [10]," p"))) + 
  ggtitle("M9 - 30 year") +  
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none")

################

plot_grid(permPlot1,
          permPlot2,
          permPlot3,
          permPlot4,
          permPlot5,
          permPlot6, nrow = 2)

saveRDS(o_F13y, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/o_F13y.rds")
saveRDS(o_F9y, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/o_F9y.rds")
saveRDS(o_M9y, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/o_M9y.rds")
saveRDS(o_F13d, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/o_F13d.rds")
saveRDS(o_F9d, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/o_F9d.rds")
saveRDS(o_M9d, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/o_M9d.rds")

saveRDS(emp_pr_F13_y, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/emp_pr_F13_y.rds")
saveRDS(emp_pr_F9_y, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/emp_pr_F9_y.rds")
saveRDS(emp_pr_M9_y, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/emp_pr_M9_y.rds")
saveRDS(emp_pr_F13_d, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/emp_pr_F13_d.rds")
saveRDS(emp_pr_F9_d, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/emp_pr_F9_d.rds")
saveRDS(emp_pr_M9_d, file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/emp_pr_M9_d.rds")

save.image(file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/drosEU_permutation_results.RData")

load(file = "C:/Users/ewanh/Dropbox/Barcelona_IBE/DrosEU/June2024_final/drosEU_permutation_results.RData")
