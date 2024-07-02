### libraries
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(patchwork)
  library(ggpubr)
  library(doMC)
  registerDoMC(4)

### F9max
  ### Load FMax Trait PC
    load("/Users/alanbergland/Documents/GitHub/DrosEU_PhenotypingWG/EnvironmentalAnalysis/PCA_Results/Fmax_drosEU.RData")

    ### trait PCA loadings
      FMax_traitPC <- as.data.table(pca3$var$coord)
      FMax_traitPC[,trait:=rownames(pca3$var$coord)]

      FMax_traitPC[,ordDim1:=rank(Dim.1)]
      FMax_traitPC[,ordDim2:=rank(Dim.2)]

  ### Load M9 Trait PC
    load("/Users/alanbergland/Documents/GitHub/DrosEU_PhenotypingWG/EnvironmentalAnalysis/PCA_Results/M9_drosEU.RData")

    ### trait PCA loadings
      M9_traitPC <- as.data.table(pca1$var$coord)
      M9_traitPC[,trait:=rownames(pca1$var$coord)]

      M9_traitPC[,ordDim1:=rank(Dim.1)]
      M9_traitPC[,ordDim2:=rank(Dim.2)]

  ### Load Day 30
    load("/Users/alanbergland/Documents/GitHub/DrosEU_PhenotypingWG/EnvironmentalAnalysis/PCA_clim/TraitsCombined/data/all_traits_30d_WS.RData")
    d30_F9max_data <- F9max_data
    d30_M9_data <- M9_data

    ### env PCA loadings
      d30_envPC <- as.data.table(d30_pca$var$coord)
      d30_envPC[,trait:=rownames(d30_pca$var$coord)]
      d30_envPC[,ordDim1:=rank(Dim.1)]
      d30_envPC[,ordDim2:=rank(Dim.2)]

  ### load y30
      load("/Users/alanbergland/Documents/GitHub/DrosEU_PhenotypingWG/EnvironmentalAnalysis/PCA_clim/TraitsCombined/data/all_traits_30y_WS.RData")
      y30_F9max_data <- F9max_data
      y30_M9_data <- M9_data

    ### env PCA loadings
      y30_envPC <- as.data.table(all_ann_clim_pca$var$coord)
      y30_envPC[,trait:=rownames(all_ann_clim_pca$var$coord)]
      y30_envPC[,ordDim1:=rank(Dim.1)]
      y30_envPC[,ordDim2:=rank(Dim.2)]



### combinations & permutations
  d30_F9max_data <- as.data.table(d30_F9max_data)
  o <- foreach(trait_pc=c("PC1_F9max", "PC2_F9max", "PC3_F9max"), .combine="rbind")%do%{
    foreach(env_pc=c("PC1_clim", "PC2_clim"), .combine="rbind")%do%{
      # trait_pc="PC2_F9max"; env_pc="PC2_clim"
      tmp <- d30_F9max_data[,c("Location", "Population", "Line", trait_pc, env_pc), with=F]
      setnames(tmp, c(trait_pc, env_pc), c("trait_pc_value", "env_pc_value"))

      foreach(perm=0:100, .combine="rbind")%dopar%{
        message(paste(trait_pc, env_pc, perm, sep=" / "))
        tmp.ag <- tmp[,list(env_pc_value=unique(env_pc_value)), list(Location)]
        if(perm==0) {
          tmp.ag[,env_pc_value_perm:=env_pc_value]
        } else {
          tmp.ag[,env_pc_value_perm:=sample(env_pc_value, replace=F)]
        }
        tmp <- merge(tmp, tmp.ag, by="Location")

        t1 <- lm(trait_pc_value~env_pc_value_perm, tmp)
        data.table(traitPC=trait_pc, envPC=env_pc, perm=perm,
                   p=summary(t1)$coef[2,4], r2=summary(t1)$r.squared)
      }

    }
  }

  o.ag <- o[,list(p_thr=quantile(p, .05), r2_thr=quantile(r2, .95)), list(traitPC, envPC, perm=as.factor(perm!=0))]


  ggplot(data=o, aes(x=-log10(p), y=r2, color=as.factor(perm==0))) +
  facet_grid(traitPC~envPC) + geom_point()

  d30_M9_data <- as.data.table(d30_M9_data)
  o_male <- foreach(trait_pc=c("PC1_M9", "PC2_M9", "PC3_M9"), .combine="rbind")%do%{
    foreach(env_pc=c("PC1_clim", "PC2_clim"), .combine="rbind")%do%{
      # trait_pc="PC2_F9max"; env_pc="PC2_clim"
      tmp <- d30_M9_data[,c("Location", "Population", "Line", trait_pc, env_pc), with=F]
      setnames(tmp, c(trait_pc, env_pc), c("trait_pc_value", "env_pc_value"))

      foreach(perm=0:100, .combine="rbind")%dopar%{
        message(paste(trait_pc, env_pc, perm, sep=" / "))
        tmp.ag <- tmp[,list(env_pc_value=unique(env_pc_value)), list(Location)]
        if(perm==0) {
          tmp.ag[,env_pc_value_perm:=env_pc_value]
        } else {
          tmp.ag[,env_pc_value_perm:=sample(env_pc_value, replace=F)]
        }
        tmp <- merge(tmp, tmp.ag, by="Location")

        t1 <- lm(trait_pc_value~env_pc_value_perm, tmp)
        data.table(traitPC=trait_pc, envPC=env_pc, perm=perm,
                   p=summary(t1)$coef[2,4], r2=summary(t1)$r.squared)
      }

    }
  }

  o_male.ag <- o_male[,list(p_thr=quantile(p, .05), r2_thr=quantile(r2, .95)), list(traitPC, envPC, perm=as.factor(perm!=0))]




  f9max_plot <- ggplot() +
  geom_point(data=o.ag[perm==F], aes(x=-log10(p_thr), y=r2_thr, color=interaction(traitPC, envPC))) +
  geom_hline(data=o.ag[perm==T], aes(yintercept=r2_thr, color=interaction(traitPC, envPC))) +
  geom_vline(data=o.ag[perm==T], aes(xintercept=-log10(p_thr), color=interaction(traitPC, envPC))) +
  ylab("R^2") + xlab("-log10(p)") + ggtitle("F9max")



  m9plot <- ggplot() +
  geom_point(data=o_male.ag[perm==F], aes(x=-log10(p_thr), y=r2_thr, color=interaction(traitPC, envPC))) +
  geom_hline(data=o_male.ag[perm==T], aes(yintercept=r2_thr, color=interaction(traitPC, envPC))) +
  geom_vline(data=o_male.ag[perm==T], aes(xintercept=-log10(p_thr), color=interaction(traitPC, envPC))) +
  ylab("R^2") + xlab("-log10(p)") + ggtitle("M9")

f9max_plot + m9plot + plot_layout(guides = "collect")


  str(d30_F9max_data)
  t1 <- lm(PC2_F9max~PC2_clim, data=d30_F9max_data)
  summary(t1)
