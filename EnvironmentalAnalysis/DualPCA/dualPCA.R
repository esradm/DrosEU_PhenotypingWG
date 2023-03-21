### libraries
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(patchwork)
  library(ggpubr)

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

### general layout
  layout <- "
  ABBB
  ABBB
  ABBB
  #CCC"


### FMax-30d
    flip_Fmax_30d <- -1
    Fmax_30d_traitPC_load <- ggplot(data=FMax_traitPC) +
                  geom_segment( aes(x=0, xend=Dim.2, y=ordDim2, yend=ordDim2), arrow = arrow(length = unit(0.15, "cm"))) +
                  geom_text(data=FMax_traitPC[Dim.2<0], aes(y=ordDim2, x=.05, label=trait), size=2.5, angle=90, hjust=0) +
                  geom_text(data=FMax_traitPC[Dim.2>0], aes(y=ordDim2, x=-.05, label=trait), size=2.5, angle=90, hjust=1) +
                  coord_flip() + xlim(-.5, 1) + ylim(-0, 13.5) +
                  xlab("Phenotype PC2 loadings") + ylab("") +
                  theme_bw()

    Fmax_30d_phenoPC_load <- ggplot(data=d30_envPC) +
                  geom_segment(aes(x=0, xend=flip_Fmax_30d*Dim.2, y=ordDim2, yend=ordDim2), arrow = arrow(length = unit(0.15, "cm"))) +
                  geom_text(data=d30_envPC[Dim.2<0], aes(y=ordDim2, x=flip_Fmax_30d*.05, label=trait), size=2.5, hjust=0) +
                  geom_text(data=d30_envPC[Dim.2>0], aes(y=ordDim2, x=flip_Fmax_30d*-.05, label=trait), size=2.5, hjust=1) +
                  xlim(-1, 1)  + ylim(-0, 15.5) +
                  xlab("Environment PC2 loadings") + ylab("") +
                  theme_bw()

    Fmax_30d_trait_env_plot <- ggplot(data=d30_F9max_data, aes(x=flip_Fmax_30d*PC2_clim, y=PC2_F9max, color=Country)) +
      geom_point() +
      theme_bw() +
      stat_smooth(method = "lm")

    Fmax_d30_mega <-
    Fmax_30d_traitPC_load + Fmax_30d_trait_env_plot + Fmax_30d_phenoPC_load +
    plot_layout(design=layout) +
    plot_annotation(title = paste("Day30 environmental + FMax pheno + flip:", flip_Fmax_30d, sep=" "))

### FMax-30y
    flip_Fmax_30y <- 1
    Fmax_30y_traitPC_load <- ggplot(data=FMax_traitPC) +
                  geom_segment( aes(x=0, xend=Dim.2, y=ordDim2, yend=ordDim2), arrow = arrow(length = unit(0.15, "cm"))) +
                  geom_text(data=FMax_traitPC[Dim.2<0], aes(y=ordDim2, x=.05, label=trait), size=2.5, angle=90, hjust=0) +
                  geom_text(data=FMax_traitPC[Dim.2>0], aes(y=ordDim2, x=-.05, label=trait), size=2.5, angle=90, hjust=1) +
                  coord_flip() + xlim(-.5, 1) + ylim(-0, 13.5) +
                  xlab("Phenotype PC2 loadings") + ylab("") +
                  theme_bw()
    Fmax_30y_phenoPC_load <- ggplot(data=y30_envPC) +
                  geom_segment(aes(x=0, xend=flip_Fmax_30y*Dim.2, y=ordDim2, yend=ordDim2), arrow = arrow(length = unit(0.15, "cm"))) +
                  geom_text(data=y30_envPC[Dim.2<0], aes(y=ordDim2, x=flip_Fmax_30y*.05, label=trait), size=2.5, hjust=0) +
                  geom_text(data=y30_envPC[Dim.2>0], aes(y=ordDim2, x=flip_Fmax_30y*-.05, label=trait), size=2.5, hjust=1) +
                  xlim(-1, 1)  + ylim(-0, 15.5) +
                  xlab("Environment PC2 loadings") + ylab("") +
                  theme_bw()

    Fmax_30y_trait_env_plot <- ggplot(data=y30_F9max_data, aes(x=flip_Fmax_30y*PC2_clim, y=PC2_F9max, color=Country)) +
      geom_point() +
      theme_bw() +
      stat_smooth(method = "lm")

    Fmax_y30_mega <-
    Fmax_30y_traitPC_load + Fmax_30y_trait_env_plot + Fmax_30y_phenoPC_load +
    plot_layout(design=layout) +
    plot_annotation(
    title = paste("Year30 environmental + FMax pheno + flip:", flip_Fmax_30y, sep=" "))

    ggsave(Fmax_d30_mega, file="/Users/alanbergland/Documents/GitHub/DrosEU_PhenotypingWG/EnvironmentalAnalysis/DualPCA/DualPCA_Fmax_d30.pdf", h=8, w=8)
    ggsave(Fmax_y30_mega, file="/Users/alanbergland/Documents/GitHub/DrosEU_PhenotypingWG/EnvironmentalAnalysis/DualPCA/DualPCA_Fmax_y30.pdf", h=8, w=8)


### M9-30d
    flip_M9_30d <- -1
    M9_30d_traitPC_load <- ggplot(data=M9_traitPC) +
                  geom_segment( aes(x=0, xend=Dim.2, y=ordDim2, yend=ordDim2), arrow = arrow(length = unit(0.15, "cm"))) +
                  geom_text(data=M9_traitPC[Dim.2<0], aes(y=ordDim2, x=.05, label=trait), size=2.5, angle=90, hjust=0) +
                  geom_text(data=M9_traitPC[Dim.2>0], aes(y=ordDim2, x=-.05, label=trait), size=2.5, angle=90, hjust=1) +
                  coord_flip() + xlim(-.5, 1) + ylim(-0, 13.5) +
                  xlab("Phenotype PC2 loadings") + ylab("") +
                  theme_bw()
    M9_30d_phenoPC_load <- ggplot(data=d30_envPC) +
                  geom_segment(aes(x=0, xend=flip_M9_30d*Dim.2, y=ordDim2, yend=ordDim2), arrow = arrow(length = unit(0.15, "cm"))) +
                  geom_text(data=d30_envPC[Dim.2<0], aes(y=ordDim2, x=flip_M9_30d*.05, label=trait), size=2.5, hjust=0) +
                  geom_text(data=d30_envPC[Dim.2>0], aes(y=ordDim2, x=flip_M9_30d*-.05, label=trait), size=2.5, hjust=1) +
                  xlim(-1, 1)  + ylim(-0, 15.5) +
                  xlab("Environment PC2 loadings") + ylab("") +
                  theme_bw()

    M9_30d_trait_env_plot <- ggplot(data=d30_M9_data, aes(x=flip_M9_30d*PC2_clim, y=PC2_M9, color=Country)) +
      geom_point() +
      theme_bw() +
      stat_smooth(method = "lm")


    M9_d30_mega <-
    M9_30d_traitPC_load + M9_30d_trait_env_plot + M9_30d_phenoPC_load +
    plot_layout(design=layout) +
    plot_annotation(
    title = paste("Day30 environmental + M9 pheno + flip:", flip_M9_30d, sep=" "))

### M9-30y
    flip_M9_30y <- 1
    M9_30y_traitPC_load <- ggplot(data=M9_traitPC) +
                  geom_segment( aes(x=0, xend=Dim.2, y=ordDim2, yend=ordDim2), arrow = arrow(length = unit(0.15, "cm"))) +
                  geom_text(data=M9_traitPC[Dim.2<0], aes(y=ordDim2, x=.05, label=trait), size=2.5, angle=90, hjust=0) +
                  geom_text(data=M9_traitPC[Dim.2>0], aes(y=ordDim2, x=-.05, label=trait), size=2.5, angle=90, hjust=1) +
                  coord_flip() + xlim(-.5, 1) + ylim(-0, 13.5) +
                  xlab("Phenotype PC2 loadings") + ylab("") +
                  theme_bw()
    M9_30y_phenoPC_load <- ggplot(data=y30_envPC) +
                  geom_segment(aes(x=0, xend=flip_M9_30y*Dim.2, y=ordDim2, yend=ordDim2), arrow = arrow(length = unit(0.15, "cm"))) +
                  geom_text(data=y30_envPC[Dim.2<0], aes(y=ordDim2, x=flip_M9_30y*.05, label=trait), size=2.5, hjust=0) +
                  geom_text(data=y30_envPC[Dim.2>0], aes(y=ordDim2, x=flip_M9_30y*-.05, label=trait), size=2.5, hjust=1) +
                  xlim(-1, 1)  + ylim(-0, 15.5) +
                  xlab("Environment PC2 loadings") + ylab("") +
                  theme_bw()

    M9_30y_trait_env_plot <- ggplot(data=y30_M9_data, aes(x=flip_M9_30y*PC2_clim, y=PC2_M9, color=Country)) +
      geom_point() +
      theme_bw() +
      stat_smooth(method = "lm")

  M9_y30_mega <-
  M9_30y_traitPC_load + M9_30y_trait_env_plot + M9_30y_phenoPC_load +
  plot_layout(design=layout) +
  plot_annotation(
  title = paste("Year30 environmental + M9 pheno + flip:", flip_M9_30y, sep=" "))


### save mega plots
  ggsave(Fmax_d30_mega, file="/Users/alanbergland/Documents/GitHub/DrosEU_PhenotypingWG/EnvironmentalAnalysis/DualPCA/DualPCA_Fmax_d30.pdf", h=8, w=8)
  ggsave(Fmax_y30_mega, file="/Users/alanbergland/Documents/GitHub/DrosEU_PhenotypingWG/EnvironmentalAnalysis/DualPCA/DualPCA_Fmax_y30.pdf", h=8, w=8)
  ggsave(M9_d30_mega,   file="/Users/alanbergland/Documents/GitHub/DrosEU_PhenotypingWG/EnvironmentalAnalysis/DualPCA/DualPCA_M9_d30.pdf", h=8, w=8)
  ggsave(M9_y30_mega,   file="/Users/alanbergland/Documents/GitHub/DrosEU_PhenotypingWG/EnvironmentalAnalysis/DualPCA/DualPCA_M9_y30.pdf", h=8, w=8)


  ggsave(Fmax_d30_mega, file="/Users/alanbergland/Documents/GitHub/DrosEU_PhenotypingWG/EnvironmentalAnalysis/DualPCA/DualPCA_Fmax_d30.png", h=8, w=8)
  ggsave(Fmax_y30_mega, file="/Users/alanbergland/Documents/GitHub/DrosEU_PhenotypingWG/EnvironmentalAnalysis/DualPCA/DualPCA_Fmax_y30.png", h=8, w=8)
  ggsave(M9_d30_mega,     file="/Users/alanbergland/Documents/GitHub/DrosEU_PhenotypingWG/EnvironmentalAnalysis/DualPCA/DualPCA_M9_d30.png", h=8, w=8)
  ggsave(M9_y30_mega,     file="/Users/alanbergland/Documents/GitHub/DrosEU_PhenotypingWG/EnvironmentalAnalysis/DualPCA/DualPCA_M9_y30.png", h=8, w=8)









F9Max_D30_mega + F9Max_D30_mega
  ggsave("/Users/alanbergland/Documents/GitHub/DrosEU_PhenotypingWG/EnvironmentalAnalysis/test_PCA_plot.pdf", h=8, w=8)
