library(tidyverse)
library(Seurat)
library(ggplot2)
library(reshape2)
library(ggpubr)

load(file = "comb.Robj")

#function to summarize data
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

#funtion to plot individual data points
plot_effects <- function(data,comparision_group,comparisions,values,sample_id,
                         patient_col_dict) {
  geneset <- values
  
  Data = data[,c(sample_id,comparision_group,geneset)]
  
  plot <- Data
  
  plot_summary <- data_summary(plot, varname=geneset, 
                               groupnames=c(comparision_group, sample_id))
  
  plot_summary$geneset <- plot_summary[[geneset]]
  plot_summary$comparision_group <- plot_summary[[comparision_group]]
  
  p1 <- ggplot(data = plot_summary, aes_string(x=comparision_group, y= geneset), fill = sample_id) +
    geom_pointrange(aes(ymin=geneset-sd, 
                        ymax=geneset+sd, color = sample_id),
                    position = position_dodge(0.5), width = 0.2)+ 
    scale_color_manual(values = patient_col_dict) +
    stat_compare_means(comparisons = comparisions,label = "p.signif", paired = F)+
    theme_light() +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    xlab("") + ylab(paste(geneset, " scores")) +
    theme(strip.background =element_rect(fill="white"))+
    theme(strip.text = element_text(colour = 'black'), ) +
    theme(legend.position="none") 
  
  return(p1)
}

#adding score for each gene separately, ALDH1 has all isotypes
Seurat_obj <- comb
DefaultAssay(Seurat_obj) <- "SCT"
Seurat_obj <- AddModuleScore(Seurat_obj, features = list(c("ALDH1A1","ALDH1A2","ALDH1A3")),
                             name = "ALDH1")
Seurat_obj <- AddModuleScore(Seurat_obj, features = list(c("VIM")),
                             name = "VIM")
Seurat_obj <- AddModuleScore(Seurat_obj, features = list(c("EPCAM")),
                             name = "EPCAM")
Seurat_obj <- AddModuleScore(Seurat_obj, features = list(c("CDH1")),
                             name = "CDH1")
Seurat_obj <- AddModuleScore(Seurat_obj, features = list(c("CD24")),
                             name = "CD24")
Seurat_obj <- AddModuleScore(Seurat_obj, features = list(c("CD44")),
                             name = "CD44")

#MET like CSC score -VIM + EPCAM +CDH1 +ALDH1
Seurat_obj@meta.data$MET_like_CSC <- -1*Seurat_obj@meta.data$VIM1 + Seurat_obj@meta.data$EPCAM1+ Seurat_obj@meta.data$CDH1 + Seurat_obj@meta.data$ALDH1

#EMT score VIM - EPCAM - CDH1 + CD44 - CD24
Seurat_obj@meta.data$EMT_like_CSC <- Seurat_obj@meta.data$VIM1 + -1*Seurat_obj@meta.data$EPCAM1+ -1*Seurat_obj@meta.data$CDH1 + Seurat_obj@meta.data$CD441 + -1*Seurat_obj@meta.data$CD241 

#add broad score in
Seurat_obj <- AddModuleScore(Seurat_obj, features = list(c("ALDH1A1","ALDH1A2","ALDH1A3","NANOG","PROM1","POU5F1","SOX2","STAT3","ATP6AP2","MSI1","CD44")),
                             name = "CSC_broad")

samples <- Seurat_obj@meta.data

patient_col_dict <- c("sample_1" = "#66D7D1", "sample_2" = "#BEE3DB","sample_3" = "#D3D3D3",
                      "sample_4" = "#FFD6BA","sample_5" = "#274C77","sample_6" = "#A3CEF1", 
                      "sample_7" = "#F7CE5B", 
                      "sample_8" = "#7CEA9C","sample_9" = "#7C90A0",
                      "sample_10" = "#DCCCFF","sample_11" = "#EE6055",
                      "sample_12" = "#C492B1")

samples$sample_id <- factor(samples$sample_id, levels = 
                           c(names(patient_col_dict)))

samples$cluster_annotations <- factor(samples$cluster_annotations, levels = 
                                   c("core","transitory","edge"))

comparision_group <- "cluster_annotations"
comparisions <- list(c("core","edge"))
sample <- "sample_id"

#CSC score plots for Fig 3
for (values in c("MET_like_CSC","EMT_like_CSC","CSC_broad1")) {
  
  png(paste(values,"stat_plot.png",sep = "_"), units="in", width=5, height=4, res=300, type = "cairo")
  print(plot_effects(samples,comparision_group, comparisions,values,sample,patient_col_dict))
  dev.off()
  
}

#nebulosa plots for Fig 3

png("CSC_broad_nebulosa.png", units="in", width=5.5, height=6, res=300, type = "cairo")
print(SCpubr::do_NebulosaPlot(Seurat_obj, features = "CSC_broad1", method = "wkde"))
dev.off()

png("ALDH1A1and2and3_feature.png", units="in", width=6, height=6, res=300, type = "cairo")
SCpubr::do_FeaturePlot(Seurat_obj, features = c("ALDH11"), order = TRUE)
dev.off()

png("ALDH1A1and2and3_density.png", units="in", width=6, height=6, res=300, type = "cairo")
SCpubr::do_NebulosaPlot(Seurat_obj, features = c("ALDH11"))
dev.off()

png("CD44_feature.png", units="in", width=6, height=6, res=300, type = "cairo")
SCpubr::do_FeaturePlot(Seurat_obj, features = "CD44", order = TRUE)
dev.off()

png("CD44_density.png", units="in", width=6, height=6, res=300, type = "cairo")
SCpubr::do_NebulosaPlot(Seurat_obj, features = "CD44")
dev.off()

png("CD24_feature.png", units="in", width=6, height=6, res=300, type = "cairo")
SCpubr::do_FeaturePlot(Seurat_obj, features = "CD24", order = TRUE)
dev.off()

png("CD24_density.png", units="in", width=6, height=6, res=300, type = "cairo")
SCpubr::do_NebulosaPlot(Seurat_obj, features = "CD24")
dev.off()

png("CD24_inverse_density.png", units="in", width=6, height=6, res=300, type = "cairo")
SCpubr::do_NebulosaPlot(Seurat_obj, features = "CD24", viridis_direction = -1)
dev.off()
