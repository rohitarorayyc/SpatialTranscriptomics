library(reshape2)
library(patchwork)
library(ggplot2)
library(dplyr)
library(Seurat)
library(RColorBrewer)
library(sctransform)

load(file = "comb.Robj")

Idents(comb) = "cluster_annotations"
comb_core_edge <- subset(comb, idents = c("core","edge"))

comb_core_edge@meta.data$sample_group <- paste(comb_core_edge@meta.data$sample_id, 
                                                      comb_core_edge@meta.data$cluster_annotations, sep = "_")

object <- comb_core_edge
Idents(object) = "sample_group"
object@meta.data$sample_group <- factor(object@meta.data$sample_group, levels = c(
  "sample_1_edge",
  "sample_2_edge",
  "sample_3_edge",
  "sample_4_edge",
  "sample_5_edge",
  "sample_6_edge",
  "sample_7_edge",
  "sample_8_edge",
  "sample_9_edge",
  "sample_10_edge",
  "sample_11_edge",
  "sample_12_edge",
  "sample_1_core",
  "sample_2_core",
  "sample_3_core",
  "sample_4_core",
  "sample_5_core",
  "sample_6_core",
  "sample_7_core",
  "sample_8_core",
  "sample_9_core",
  "sample_10_core",
  "sample_11_core",
  "sample_12_core"))
Idents(object) = "sample_group"

av.exp <- AverageExpression(object)$SCT
av.exp <- na.omit(av.exp)
cor.exp <- as.data.frame(cor(av.exp))
cor.exp$x <- rownames(cor.exp)
cor.df <- tidyr::gather(data = cor.exp, y, correlation, unique(comb_core_edge@meta.data$sample_group))

setwd("corplot/")

group = "cluster_annotations"
cor.exp$x = NULL
cor.exp <- as.matrix(cor.exp)
cols <- RColorBrewer::brewer.pal(11, "Spectral")
cols <- rev(cols)

group1 <- unique(object@meta.data[[group]])[1]
group2 <- unique(object@meta.data[[group]])[2]


cor.exp2 = data.frame(group = c(rep(group1, length(unique(object@meta.data$sample_group))/2), 
                                rep(group2, length(unique(object@meta.data$sample_group))/2)))


list1 = list(group = c(group1 = '#E64B35FF',
                       group2 = '#4DBBD5FF'))

list2 = list(setNames(list1$group, c(group1, group2)))

names(list2) = "group"
library(ComplexHeatmap)
ha = HeatmapAnnotation(df = cor.exp2, col = list2)
hab = rowAnnotation(df = cor.exp2, col = list2)

#generate correlation plot for Figure 2
png("heatmap_corr.png",units = "in", res = 300, width = 10, height = 9,type = "cairo")

print(ComplexHeatmap::Heatmap(cor.exp, col = cols, bottom_annotation = ha, 
                              right_annotation = hab, column_split = c(rep(group1, length(unique(object@meta.data$sample_group))/2), 
                                                                       rep(group2, length(unique(object@meta.data$sample_group))/2)), 
                              
                              row_split = c(rep(group1, length(unique(object@meta.data$sample_group))/2), 
                                            rep(group2, length(unique(object@meta.data$sample_group))/2)),
                              cluster_rows = T,
                              cluster_columns = T
))
dev.off()
