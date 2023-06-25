library(Matrix)
library(plyr)
library(dplyr)
library(Seurat)
library(igraph)
library(ComplexHeatmap)
library(circlize)
require(dplyr)
require(ggplot2)
library(ggpubr)
require(cowplot)
library(data.table)
library(RColorBrewer)
library(scater)
library(tidyverse)
library(AUCell)
library(SCENIC)
library(tibble)

custom_fill_colors <- c("sample_1" = "#66D7D1", "sample_2" = "#BEE3DB","sample_3" = "#D3D3D3",
                                            "sample_4" = "#FFD6BA","sample_5" = "#274C77","sample_6" = "#A3CEF1", 
                                            "sample_7" = "#F7CE5B", 
                                            "sample_8" = "#7CEA9C","sample_9" = "#7C90A0",
                                            "sample_10" = "#DCCCFF","sample_11" = "#EE6055",
                                            "sample_12" = "#C492B1")

#load file where TF AUC has been added as an assay to a Seurat object of "comb"

SCENIC <- readRDS("x_AUC.rds")

Idents(SCENIC) = "cluster_annotations"
DefaultAssay(SCENIC) <- "AUC"
SCENIC_proj <- subset(SCENIC, idents = c("core","edge"))

sub <- SCENIC_proj
sample_ids.2 <- list(unique(sub$sample_id))

markers_list <- list()

for (k in 1:length(sample_ids.2[[1]])) {
  
      markers <- subset(sub, cells = c(which(sub$sample_id == sample_ids.2[[1]][k]))) %>%
        FindMarkers(ident.1 ="core", group.by = "cluster_annotations",logfc.threshold = 0,min.pct = 0) %>%
        rownames_to_column(var = "gene")%>%
        add_column(sample_id = sample_ids.2[[1]][k])
      
      markers_list[[k]] <- markers
  
  
}

markers.list <- plyr::ldply(markers_list, data.frame)
#pval adj 0.05
markers.list.mat <- markers.list[markers.list$p_val_adj<0.05,]
markers.list.mat <- markers.list.mat[,c("gene", "avg_log2FC", "sample_id")]

markers.list.mat <- reshape2::dcast(markers.list.mat, formula = gene~sample_id, value.var = "avg_log2FC")

markers.list.mat[is.na(markers.list.mat)] = 0
markers.list.mat <- markers.list.mat %>% tibble::column_to_rownames(var = "gene")


markers.matrix = markers.list.mat

markers.matrix.m <- markers.matrix %>% tibble::rownames_to_column(var = "gene")
markers.matrix.m <- reshape2::melt(markers.matrix.m)
order.markers <- rownames(markers.matrix[order(rowSums(-markers.matrix)),])
markers.sum <- as.data.frame(ifelse(markers.matrix != 0, 1, 0))
markers.sum$total <- rowSums(markers.sum)
#at least 9 samples
keep <- rownames(markers.sum[markers.sum$total>9,])
markers.matrix.m <- markers.matrix.m[markers.matrix.m$gene %in% keep,]

#genes is only used as this code is repurposed from the degs script, these are not genes
#these are TFs

core_genes <- markers.matrix.m[markers.matrix.m$value > 0,]
edge_genes <- markers.matrix.m[markers.matrix.m$value < 0,]
core_genes <- unique(core_genes$gene)
edge_genes <- unique(edge_genes$gene)

saveRDS(core_genes,"TF_lists//core_TFS.RDS")
saveRDS(edge_genes,"TF_lists/edge_TFS.RDS")

markers.matrix.m <- markers.matrix.m[markers.matrix.m$gene %in% union(edge_genes,core_genes),]

markers_table <- reshape2::dcast(markers.matrix.m, formula = gene~variable, value.var = "value")
markers_table <- column_to_rownames(markers_table, "gene")
markers_table$sum <- rowSums(markers_table)
write.csv(markers_table, file = "markers_table.csv")


png("diff_exp.png", units="in", width=15, height=4, res=300, type = "cairo")
print(ggplot(markers.matrix.m, aes(x = factor(gene, level = order.markers),
                                 y = value, fill = variable)) +
        geom_bar(stat="identity", color = "black", size = 0.25) + theme_minimal() +
        labs(x = "", y = "Cumulative average log(fold-change)", fill = "sample_id") +
        ggpubr::rotate_x_text() + scale_fill_manual(values = custom_fill_colors)+ 
        theme(axis.text = element_text(face="bold"), legend.text = element_text(face="bold"), 
              axis.text.y = element_text(face="bold")))
dev.off()

#only plot top 30 of each for supplementary fig
topweights <-aggregate(value~gene,markers.matrix.m,sum)

topweights <- topweights[order(-topweights$value),] 

markers_top_25 <- markers.matrix.m[markers.matrix.m$gene %in% c(head(topweights,25)$gene,tail(topweights,25)$gene),]

png("diff_exp_top.png", units="in", width=8, height=4, res=300, type = "cairo")
print(ggplot(markers_top_25, aes(x = factor(gene, level = order.markers),
                                   y = value, fill = variable)) +
        geom_bar(stat="identity", color = "black", size = 0.25) + theme_minimal() +
        labs(x = "", y = "Cumulative average log(fold-change)", fill = "sample_id") +
        ggpubr::rotate_x_text() + scale_fill_manual(values = custom_fill_colors)+ 
        theme(axis.text = element_text(face="bold"), legend.text = element_text(face="bold"), 
              axis.text.y = element_text(face="bold")))
dev.off()

setwd("..")

#reformat for interpretability
markers.list <- plyr::ldply(markers_list, data.frame)

markers.list.mat <- markers.list[markers.list$p_val_adj<0.05,]
markers.list.mat <- markers.list.mat[,c("gene", "avg_log2FC","p_val_adj","sample_id")]

markers.list.mat$enrichment <- ifelse(markers.list.mat$gene %in% core_genes, "core",
                                      "edge")

write.csv(markers.list.mat, file = "markers.list.mat.csv")
