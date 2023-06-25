library(Seurat)
library(dplyr)
library(tibble)
library(cowplot)
library(patchwork)
load(file = "annotated_objects.Robj")
load(file = "comb.Robj")
 
comb_meta <- comb@meta.data
#transfer annotations of core/transitory/edge back to objects with spatial positions
for (object in annotated_objects) {
  meta <- object@meta.data
  
  meta_sample_id <- left_join(meta,comb_meta, by = "Row.names")%>%column_to_rownames("Row.names")
  
  meta_sample_id$cluster_annotations <- factor(meta_sample_id$cluster_annotations, levels = c("core","edge","transitory","nc"))
  meta_sample_id$cluster_annotations[is.na(meta_sample_id$cluster_annotations)] <- "nc"
  
  annotated_objects[[unique(object$sample_id)]]@meta.data <- meta_sample_id
}

#barplot of state for core transtiroy adn edge regions 
png("core_transitory_edge/barplot.png",units = "in", res = 300, width = 5, height = 5, type = "cairo")
ggplot(comb@meta.data, aes(x = factor(comb@meta.data$sample_id, levels = c("sample_1","sample_2","sample_3","sample_4","sample_5","sample_6","sample_7","sample_8",
                                                       "sample_9","sample_10", "sample_11","sample_12")), fill = cluster_annotations))  + 
  geom_bar(position="fill")+
  theme_bw() + 
  labs(x=c("spatial annotation"), y="Proportion")+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  guides(fill = guide_legend(title = 'cancer cell state')) +
  scale_fill_manual(values = c('edge' = '#E64B35FF',
                               'transitory' = '#F9E076',
                               'core' = '#4DBBD5FF'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()

#core and edge annoations for each sample (used in Figure 2)
for(object in annotated_objects) {
  png(paste0("core_transitory_edge/",unique(object@meta.data$sample_id.x),"core_edge_anno.png"),units = "in", res = 300, width = 5, height = 5, type = "cairo")
  print(SpatialDimPlot(object, group.by = "cluster_annotations", images = "tumor",image.alpha = 0.5,
                       stroke = 0,
                       cols = c('edge' = '#E64B35FF',
                                'transitory' = '#F9E076',
                                'core' = '#4DBBD5FF')) + NoLegend())
  dev.off()
  
}

# lapply(
#   annotated_objects,
#   function(sample) {
#     SpatialDimPlot(sample,group.by  = "cluster_annotations")
#   }
# ) %>% 
#   wrap_plots(guides = 'collect')

all_samples <- annotated_objects

save(all_samples, file = "all_samples.Robj")

load(file = "all_samples.Robj")

#### DEG analysis

#define colours
patient_col_dict <- c("sample_1" = "#66D7D1", "sample_2" = "#BEE3DB","sample_3" = "#D3D3D3",
                      "sample_4" = "#FFD6BA","sample_5" = "#274C77","sample_6" = "#A3CEF1", 
                      "sample_7" = "#F7CE5B", 
                      "sample_8" = "#7CEA9C","sample_9" = "#7C90A0",
                      "sample_10" = "#DCCCFF","sample_11" = "#EE6055",
                      "sample_12" = "#C492B1")
markers_list = c()
i = 1
for (Seurat_obj in all_samples) { 

  Idents(Seurat_obj) <- "cluster_annotations"
  #scale and normalize each sample individually
  Seurat_obj <- ScaleData(Seurat_obj)
  Seurat_obj <- NormalizeData(Seurat_obj)
  sub_obj <- subset(Seurat_obj, idents = c("core","edge"))
  
  #core is positive, edge is negative
  results <- FindMarkers(sub_obj,ident.1  = "core",ident.2 = "edge", group.by = "cluster_annotations")%>%
    tibble::rownames_to_column(var = "gene")%>%
    tibble::add_column("sample_id.x" = unique(sub_obj@meta.data$sample_id.x))
  
  write.csv(results, file = paste("differential_expression/raw_csvs/",unique(sub_obj@meta.data$sample_id.x),"_raw_deg.csv", sep = ""))
  markers_list[[i]] <- results
  i = i+1
}

#code adapted to make consensus plots
library(dplyr)

markers_list_saved <- markers_list

markers.list <- plyr::ldply(markers_list, data.frame)

#pval cutoff of 0.001

markers.list.mat <- markers.list[markers.list$p_val_adj<0.001,]

markers.list.mat <- markers.list.mat[,c("gene", "avg_log2FC", "sample_id.x")]

markers.list.mat <- reshape2::dcast(markers.list.mat, formula = gene~sample_id.x, value.var = "avg_log2FC")

# markers.list.mat <- reshape2::dcast(markers.list.mat, formula = gene~sample_id.x, value.var = "log2FC")

markers.list.mat[is.na(markers.list.mat)] = 0
markers.list.mat <- markers.list.mat %>% tibble::column_to_rownames(var = "gene")

markers.matrix = markers.list.mat

markers.matrix.m <- markers.matrix %>% tibble::rownames_to_column(var = "gene")
markers.matrix.m <- reshape2::melt(markers.matrix.m)
order.markers <- rownames(markers.matrix[order(rowSums(-markers.matrix)),])
markers.sum <- as.data.frame(ifelse(markers.matrix != 0, 1, 0))
markers.sum$total <- rowSums(markers.sum)
#must be differnetially expressed in at least 9 samples
keep <- rownames(markers.sum[markers.sum$total>9,])
markers.matrix.m <- markers.matrix.m[markers.matrix.m$gene %in% keep,]

markers_table <- reshape2::dcast(markers.matrix.m, formula = gene~variable, value.var = "value")
markers_table <- column_to_rownames(markers_table, "gene")
markers_table$sum <- rowSums(markers_table)
write.csv(markers_table, file = "markers_table.csv")

core_genes <- markers.matrix.m[markers.matrix.m$value > 0,]
edge_genes <- markers.matrix.m[markers.matrix.m$value < 0,]

core_genes <- unique(core_genes$gene)
edge_genes <- unique(edge_genes$gene)

saveRDS(core_genes,"gene_lists/core_genes.RDS")
saveRDS(edge_genes,"gene_lists/edge_genes.RDS")

markers.matrix.m <- markers.matrix.m[markers.matrix.m$gene %in% union(edge_genes,core_genes),]

setwd("differential_expression/")
#create consensus plot 
png("differentially_expressed_genes_sct_seurat.png", units="in", width=20, height=7, res=300, type = "cairo")
ggplot(markers.matrix.m, aes(x = factor(gene, level = order.markers), 
                             y = value, fill = variable)) +
  geom_bar(stat="identity", color = "black", size = 0.25) + theme_minimal() + 
  scale_fill_manual( values = patient_col_dict)+
  labs(x = "", y = "Cumulative average log(fold-change)", fill = "sample_id.x") +
  ggpubr::rotate_x_text()
dev.off()

#only plot top 25 genes for Figure 2
topweights <-aggregate(value~gene,markers.matrix.m,sum)

topweights <- topweights[order(-topweights$value),] 

markers_top_25 <- markers.matrix.m[markers.matrix.m$gene %in% c(head(topweights,25)$gene,tail(topweights,25)$gene),]

png("deg_top25.png", units="in", width=9, height=3, res=300, type = "cairo")
ggplot(markers_top_25, aes(x = factor(gene, level = order.markers), 
                             y = value, fill = variable)) +
  geom_bar(stat="identity", color = "black", size = 0.25) + theme_minimal() + 
  scale_fill_manual( values = patient_col_dict)+
  labs(x = "", y = "Cumulative average log(fold-change)", fill = "sample_id.x") +
  ggpubr::rotate_x_text()
dev.off()


setwd("..")

#do similar analysis to extract raw data for each sample for IPA

markers_list = c()
i = 1
for (Seurat_obj in all_samples) { 
 
  Idents(Seurat_obj) <- "cluster_annotations"
  Seurat_obj <- ScaleData(Seurat_obj)
  Seurat_obj <- NormalizeData(Seurat_obj)
  sub_obj <- subset(Seurat_obj, idents = c("core","edge"))
  
results <- FindMarkers(sub_obj,ident.1  = "core",ident.2 = "edge", group.by = "cluster_annotations")%>%
  tibble::rownames_to_column(var = "gene")%>%
  tibble::add_column("sample_id.x" = unique(sub_obj@meta.data$sample_id.x))
 
 results <- results[results$p_val_adj < 0.001,]
 results_edge <- results[results$avg_log2FC < 0,]
 results_edge$avg_log2FC <- results_edge$avg_log2FC*-1
 results_core <- results[results$avg_log2FC > 0,]
 write.csv(results_edge, file = paste("differential_expression/IPA/",unique(sub_obj@meta.data$sample_id.x),"_edge_raw_deg.csv", sep = ""))
 write.csv(results_core, file = paste("differential_expression/IPA/",unique(sub_obj@meta.data$sample_id.x),"_core_raw_deg.csv", sep = ""))
 
 markers_list[[i]] <- results
 i = i+1
}


#check overlap with puram et al. genesets

puram_core <- scan("puram_genesets/epithelial_diff.txt", what="", sep="\n")

puram_edge <- scan("puram_genesets/p_emt.txt", what="", sep="\n")

core_genes <- readRDS("gene_lists/core_genes.RDS")

edge_genes <- readRDS("gene_lists/edge_genes.RDS")

length(intersect(puram_core,core_genes))

#47/80 overlapping for epithelial diff

length(intersect(puram_edge,edge_genes))

#only 7/100 edge overlap
           