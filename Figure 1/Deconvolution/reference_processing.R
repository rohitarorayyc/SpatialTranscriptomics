###load in and process single-cell data from GEO
library(Seurat)

all_data <- read.delim2("HNSCC_all_data.txt")

all_data$X  <- as.list(sapply(all_data$X , function(x) gsub("\'", "", x)))

rownames(all_data) <- all_data$X

all_data$X <- NULL

meta_data <- all_data[0:5,]
meta_data <- as.data.frame(t(meta_data))

expr_data <- all_data[6:23691,]

seurat_data <- CreateSeuratObject(expr_data, project = "puram_data")

seurat_data <- AddMetaData(seurat_data,metadata = meta_data)

library(reshape2)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(tidyr)
#heatmap

seurat_data <- FindVariableFeatures(seurat_data)
seurat_data <- ScaleData(seurat_data)
seurat_data <- Seurat::RunPCA(seurat_data, verbose = FALSE) %>%
  Seurat::RunUMAP(., dims = 1:30, verbose = FALSE)%>%
  FindNeighbors(., dims = 1:30)%>%
  FindClusters(., resolution = 0.5)

seurat_data@meta.data$non.cancer.cell.type <- gsub("-","",as.character(seurat_data@meta.data$non.cancer.cell.type))

#rename cancer cells
seurat_data@meta.data$non.cancer.cell.type <- gsub("0","Cancer cell",as.character(seurat_data@meta.data$non.cancer.cell.type))

Seurat::DimPlot(seurat_data,
                group.by = "RNA_snn_res.0.5",
                label = TRUE) + Seurat::NoLegend()

#identify celltypes based on published markers
cluster_ids <- read.csv(file = "cluster_identification.csv")

Idents(seurat_data) = "RNA_snn_res.0.5"
new.cluster.ids <- cluster_ids$cell_type
names(new.cluster.ids) <- levels(seurat_data)
seurat_data <- RenameIdents(seurat_data, new.cluster.ids)
seurat_data@meta.data$cellype_fine <- Idents(seurat_data)

colourCount = length(unique(seurat_data@meta.data$cellype_fine))
library(RColorBrewer)
mycolors <- colorRampPalette(brewer.pal(7, "Set1"))(colourCount)

#refactor

seurat_data@meta.data$cellype_fine <- factor(seurat_data@meta.data$cellype_fine, levels = 
                                               c("myofibroblast","cancer cell","B cell","ecm-myCAF",
                                                 "Intermediate fibroblast","detox-iCAF","macrophage",
                                                 "endothelial","dendritic ","mast","conventional CD4+ T-helper cells",
                                                 "cytotoxic CD8+ T ","Tregs","cytotoxic CD8+ T exhausted"))

Idents(seurat_data) = "cellype_fine"
png("puram_corescore_subtype.png",units = "in", res = 300, width = 12, height = 8)
VlnPlot(seurat_data, features = "core_genes1", cols = mycolors)
dev.off()

png("puram_edgescore_subtype.png",units = "in", res = 300, width = 12, height = 8)
VlnPlot(seurat_data, features = "edge_genes1", cols = mycolors)
dev.off()


png("puram_annotated.png",units = "in", res = 300, width = 8, height = 6)
Seurat::DimPlot(seurat_data,
                group.by = "cellype_fine",
                label = F, cols = mycolors) 
dev.off()
puram_data <- seurat_data
save(puram_data, file = "puram_data.Robj")                    

