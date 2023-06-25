library(CARD)
library(Seurat)
library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
#load all objects
load(file = "Seurat/obj_1.Robj")
load(file = "Seurat/obj_2.Robj")
load(file = "Seurat/obj_3.Robj")
load(file = "Seurat/obj_4.Robj")
load(file = "Seurat/obj_5.Robj")
load(file = "Seurat/obj_6.Robj")
load(file = "Seurat/obj_7.Robj")
load(file = "Seurat/obj_8.Robj")
load(file = "Seurat/obj_9.Robj")
load(file = "Seurat/obj_10.Robj")
load(file = "Seurat/obj_11.Robj")
load(file = "Seurat/obj_12.Robj")

Seurat_objs <- c(obj_1,obj_2,obj_3,obj_4,obj_5,obj_6,obj_7,obj_8,
                 obj_9,obj_10,obj_11,obj_12)

samples <- c("sample_1","sample_2","sample_3","sample_4","sample_5","sample_6","sample_7","sample_8",
             "sample_9","sample_10", "sample_11","sample_12")

names(Seurat_objs) <- samples



for (sample in Seurat_objs) {
sample@meta.data$sample <- NULL
spatial_count <- sample@assays$Spatial@counts

spatial_location <- sample@images$tumor@coordinates %>% select(row, col)
names(spatial_location) <- c("x","y")

if (!identical(colnames(spatial_count), rownames(spatial_location))){
#ensure that count matrix contains the same spots as the location matrix 
spatial_count <- spatial_count[, colnames(spatial_count) %in% rownames(spatial_location)]
}
#match order
spatial_location <- spatial_location[order(match(rownames(spatial_location), colnames(spatial_count))), , drop = FALSE]

load(file = "puram/puram_data.Robj")
sc_count <- puram_data@assays$RNA@counts
sc_meta <- puram_data@meta.data

#deconvolute using CARD
CARD_obj = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = "cellype_fine",
  ct.select = unique(sc_meta$cellype_fine),
  sample.varname = "orig.ident",
  minCountGene = 100,
  minCountSpot = 5) 

CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)

proportions <- as.data.frame(CARD_obj@Proportion_CARD)

proportions$dominant_celltype <- apply(proportions, 1, function(x) names(proportions)[which.max(x)])

sample <- Seurat::AddMetaData(sample,proportions)

Seurat_objs[unique(sample@meta.data$sample_id)] <- sample

}

# save(Seurat_objs,file = "Seurat_objs.Robj")
# 
# load(file = "Seurat_objs.Robj")

lapply(
  Seurat_objs,
  function(sample) {
    SpatialFeaturePlot(sample,features = "cancer.cell", min.cutoff = 1)
  }
) %>% 
  wrap_plots(guides = 'collect')

