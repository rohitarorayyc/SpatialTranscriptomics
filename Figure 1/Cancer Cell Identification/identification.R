library(tidyverse)
library(dplyr)
library(Seurat)
library(patchwork)
load(file = "Seurat_objs.Robj")

load(file = "nb.Robj")

samples = c('sample_1','sample_2','sample_3','sample_4','sample_5','sample_6','sample_7','sample_8','sample_9','sample_10','sample_11','sample_12')
for (sample in samples) {
  
  Seurat_meta <- Seurat_objs[[sample]]@meta.data
  Seurat_meta$barcode <- rownames(Seurat_meta)
  barcode_id <- strsplit(rownames(Seurat_meta), "-")[[1]][[2]]
  p_cnvs <- nb[[sample]]$clone_post
  
  p_cnvs$cell <- str_replace(p_cnvs$cell, "1", barcode_id)
  
  Seurat_objs[[sample]]@meta.data <- Seurat_meta %>%
    left_join(
      p_cnvs, 
      by = c('barcode' = 'cell')
    ) %>%
    column_to_rownames("barcode")
  
  #define cancer cell spots as being either p cnv > 0.99 or deconvolution probability > 0.99
  Seurat_objs[[sample]]$is_cancer_cell <- ifelse( (Seurat_objs[[sample]]$p_cnv > 0.99 |  Seurat_objs[[sample]]$cancer.cell > 0.99) & Seurat_objs[[sample]]$pathologist_anno == "SCC", "1", "0")
  
  #define non cancer cell states
  Seurat_objs[[sample]]$celltype <- ifelse(Seurat_objs[[sample]]$is_cancer_cell == "1","cancer cell","other celltype")
  
  meta <- Seurat_objs[[sample]]@meta.data

  noncancer_meta <- meta[c("myofibroblast","B.cell","ecm.myCAF","Intermediate.fibroblast","detox.iCAF","macrophage","endothelial","dendritic.","mast",
         "conventional.CD4..T.helper.cells","cytotoxic.CD8..T.","Tregs","cytotoxic.CD8..T.exhausted")]
  
  #identify the noncancer celltype as the noncancer celltype with the greatest proportion
  noncancer_meta <- noncancer_meta %>% mutate(noncancer_celltype=names(.)[max.col(.)])
  noncancer_meta$barcode <- rownames(noncancer_meta)
  
  Seurat_objs[[sample]]@meta.data <- merge(x = Seurat_objs[[sample]]@meta.data, y = noncancer_meta[ , c("noncancer_celltype","barcode")], all.x=TRUE, by = 0) %>%
    column_to_rownames("barcode")
  
  Seurat_objs[[sample]]@meta.data$final_celltype <- ifelse(Seurat_objs[[sample]]$is_cancer_cell == "1","cancer cell",Seurat_objs[[sample]]@meta.data$noncancer_celltype)
}

#check all annotations
lapply(
  samples,
  function(sample) {
    SpatialDimPlot(Seurat_objs[[sample]], group.by = "final_celltype", interactive = T)
  }
) %>% 
  wrap_plots(guides = 'collect')

annotated_objects <- Seurat_objs

save(annotated_objects, file = "annotated_objects.Robj")

#plots for Figure 1 were generated in Figure 2, Core Edge Identification code
