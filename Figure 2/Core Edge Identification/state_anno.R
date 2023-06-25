library(Seurat)
library(harmony)
library(cowplot)
library(RColorBrewer)
library(ggplot2)
load(file = "annotated_objects.Robj")

all_spots <- purrr::reduce(annotated_objects, merge)

DefaultAssay(all_spots) <- "Spatial"

#integration by sample

obj_list <- SplitObject(all_spots, split.by = "sample_id")
rm(list=setdiff(ls(), "obj_list"))

obj_list <- lapply(X = obj_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

anchors <- FindIntegrationAnchors(object.list = obj_list, dims = 1:30)

comb_as <- IntegrateData(anchorset = anchors, dims = 1:30)
comb_as <- ScaleData(comb_as, verbose = FALSE)
comb_as <- RunPCA(comb_as, npcs = 50, verbose = FALSE)

png("elbowplot_all.png",units = "in", res = 300, width = 6, height = 4, type = "cairo")
ElbowPlot(comb_as, ndims = 25)
dev.off()

comb_as <- RunUMAP(comb_as, reduction = "pca", dims = 1:15)
comb_as <- FindNeighbors(comb_as, reduction = "pca", dims = 1:15)
#comb_as <- FindClusters(comb_as, resolution = 0.5)

DefaultAssay(comb_as) <- "integrated"

save(comb_as, file = "comb_as.Robj")

comb_as <- RunPCA(comb_as, npcs = 6, verbose = TRUE, seed.use = 12)
comb_as <- RunUMAP(comb_as, reduction = "pca", dims = 1:6, seed.use = 12)
comb_as <- FindNeighbors(comb_as, reduction = "pca", dims = 1:6, seed.use = 12)

#plot for pathologist annotations used in Figure 1
png("all_samples_patho.png",units = "in", res = 300, width = 6, height = 4, type = "cairo")
DimPlot(comb_as, reduction = "umap", group.by = "pathologist_anno", 
        cols = c("SCC" = "maroon", "Lymphocyte Negative Stroma" = "darkblue", 
        "Lymphocyte Positive Stroma" = "lightblue", "Keratin" = "lightgreen", "Fold" = "grey", 
        "Artifact" = "grey", "Artery/Vein" = "orange",
        "Muscle" = "lightyellow","Non-cancerous Mucosa" = "purple",
        "Glandular Stroma" = "violet"))
dev.off()

#plot for cancer cell identification by pcnv used in Figure 1
png("all_samples_pcnv.png",units = "in", res = 300, width = 5, height = 5, type = "cairo")
FeaturePlot(comb_as,features = "p_cnv", min.cutoff = 0.99, pt.size = 0.2,cols = c("grey", "#D9381E"),
            max.cutoff = 1)
dev.off()

#plot for cancer cell identification by deconvolution used in Figure 1
png("all_samples_deconv.png",units = "in", res = 300, width = 5, height = 5, type = "cairo")
FeaturePlot(comb_as,features = "cancer.cell", min.cutoff = 0.99, pt.size = 0.2,cols = c("grey", "#DC143C"),
            max.cutoff = 1)
dev.off()

colourCount = length(unique(comb_as$final_celltype))

#plot for celltype used in Figure 1.
png("all_samples_final_celltype.png",units = "in", res = 300, width = 6.5, height = 5, type = "cairo")
DimPlot(comb_as, reduction = "umap", group.by = "final_celltype", 
        cols = c("B.cell" = cols[1], "cancer cell"= cols[2], "cytotoxic.CD8..T." = cols[3],
        "cytotoxic.CD8..T.exhausted" = cols[4],
        "dendritic." = cols[5],
        "detox.iCAF" = cols[6],
        "ecm.myCAF" = cols[7],
        "endothelial" = cols[8],
        "Intermediate.fibroblast" = cols[9],
        "macrophage" = cols[10],
        "mast" = cols[11],
        "myofibroblast" = cols[12],
        "Tregs" = cols[13],
        "NA" = cols[14]))
dev.off()

#plot for cancer cell identification used in figure 1
png("all_samples_iscancercell.png",units = "in", res = 300, width = 5, height = 5, type = "cairo")
DimPlot(comb_as, reduction = "umap", group.by = "is_cancer_cell",cols = c("grey", "red"))
dev.off()

#example plot for sample 1 of pathologist annotations used in Figure 1
png("sample_1_patho_anno.png",units = "in", res = 300, width = 5, height = 5, type = "cairo")
SpatialDimPlot(obj_list[[1]], group.by = "pathologist_anno", images = "tumor", image.alpha = 0.5,
               stroke = 0,
        cols = c("SCC" = "maroon", "Lymphocyte Negative Stroma" = "darkblue", 
                 "Lymphocyte Positive Stroma" = "lightblue", "Keratin" = "lightgreen", "Fold" = "grey", 
                 "Artifact" = "grey", "Artery/Vein" = "orange",
                 "Muscle" = "lightyellow","Non-cancerous Mucosa" = "purple",
                 "Glandular Stroma" = "violet")) + NoLegend()
dev.off()

#example plot for pcnv for sample 1 used in Figure 1
png("sample_1_pcnv.png",units = "in", res = 300, width = 5, height = 5, type = "cairo")
SpatialFeaturePlot(obj_list[[1]],features = "p_cnv", images = "tumor", min.cutoff = 0.99,stroke = 0,
                   image.alpha = 0.5) +
  scale_fill_gradient(low = "grey", high = "#D9381E")+ NoLegend()
dev.off()

#example plot for deconvolution for sample 1 used in Figure 1
png("sample_1_deconv.png",units = "in", res = 300, width = 5, height = 5, type = "cairo")
SpatialFeaturePlot(obj_list[[1]],features = "cancer.cell", images = "tumor", min.cutoff = 0.99,
                   stroke = 0,
                   image.alpha = 0.5) +
  scale_fill_gradient(low = "grey", high = "#DC143C")+ NoLegend()
dev.off()

#example plot for cancer cell status for sample 1 used in Figure 1
png("sample_1_iscancercell.png",units = "in", res = 300, width = 5, height = 5, type = "cairo")
SpatialDimPlot(obj_list[[1]], group.by = "is_cancer_cell", images = "tumor",image.alpha = 0.5,

               stroke = 0,
               cols = c("0" = "grey","1" = "red")) + NoLegend()
dev.off()

cols = colorRampPalette(brewer.pal(9, "Set1"))(colourCount)
#example plot for final celltype for sample 1 used in Figure 1
png("sample_1_final_celltype.png",units = "in", res = 300, width = 5, height = 5, type = "cairo")
SpatialDimPlot(annotated_objects[[1]], group.by = "final_celltype", images = "tumor",image.alpha = 0.5,
               stroke = 0,
              cols = c("B.cell" = cols[1], "cancer cell"= cols[2], "cytotoxic.CD8..T." = cols[3],
                       "cytotoxic.CD8..T.exhausted" = cols[4],
                       "detox.iCAF" = cols[6],
                       "ecm.myCAF" = cols[7],
                       "endothelial" = cols[8],
                       "Intermediate.fibroblast" = cols[9],
                       "macrophage" = cols[10],
                       "mast" = cols[11],
                       "myofibroblast" = cols[12],
                       "Tregs" = cols[13]))+ NoLegend()

#plot proportions of celltypes for supplementary fig 1
all_data <- all_spots@meta.data

all_data[is.na(all_data)] <- "unknown"
colourCount = length(unique(all_data$final_celltype))

library(RColorBrewer)
all_data$sample_id <- factor(all_data$sample_id, levels = c("sample_1","sample_2","sample_3","sample_4","sample_5","sample_6","sample_7","sample_8",
                                                            "sample_9","sample_10", "sample_11","sample_12"))
png("celltype_bar.png",units = "in", res = 300, width = 7, height = 5, type = "cairo")
ggplot(all_data, aes(x = sample_id, fill = final_celltype))  + 
  geom_bar(position="fill")+
  theme_bw() + 
  labs(x=c("spatial annotation"), y="Proportion")+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  guides(fill = guide_legend(title = 'cell type')) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(colourCount))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#plot proportions of pathologist regions for supplementary fig 1
all_data <- all_spots@meta.data

all_data[is.na(all_data)] <- "unknown"
colourCount = length(unique(all_data$pathologist_anno))

library(RColorBrewer)
all_data$sample_id <- factor(all_data$sample_id, levels = c("sample_1","sample_2","sample_3","sample_4","sample_5","sample_6","sample_7","sample_8",
                                                            "sample_9","sample_10", "sample_11","sample_12"))
png("celltype_bar_pathologist_anno.png",units = "in", res = 300, width = 7, height = 5, type = "cairo")
ggplot(all_data, aes(x = sample_id, fill = pathologist_anno))  + 
  geom_bar(position="fill")+
  theme_bw() + 
  labs(x=c("spatial annotation"), y="Proportion")+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  guides(fill = guide_legend(title = 'pathologist annotations')) +
  scale_fill_manual(values = c("SCC" = "maroon", "Lymphocyte Negative Stroma" = "darkblue", 
                               "Lymphocyte Positive Stroma" = "lightblue", "Keratin" = "lightgreen", "Fold" = "grey", 
                               "Artifact" = "grey", "Artery/Vein" = "orange",
                               "Muscle" = "lightyellow","Non-cancerous Mucosa" = "purple","Glandular Stroma" = "violet"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#plot dimplot of pathologist region for each sample and cancer cell annotation for each sample for supplementary fig 1
for(object in annotated_objects) {
  #image = paste0("tumor.",as.numeric(stringr::str_split(unique(object@meta.data$sample_id),"_")[[1]][2])-1)
  png(paste0("pathologist_anno/",unique(object@meta.data$sample_id),"patho_anno.png"),units = "in", res = 300, width = 5, height = 5, type = "cairo")
  print(SpatialDimPlot(object, group.by = "pathologist_anno", images = "tumor", image.alpha = 0.5,
                       stroke = 0,
                       cols = c("SCC" = "maroon", "Lymphocyte Negative Stroma" = "darkblue", 
                                "Lymphocyte Positive Stroma" = "lightblue", "Keratin" = "lightgreen", "Fold" = "grey", 
                                "Artifact" = "grey", "Artery/Vein" = "orange",
                                "Muscle" = "lightyellow","Non-cancerous Mucosa" = "purple",
                                "Glandular Stroma" = "violet")) + NoLegend())
  dev.off()
  
  png(paste0("cancer_cell_anno/",unique(object@meta.data$sample_id),"iscancercell.png"),units = "in", res = 300, width = 5, height = 5, type = "cairo")
  print(SpatialDimPlot(object, group.by = "is_cancer_cell", images = "tumor",image.alpha = 0.5,
                       stroke = 0,
                       cols = c("0" = "grey","1" = "red")) + NoLegend())
  dev.off()
  
}

all_data <- comb_as@meta.data

#Only unknown spots are in sample 3 and this is because there are a few low quality spots that are not on the tissue section
#Due to this, the "unknown status" will be assumed to be noncancer 
View(all_data$is_cancer_cell[is.na(all_data$is_cancer_cell)])
all_data[is.na(all_data)] <- "0"
all_data$sample_id <- factor(all_data$sample_id, levels = c("sample_1","sample_2","sample_3","sample_4","sample_5","sample_6","sample_7","sample_8",
                                                            "sample_9","sample_10", "sample_11","sample_12"))

png("malignant_barplot.png",units = "in", res = 300, width = 7, height = 5, type = "cairo")
ggplot(all_data, aes(x = sample_id, fill = is_cancer_cell))  + 
  geom_bar(position="fill")+
  theme_bw() + 
  labs(x=c("spatial annotation"), y="Proportion")+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  guides(fill = guide_legend(title = 'malignant status')) +
  scale_fill_manual(values = c("0" = "grey","1" = "red"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

comb_d2c <- comb
DefaultAssay(comb_d2c) <- "Spatial"

comb_d2c@assays$SCT = NULL

comb_d2c$cluster_annotations <- as.character(comb_d2c$cluster_annotations)
SaveH5Seurat(comb_d2c, filename = "comb_d2c.h5Seurat")
Convert("comb_d2c.h5Seurat", dest = "h5ad")

#integrate only cancer cells together
all_spots_c <- subset(all_spots, subset = final_celltype == "cancer cell")

DefaultAssay(all_spots_c) <- "Spatial"
#integration
obj_list <- SplitObject(all_spots_c, split.by = "sample_id")
rm(list=setdiff(ls(), "obj_list"))
obj_list <- lapply(X = obj_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

anchors <- FindIntegrationAnchors(object.list = obj_list, dims = 1:30)
comb <- IntegrateData(anchorset = anchors, dims = 1:30)
comb <- ScaleData(comb, verbose = FALSE)
comb <- RunPCA(comb, npcs = 50, verbose = FALSE)
comb <- RunUMAP(comb, reduction = "pca", dims = 1:50)

ElbowPlot(comb, ndims = 25)

comb <- RunUMAP(comb, reduction = "pca", dims = 1:10)
comb <- FindNeighbors(comb, reduction = "pca", dims = 1:10)

save(comb, file = "comb.Robj")

library(Seurat)

load(file = "comb.Robj")

comb <- Seurat::FindClusters(object = comb, resolution = 1)
# pull the tree
data.tree <- Tool(object = comb, slot = "BuildClusterTree")

#plot of clusters for fig 2
png("cluster_res1.png",units = "in", res = 300, width = 5, height = 4, type = "cairo")
Seurat::DimPlot(comb)
dev.off()

#plot of phylogenetic tree for fig 2
png("cluster_res1_phylo.png",units = "in", res = 300, width = 5, height = 5, type = "cairo")
ape::plot.phylo(x = data.tree, direction = "right")
dev.off()

comb <- Seurat::FindClusters(object = comb, resolution = 1)

trial_ids <- c("n3", "n2","n3","n2","n3", "n1","n3","n3","n2","n1","n3","n2","n2","n2","n1")

names(trial_ids) <- levels(comb)

comb <- RenameIdents(comb, trial_ids)

library(SCpubr)
library(ggplot2)

#plot of nodes for Fig 2
png("nodal_definitions.png",units = "in", res = 300, width = 5, height = 4, type = "cairo")
DimPlot(comb, reduction = "umap", label = FALSE, pt.size = 0.5, 
        cols = c("n1" = "#4DBBD5FF",
                 "n2" = "#F9E076", "n3" = "#E64B35FF")) 
dev.off()

load(file = "comb.Robj")
#heatmap of top markers
DefaultAssay(comb) <- "integrated"

comb <- NormalizeData(comb)
top_markers <- FindAllMarkers(comb, only.pos = TRUE)

#plot of top markers in each nodal group 
png("logfc_top_markers.png",units = "in", res = 300, width = 5, height = 5, type = "cairo")
SCpubr::do_GroupwiseDEPlot(sample = comb,
                           de_genes = tibble::tibble(top_markers),
                           top_genes = 5)
dev.off()

write.csv(top_markers, file = "node_markers.csv")

#nebulosa plot of validated markers for fig 2 
png("nebulosa_CLDN4.png",units = "in", res = 300, width = 4.5, height = 5, type = "cairo")
SCpubr::do_NebulosaPlot(comb, "CLDN4")
dev.off()

png("nebulosa_SPRR1B.png",units = "in", res = 300, width = 4.5, height = 5, type = "cairo")
SCpubr::do_NebulosaPlot(comb, "SPRR1B")
dev.off()

png("nebulosa_LAMC2.png",units = "in", res = 300, width = 4.5, height = 5, type = "cairo")
SCpubr::do_NebulosaPlot(comb, "LAMC2")
dev.off()

png("nebulosa_ITGA5.png",units = "in", res = 300, width = 4.5, height = 5, type = "cairo")
SCpubr::do_NebulosaPlot(comb, "ITGA5")
dev.off()

comb <- Seurat::FindClusters(object = comb, resolution = 1)

trial_ids <- c("edge", "transitory","edge","transitory","edge", "core","edge","edge","transitory","core","edge","transitory","transitory","transitory","core")

names(trial_ids) <- levels(comb)

comb <- RenameIdents(comb, trial_ids)

library(SCpubr)
library(ggplot2)

#renaming regions as core and edge
png("core_edge_transitory_anno.png",units = "in", res = 300, width = 5.5, height = 4, type = "cairo")
DimPlot(comb, reduction = "umap", label = FALSE, pt.size = 0.5, 
        cols = c("core" = "#4DBBD5FF",
                 "transitory" = "#F9E076", "edge" = "#E64B35FF")) 
dev.off()

##plot of clonotypes by state for supplementary fig 3
library(ggplot2)
library(RColorBrewer)
comb@meta.data$core_edge_anno
meta_df <- comb@meta.data
meta_df$clone_opt <- as.factor(meta_df$clone_opt)

png("cancer_cell_clone_byregion.png",units = "in", res = 300, width = 5, height = 5, type = "cairo")
ggplot(meta_df, aes(x = core_edge_anno, fill = clone_opt))  + 
  geom_bar(position="fill")+
  facet_wrap(~sample_id)+
  theme_bw() + 
  labs(x=c("sample ID"), y="Proportion of each clone in cancer cell regions")+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(6, "Set2"))(6))
dev.off()

# top_markers <- FindAllMarkers(comb, only.pos = TRUE)
# write.csv(top_markers, file = "node_markers.csv")



