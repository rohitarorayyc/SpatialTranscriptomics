library(reshape2)
library(patchwork)
library(ggplot2)
library(gprofiler2)
library(dplyr)
library(Seurat)
library(RColorBrewer)
library(sctransform)
library(CellChat)
library(patchwork)

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB

library(Seurat)
library(dplyr)
library(tibble)
library(cowplot)
library(patchwork)
load(file = "annotated_objects.Robj")
load(file = "comb.Robj")

comb_meta <- comb@meta.data

for (object in annotated_objects) {
  meta <- object@meta.data
  
  meta_sample_id <- left_join(meta,comb_meta, by = "Row.names")%>%column_to_rownames("Row.names")
  
  meta_sample_id$cluster_annotations <- as.character(meta_sample_id$cluster_annotations)
  for(row in rownames(meta_sample_id)){
    if(is.na(meta_sample_id[row,]$cluster_annotations)) {
      meta_sample_id[row,]$cluster_annotations <- meta_sample_id[row,]$noncancer_celltype.x
    }
  }
  annotated_objects[[unique(object$sample_id)]]@meta.data <- meta_sample_id
}

merged_obj <- purrr::reduce(annotated_objects,merge)

Idents(merged_obj) <- "cluster_annotations"

#Analysis only with core, edge, and ecm myCAFs

merged_obj <- subset(merged_obj, idents = c("core","edge","ecm.myCAF"))
merged_obj@images <- list()

sp_cc_full <- createCellChat(merged_obj)
sp_cc_full@DB <- CellChatDB.use
sp_cc_full <- subsetData(sp_cc_full)
sp_cc_full <- identifyOverExpressedGenes(sp_cc_full)
sp_cc_full <- identifyOverExpressedInteractions(sp_cc_full)

sp_cc_full <- computeCommunProb(sp_cc_full)
sp_cc_full <- filterCommunication(sp_cc_full, min.cells = 10)
sp_cc_full <- computeCommunProbPathway(sp_cc_full)
sp_cc_full <- aggregateNet(sp_cc_full)
cellchat <- sp_cc_full
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")


count_interactions <- as.data.frame(cellchat@net[["count"]])
count_interactions <- melt(count_interactions)
library(dplyr)
count_interactions <- count_interactions[c(1,5,9),]

gg1 <- ggplot(count_interactions, aes(x=variable, y=value, fill=variable))+
  geom_bar(stat="identity", color="black")+
  scale_fill_manual(values=c("core" = "#4DBBD5FF", 
                             "edge" = "#E64B35FF",
                             "ecm.myCAF" = "#3f908d"))+
  labs(y = "Number of intracellular interactions") +
  theme_minimal()

color_map = c("core" = "#4DBBD5FF", 
              "edge" = "#E64B35FF",
              "ecm.myCAF" = "#3f908d")

#intercellular interaction strength

weight <- as.data.frame(cellchat@net[["weight"]])

weight_interactions <- melt(weight)
library(dplyr)
weight_interactions <- weight_interactions[c(1,5,9),]

gg2 <- ggplot(weight_interactions, aes(x=variable, y=value, fill=variable))+
  geom_bar(stat="identity", color="black")+
  scale_fill_manual(values=c("core" = "#4DBBD5FF", 
                             "edge" = "#E64B35FF",
                             "ecm.myCAF" = "#3f908d"))+
  labs(y = "Intracellular interaction strength") + 
  theme_minimal()

png("intracellular_weighted_caf.png", units="in", width=10, height=5, res=300, type = "cairo")
gg1 + gg2
dev.off()

#count of interactions

count_interactions <- as.data.frame(cellchat@net[["count"]])
count_interactions <- melt(count_interactions)
library(dplyr)
count_interactions <- count_interactions[c(4,5,9),]
count_interactions[1,] <- c('ecm.myCAF',count_interactions[1,]$value)

gg1 <- ggplot(count_interactions, aes(x=variable, y=value, fill=variable))+
  geom_bar(stat="identity", color="black")+
  scale_fill_manual(values=c("core" = "#4DBBD5FF", 
                             "edge" = "#E64B35FF",
                             "ecm.myCAF" = "#3f908d"))+
  labs(y = "cell-cancer interactions") +
  theme_minimal()

#
weight <- as.data.frame(cellchat@net[["weight"]])

weight_interactions <- melt(weight)
library(dplyr)
weight_interactions <- weight_interactions[c(4,5,9),]
weight_interactions[1,] <- c('ecm.myCAF',weight_interactions[1,]$value)

gg2 <- ggplot(weight_interactions, aes(x=variable, y=value, fill=variable))+
  geom_bar(stat="identity", color="black")+
  scale_fill_manual(values=c("core" = "#4DBBD5FF", 
                             "edge" = "#E64B35FF",
                             "ecm.myCAF" = "#3f908d"))+
  labs(y = "cell-cancer interaction strength") + 
  theme_minimal()

png("cell_to_cancer.png", units="in", width=10, height=5, res=300, type = "cairo")
gg1 + gg2
dev.off()

color_map = c("core" = "#4DBBD5FF", 
              "edge" = "#E64B35FF",
              "ecm.myCAF" = "#3f908d")

#cell signaling from edge cells to edge cells
png("e2e.png", units="in", width=20, height=10, res=300, type = "cairo")
netVisual_chord_gene(cellchat, sources.use = c(2), targets.use = c(2),legend.pos.x = 8,
                     color.use = c("edge" = "#E64B35FF"), lab.cex = 1.2)
dev.off()

#cell signaling from core cells to core cells
png("c2c.png", units="in", width=20, height=10, res=300, type = "cairo")
netVisual_chord_gene(cellchat, sources.use = c(3), targets.use = c(3),legend.pos.x = 8, lab.cex = 1.2,
                     color.use = c("core" = "#4DBBD5FF"))
dev.off()

#cell signaling between CAF and edge cells
png("caf_and_edge.png", units="in", width=20, height=10, res=300, type = "cairo")
netVisual_chord_gene(cellchat, sources.use = c(1), targets.use = c(1,2),legend.pos.x = 8, lab.cex = 1.2,
                     color.use = c("edge" = "#E64B35FF",
                                   "ecm.myCAF" = "#3f908d"))
dev.off()

#subset for only top ligand recceptor pairs 
object <- cellchat
prob <- slot(object, "net")$prob
pval <- slot(object, "net")$pval
prob[pval > 0.05] <- 0
net <- reshape2::melt(prob, value.name = "prob")
colnames(net)[1:3] <- c("source", "target", "interaction_name")
pairLR = dplyr::select(object@LR$LRsig, c("interaction_name_2", 
                                          "pathway_name", "ligand", "receptor", "annotation", 
                                          "evidence"))
idx <- match(net$interaction_name, rownames(pairLR))
temp <- pairLR[idx, ]
net <- cbind(net, temp)

net <- net[net$prob > 1.0e-03,]

# #cell signaling between CAF and edge cells with a more stringent cutoff
# png("caf_and_edge_sub.png", units="in", width=20, height=10, res=300, type = "cairo")
# netVisual_chord_gene(object, net = net, sources.use = c(1), targets.use = c(1,2),legend.pos.x = 8,
#                      color.use = c("edge" = "#E64B35FF",
#                                    "ecm.myCAF" = "#3f908d"), thresh = 0.001,lab.cex = 1.2)
# dev.off()


# png("e2e_sub.png", units="in", width=20, height=10, res=300, type = "cairo")
# netVisual_chord_gene(object, net = net, sources.use = c(2), targets.use = c(2),legend.pos.x = 8,
#                      color.use = c("edge" = "#E64B35FF"), thresh = 0.001,lab.cex = 1.2)
# dev.off()

#cell signaling between core and core cells with a more stringent cutoff of 0.001
png("c2c_0.001.png", units="in", width=20, height=10, res=300, type = "cairo")
netVisual_chord_gene(object, net = net, sources.use = c(3), targets.use = c(3),legend.pos.x = 8,
                     color.use = c("core" = "#4DBBD5FF"), thresh = 0.001,lab.cex = 1.2)
dev.off()

net <- net[net$prob > 5.0e-03,]
#cell signaling between edge and edge cells with a more stringent cutoff of 0.005
png("e2e_0.005.png", units="in", width=20, height=10, res=300, type = "cairo")
netVisual_chord_gene(object, net = net, sources.use = c(2), targets.use = c(2),legend.pos.x = 8,
                     color.use = c("edge" = "#E64B35FF"), thresh = 0.001,lab.cex = 1.2)
dev.off()

# net <- net[net$prob > 1e-02,]
# 
# png("caf_and_edge_0.01.png", units="in", width=20, height=10, res=300, type = "cairo")
# netVisual_chord_gene(object, net = net, sources.use = c(1), targets.use = c(1,2),legend.pos.x = 8,
#                      color.use = c("edge" = "#E64B35FF",
#                                    "ecm.myCAF" = "#3f908d"), thresh = 0.001,lab.cex = 1.2)
# dev.off()

net <- net[net$prob > 5e-02,]
#cell signaling between edge and edge cells with a more stringent cutoff of 0.05
png("caf_and_edge_0.05.png", units="in", width=20, height=10, res=300, type = "cairo")
netVisual_chord_gene(object, net = net, sources.use = c(1), targets.use = c(1,2),legend.pos.x = 8,
                     color.use = c("edge" = "#E64B35FF",
                                   "ecm.myCAF" = "#3f908d"), thresh = 0.001,lab.cex = 1.2)
dev.off()

# extract signaling info
object <- cellchat
prob <- slot(object, "net")$prob
pval <- slot(object, "net")$pval
prob[pval > 0.05] <- 0
net <- reshape2::melt(prob, value.name = "prob")
colnames(net)[1:3] <- c("source", "target", "interaction_name")
pairLR = dplyr::select(object@LR$LRsig, c("interaction_name_2", 
                                          "pathway_name", "ligand", "receptor", "annotation", 
                                          "evidence"))
idx <- match(net$interaction_name, rownames(pairLR))
temp <- pairLR[idx, ]
net <- cbind(net, temp)

write.csv(net, file ="cellchat_signaling.csv")

#analyze only core and edge cells

merged_obj <- subset(merged_obj, idents = c("core", 
                                            "edge"))

sp_cc_full <- createCellChat(merged_obj)
sp_cc_full@DB <- CellChatDB.use
sp_cc_full <- subsetData(sp_cc_full)
sp_cc_full <- identifyOverExpressedGenes(sp_cc_full)
sp_cc_full <- identifyOverExpressedInteractions(sp_cc_full)

sp_cc_full <- computeCommunProb(sp_cc_full)
sp_cc_full <- filterCommunication(sp_cc_full, min.cells = 10)
sp_cc_full <- computeCommunProbPathway(sp_cc_full)
sp_cc_full <- aggregateNet(sp_cc_full)
cellchat <- sp_cc_full
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

color_map = c( 
  "edge" = "#E64B35FF", "core" = "#4DBBD5FF")

png("signaling_pathways_onlycancercell.png", units="in", width=15, height=7, res=300, type = "cairo")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", color.use = color_map,
                                         width = 3)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", color.use = color_map,
                                         width = 3)
ht1 + ht2
dev.off()

#analyze other celltypes and how they interact with core/edge cancer cells

rm(list=setdiff(ls(), "annotated_objects"))
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB

merged_obj <- purrr::reduce(annotated_objects,merge)
merged_obj@images <- list()
Idents(merged_obj) <- "cluster_annotations"

#cellchat  with core cancer cells and macrophages
core_obj <- subset(merged_obj, idents = c("core","macrophage"), invert = FALSE)

core_cc <- createCellChat(core_obj)
core_cc@DB <- CellChatDB.use
core_cc <- subsetData(core_cc)
core_cc <- identifyOverExpressedGenes(core_cc)
core_cc <- identifyOverExpressedInteractions(core_cc)

core_cc <- computeCommunProb(core_cc)
core_cc <- filterCommunication(core_cc, min.cells = 10)
core_cc <- computeCommunProbPathway(core_cc)
core_cc <- aggregateNet(core_cc)
core_cc <- netAnalysis_computeCentrality(core_cc, slot.name = "netP")

#make a plot for core signaling to macrophages, in prominent core pathways
netVisual_individual(core_cc, signaling = c(core_cc@netP$pathways), 
                     pairLR.use = c(2), layout = "chord")

png("other_celltypes/core_mac_demosome.png", units="in", width=6, height=5, res=300, type = "cairo")
netVisual_aggregate(core_cc, signaling = c("DESMOSOME"), layout = "chord",color.use = c("core" = "#4DBBD5FF",
                                                                                        "macrophage" = "#cfa42d"))
dev.off()

png("other_celltypes/core_mac_cdh.png", units="in", width=6, height=5, res=300, type = "cairo")
netVisual_aggregate(core_cc, signaling = c("CDH"), layout = "chord",color.use = c("core" = "#4DBBD5FF",
                                                                                  "macrophage" = "#cfa42d"))
dev.off()


#cellchat  with edge cancer cells, macropahges, ecmmyCAFs and cytotoxic CD8T cells
edge_obj <- subset(merged_obj, idents = c("edge", "macrophage","cytotoxic.CD8..T.",
                                          "ecm.myCAF"))
edge_cc <- createCellChat(edge_obj)
edge_cc@DB <- CellChatDB.use
edge_cc <- subsetData(edge_cc)
edge_cc <- identifyOverExpressedGenes(edge_cc)
edge_cc <- identifyOverExpressedInteractions(edge_cc)

edge_cc <- computeCommunProb(edge_cc)
edge_cc <- filterCommunication(edge_cc, min.cells = 10)
edge_cc <- computeCommunProbPathway(edge_cc)
edge_cc <- aggregateNet(edge_cc)
edge_cc <- netAnalysis_computeCentrality(edge_cc, slot.name = "netP")

#make a plot for edge signaling to macrophages, CD8 T cells, and ECM myCAFs, in prominent edge pathways
png("other_celltypes/edge_APP.png", units="in", width=10, height=10, res=300, type = "cairo")
netVisual_aggregate(edge_cc, signaling = c("APP"), layout = "chord",color.use = c("edge" = "#E64B35FF",
                                                                                        "macrophage" = "#cfa42d",
                                                                                       "cytotoxic.CD8..T." = "#79577c",
                                                                                       "ecm.myCAF" = "#3f908d"))
dev.off()

png("other_celltypes/edge_LAMININ.png", units="in", width=10, height=10, res=300, type = "cairo")
netVisual_aggregate(edge_cc, signaling = c("LAMININ"), layout = "chord",color.use = c("edge" = "#E64B35FF",
                                                                                  "macrophage" = "#cfa42d",
                                                                                  "cytotoxic.CD8..T." = "#79577c",
                                                                                  "ecm.myCAF" = "#3f908d"))
dev.off()

png("other_celltypes/edge_COLLAGEN.png", units="in", width=10, height=10, res=300, type = "cairo")
netVisual_aggregate(edge_cc, signaling = c("COLLAGEN"), layout = "chord",color.use = c("edge" = "#E64B35FF",
                                                                                      "macrophage" = "#cfa42d",
                                                                                      "cytotoxic.CD8..T." = "#79577c",
                                                                                      "ecm.myCAF" = "#3f908d"))
dev.off()
