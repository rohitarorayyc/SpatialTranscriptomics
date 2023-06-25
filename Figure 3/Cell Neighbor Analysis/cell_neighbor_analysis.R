library(Seurat)
library(dplyr)
library(tibble)
library(cowplot)
library(patchwork)
load(file = "annotated_objects.Robj")
load(file = "comb.Robj")

comb_meta <- comb@meta.data

#transfer annotations over from combined object to individual objects with spatial positions
#in similar fashion to degs script 

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

#now find neighboring cells and quanitfy neighbors for core cells and neighbors for edge cells 

annotations_bcs<- list()
i = 1
for (object in annotated_objects) {

#check if a spot is an edge spot,# if the spot is an edge spot check its neighbors,note down the barcodes
#of its neighbors, quantify the neighboring cells based on the barcode matrix

BC_mtx_neighbors <- object@meta.data[c("sample_id.x","cluster_annotations")]


offsets <- data.frame(x.offset=c(-2, 2, -1,  1, -1, 1),
                      y.offset=c( 0, 0, -1, -1,  1, 1))

spot.positions <- object@images$tumor@coordinates[,c("row","col")]
spot.positions$spot.idx <- rownames(spot.positions)

## Compute coordinates of each possible spot neighbor
neighbor.positions <- merge(spot.positions, offsets)
neighbor.positions$x.pos <- neighbor.positions$col + neighbor.positions$x.offset
neighbor.positions$y.pos <- neighbor.positions$row + neighbor.positions$y.offset

## Select spots that exist at neighbor coordinates
neighbors <- merge(as.data.frame(neighbor.positions), 
                   as.data.frame(spot.positions), 
                   by.x=c("x.pos", "y.pos"), by.y=c("col", "row"),
                   suffixes=c(".primary", ".neighbor"),
                   all.x=TRUE)

neighbors_anno <- neighbors[complete.cases(neighbors), ]

neighbors_anno$celltype_primary <- object@meta.data[neighbors_anno$spot.idx.primary,]$cluster_annotations
neighbors_anno$celltype_neighbor <- object@meta.data[neighbors_anno$spot.idx.neighbor,]$cluster_annotations
annotations_bcs[[i]] <- neighbors_anno
#based on neighbor graph,identify neighboring celltypes for each cell state
object@meta.data$neighbor_clust_anno <- "none"
object@meta.data$neighbor_bcs <- "none"

for(bc in rownames(object@meta.data)) {
  nbs <- neighbors[neighbors$spot.idx.primary == bc,]
  nb_bc <- na.omit(nbs$spot.idx.neighbor)
  celltypes <- object@meta.data[nb_bc,]$cluster_annotations
  
  object@meta.data[bc,]$neighbor_clust_anno <- paste(celltypes, collapse = ",")
  object@meta.data[bc,]$neighbor_bcs <- paste(nb_bc, collapse = ",")
}
annotated_objects[[i]] <- object
i =i+1

}

nb_list <- list()
for (object in annotated_objects) {
  # Convert to a data frame
  ct_nbs_df <- object@meta.data[,c("cluster_annotations","neighbor_clust_anno", "neighbor_bcs")]
  
  all_ct_nbs <- strsplit(ct_nbs_df$neighbor_clust_anno, ",")
  all_ct_nbs_barcodes <- strsplit(ct_nbs_df$neighbor_bcs, ",")
  
  neigh_df <- data.frame(primary = rep(ct_nbs_df$cluster_annotations, sapply(all_ct_nbs, length)),
                         neighbor = unlist(all_ct_nbs), neighbor_barcodes = unlist(all_ct_nbs_barcodes))
  
  cne_neighbors <- neigh_df[neigh_df$primary %in% c("core","edge"),]
  `%!in%` = Negate(`%in%`)
  cne_neighbors <- as.data.frame(cne_neighbors[cne_neighbors$neighbor %!in% c("core","edge","transitory"),])
  
  cne_neighbors <- unique(subset(cne_neighbors, !duplicated(cne_neighbors)))
  
  cne_nbs <- as.data.frame(table(cne_neighbors$primary, cne_neighbors$neighbor))
  
  colnames(cne_nbs) <- c("primary","neighbor","#neighbors")
  `%!in%` = Negate(`%in%`)
  
  cne_nbs$sample <- unique(object$sample_id.x)
  
  nb_list <- append(nb_list, list(cne_nbs), 0)
}

allc2c_net <- plyr::ldply(nb_list)

library(rstatix)

signif_abbreviation <- function(p_values) {
  sapply(p_values, function(x) {
    if (is.na(x)) {
      return(NA)
    } else if (x < 0.0001) {
      return("****")
    } else if (x < 0.001) {
      return("***")
    } else if (x < 0.01) {
      return("**")
    } else if (x < 0.05) {
      return("*")
    } else if (x > 0.05) {
      return("NS")
    } else {
      return("NS")
    }
  })
}

allc2c_net_sig <- allc2c_net[allc2c_net$neighbor %in% c("ecm.myCAF","macrophage","cytotoxic.CD8..T.","Intermediate.fibroblast"),]
pairwise_comparisons <- allc2c_net_sig %>%
  group_by(neighbor) %>%
  filter(n_distinct(primary) > 1) %>%# skip groups with only one observation
  wilcox_test(`#neighbors` ~ primary, p.adjust.method = "bonferroni")

#adjust p value
pairwise_comparisons$p.adj <- p.adjust(pairwise_comparisons$p, method = "BH")
pairwise_comparisons <- pairwise_comparisons %>%
  mutate(p.adj.signif = as.character(signif_abbreviation(p.adj)))

allc2c_net <- left_join(allc2c_net, pairwise_comparisons, by = "neighbor")

summary_data <-allc2c_net[allc2c_net$primary %in% c("core"),]

png("signaling_cellneighbor_sig_adj.png", units="in", width=12, height=5, res=300, type = "cairo")
ggplot(allc2c_net, aes(x = primary, y = `#neighbors`, color = factor(primary), fill = primary)) +
  geom_boxplot(outlier.shape = NA, aes(fill = factor(primary)), alpha = 0.5) +
  scale_fill_manual(values = c('edge' = '#E64B35FF',
                               'core' = '#4DBBD5FF')) +
  facet_wrap(~neighbor, scales = "free", nrow = 1) +
  geom_point(alpha = 1, aes(color = factor(sample), group = primary, fill = primary)) +
  scale_color_manual(values = c("sample_1" = "#66D7D1", "sample_2" = "#BEE3DB","sample_3" = "#D3D3D3",
                                "sample_4" = "#FFD6BA","sample_5" = "#274C77","sample_6" = "#A3CEF1",
                                "sample_7" = "#F7CE5B",
                                "sample_8" = "#7CEA9C","sample_9" = "#7C90A0",
                                "sample_10" = "#DCCCFF","sample_11" = "#EE6055",
                                "sample_12" = "#C492B1")) +
  ggtitle("Number of Neighboring Spots") +
  xlab("Neighbor") +
  ylab("Number of Neighbors") +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(title = "Sample")) +
  geom_text(data = summary_data, aes(x = primary, y = Inf, label = p.adj.signif), vjust = 1.5, hjust = 0.5, fontface = "bold") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),   # Rotate and adjust axis labels
        plot.margin = unit(c(1, 5, 1, 1), "lines"),
        axis.title = element_text(face = "bold"),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 10, face = "bold"),
        axis.ticks.y = element_blank())
dev.off()
