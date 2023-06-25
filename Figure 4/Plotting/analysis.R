library(Seurat)
nm <- list.files(path="processed_data/")

df_full <- data.frame()

#create a full 3x13 figure

for (file in nm) { 
  load(paste("processed_data/",file, sep = ""))
  df <- query@meta.data
  df <- data.frame(df$scpred_prediction,stringr::str_split(file,".R")[[1]][1])
  df_full <- rbind(df_full, df)
}

names(df_full) = c("core_edge","dataset")
library(ggplot2)

library(tidyr)
library(dplyr)
df_full$core_edge <- factor(df_full$core_edge, levels = c("core","transitory","edge","noncancer"))

library(tidyr)
library(dplyr)
df_count <- df_full%>%
  group_by( dataset,core_edge) %>%
  summarize(count = n()) 

png("scpred.png",units = "in", res = 300, width = 8, height = 5, type = "cairo")
ggplot(df_count,aes(x = dataset, y = count, fill = core_edge)) +
  geom_bar(position="fill", stat="identity")+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  guides(fill = guide_legend(title = 'cell type')) +
  scale_fill_manual(values = c("core" = "#4DBBD5FF",
                               "transitory" = "#F9E076", "edge" = "#E64B35FF", "noncancer" = "lightgray"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#cluster on composition simillarity 
library(tidyverse)
library(ggdendro)
library(vegan)
library(colorspace)
library(cowplot)
t <- reshape2::dcast(df_count, dataset ~ core_edge, fun.aggregate = sum) 
t <- t %>%
  group_by(dataset) %>% 
  mutate(sum = as.character(rowSums(select(cur_data(), is.numeric)))) %>%
  summarise_if(is.numeric, ~ . / as.numeric(sum)) %>% tibble::column_to_rownames("dataset")
hc=hclust(dist(t),method="ward.D2")
hc=reorder(hc,wts=-as.matrix(t)%*%seq(ncol(t))^2) # vegan::reorder.hclust
tree=ggdendro::dendro_data(as.dendrogram(hc),type="rectangle")

p1=ggplot(ggdendro::segment(tree))+
  geom_segment(aes(x=y,y=x,xend=yend,yend=xend),lineend="round",size=.4)+
  scale_x_continuous(expand=expansion(add=c(0,.01)))+ # don't crop half of line between top-level nodes
  scale_y_continuous(limits=.5+c(0,nrow(t)),expand=c(0,0))+
  theme(
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    axis.ticks.length=unit(0,"pt"), # remove extra space occupied by ticks
    axis.title=element_blank(),
    panel.background=element_rect(fill="white"),
    panel.grid=element_blank(),
    plot.margin=margin(5,5,5,0)
  )

t=t[hc$labels[hc$order],]
t2=data.frame(V1=rownames(t)[row(t)],V2=colnames(t)[col(t)],V3=unname(do.call(c,t)))
lab=round(100*t2$V3)
lab[lab==0]=""

t2$V2 <- factor(t2$V2, levels = c("core","transitory","edge","noncancer"))
p2=ggplot(t2,aes(x=factor(V1,level=rownames(t)),y=V3,fill=V2))+
  geom_bar(stat="identity",width=1,position=position_fill(reverse=T))+
  geom_text(aes(label=lab),position=position_stack(vjust=.5,reverse=T),size=3.5)+
  coord_flip()+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  scale_fill_manual(values = c("core" = "#4DBBD5FF",
                               "transitory" = "#F9E076", "edge" = "#E64B35FF", "noncancer" = "lightgray"))+
  theme(
    axis.text=element_text(color="black",size=11),
    axis.text.x=element_blank(),
    axis.ticks=element_blank(),
    axis.title=element_blank(),
    legend.position="none",
    plot.margin=margin(5,0,5,5)
  )

png("scpred_clustered.png",units = "in", res = 300, width = 7, height = 10, type = "cairo")
cowplot::plot_grid(p2,p1,rel_widths=c(1,.2))
dev.off()

# #cluster plot for only core, transitory, and edge
# library(tidyverse)
# library(ggdendro)
# library(vegan)
# library(colorspace)
# library(cowplot)
# df_full_cc <- subset(df_full, df_full$core_edge %in% c("core","edge","transitory"))
# df_count_cc <- df_full_cc%>%
#   group_by( dataset,core_edge) %>%
#   summarize(count = n()) 
# 
# t <- reshape2::dcast(df_count_cc, dataset ~ core_edge, fun.aggregate = sum) 
# t <- t %>%
#   group_by(dataset) %>% 
#   mutate(sum = as.character(rowSums(select(cur_data(), is.numeric)))) %>%
#   summarise_if(is.numeric, ~ . / as.numeric(sum)) %>% tibble::column_to_rownames("dataset")
# hc=hclust(dist(t),method="ward.D2")
# hc=reorder(hc,wts=-as.matrix(t)%*%seq(ncol(t))^2) # vegan::reorder.hclust
# tree=ggdendro::dendro_data(as.dendrogram(hc),type="rectangle")
# 
# p1=ggplot(ggdendro::segment(tree))+
#   geom_segment(aes(x=y,y=x,xend=yend,yend=xend),lineend="round",size=.4)+
#   scale_x_continuous(expand=expansion(add=c(0,.01)))+ # don't crop half of line between top-level nodes
#   scale_y_continuous(limits=.5+c(0,nrow(t)),expand=c(0,0))+
#   theme(
#     axis.text=element_blank(),
#     axis.ticks=element_blank(),
#     axis.ticks.length=unit(0,"pt"), # remove extra space occupied by ticks
#     axis.title=element_blank(),
#     panel.background=element_rect(fill="white"),
#     panel.grid=element_blank(),
#     plot.margin=margin(5,5,5,0)
#   )
# 
# t=t[hc$labels[hc$order],]
# t2=data.frame(V1=rownames(t)[row(t)],V2=colnames(t)[col(t)],V3=unname(do.call(c,t)))
# lab=round(100*t2$V3)
# lab[lab==0]=""
# 
# t2$V2 <- factor(t2$V2, levels = c("core","transitory","edge","noncancer"))
# p2=ggplot(t2,aes(x=factor(V1,level=rownames(t)),y=V3,fill=V2))+
#   geom_bar(stat="identity",width=1,position=position_fill(reverse=T))+
#   geom_text(aes(label=lab),position=position_stack(vjust=.5,reverse=T),size=3.5)+
#   coord_flip()+
#   scale_x_discrete(expand=c(0,0))+
#   scale_y_discrete(expand=c(0,0))+
#   scale_fill_manual(values = c("core" = "#4DBBD5FF",
#                                "transitory" = "#F9E076", "edge" = "#E64B35FF", "noncancer" = "lightgray"))+
#   theme(
#     axis.text=element_text(color="black",size=11),
#     axis.text.x=element_blank(),
#     axis.ticks=element_blank(),
#     axis.title=element_blank(),
#     legend.position="none",
#     plot.margin=margin(5,0,5,5)
#   )
# 
# cowplot::plot_grid(p2,p1,rel_widths=c(1,.4))


#plot classification of every single sample 
library(ggplot2)
library(stringr)
p <- list()
i =1 

for (file in nm) { 
  load(paste("processed_data/",file, sep = ""))
  
  p1 = SpatialDimPlot(query, alpha = c(0,0)) +theme(legend.position="none") + ggtitle(label = strsplit(file[1], ".", fixed = TRUE)[[1]][1])
  
  p2 = SpatialDimPlot(query, image.alpha = 0.5, group.by = "scpred_prediction", pt.size.factor = 2,
                      cols = c("core" = "#4DBBD5FF",
                               "transitory" = "#F9E076", "edge" = "#E64B35FF", "noncancer" = "lightgray"), stroke = 0) +theme(legend.position="none")
  
  #find clusters
  query <- FindVariableFeatures(query)
  query <- ScaleData(query)
  query <- RunPCA(query)
  query <- RunUMAP(query,dims = 1:15,reduction = "pca")
  p3 = DimPlot(query, group.by = "scpred_prediction",
               cols = c("core" = "#4DBBD5FF",
                        "transitory" = "#F9E076", "edge" = "#E64B35FF", "noncancer" = "lightgray"))+theme(legend.position="none")+ theme(
                          plot.title = element_blank(),axis.title.x = element_blank(),
                          axis.title.y = element_blank())+ 
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) + 
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
  all_plots <- p1 + p2 + p3
  p[[i]] = all_plots
  i = i + 1
  
  
}

library(gridExtra)
library(patchwork)
png("scpred_full.png",units = "in", res = 300, width = 12, height = 20, type = "cairo")
wrap_plots(p, ncol = 3)
dev.off()

#plot some examples
png("scpred_example_cscc.png",units = "in", res = 300, width = 5, height = 4, type = "cairo")
p[[20]]
dev.off()

png("scpred_example_cervicalscc.png",units = "in", res = 300, width = 5, height = 4, type = "cairo")
p[[1]]
dev.off()

png("scpred_example_melanoma.png",units = "in", res = 300, width = 5, height = 4, type = "cairo")
p[[18]]
dev.off()

png("scpred_example_colorectal.png",units = "in", res = 300, width = 5, height = 4, type = "cairo")
p[[4]]
dev.off()
