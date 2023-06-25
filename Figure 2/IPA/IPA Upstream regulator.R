library(remotes)
library(clusterProfiler)
library(rlang)
library(multienrichjam);
library(jamba)
library(amap)
#> 
#> Attaching package: 'jamba'
#> The following objects are masked from 'package:multienrichjam':
#> 
#>     heatmap_column_order, heatmap_row_order
library(colorjam);
suppressPackageStartupMessages(library(ComplexHeatmap));
options("stringsAsFactors"=FALSE, "warn"=-1);
knitr::opts_chunk$set(
  fig.height=10,
  fig.width=10,
  fig.align="center"
)
ragg_png = function(..., res = 192) {
  ragg::agg_png(..., res = res, units = "in")
}
knitr::opts_chunk$set(dev = "ragg_png", fig.ext = "png")


#Importing IPA Text file
ipa_files <- c(Patient_1_Core='patient1-core.txt',
               Patient_2_Core='patient2-core.txt',
               Patient_3_Core='patient3-core.txt',
               Patient_4_Core='patient4-core.txt',
               Patient_5_Core='patient5-core.txt',
               Patient_6_Core='patient6-core.txt',
               Patient_7_Core='patient7-core.txt',
               Patient_8_Core='patient8-core.txt',
               Patient_9_Core='patient9-core.txt',
               Patient_10_Core='patient10-core.txt',
               Patient_11_Core='patient11-core.txt',
               Patient_12_Core='patient12-core.txt',
               Patient_1_Edge='patient1-edge.txt',
               Patient_2_Edge='patient2-edge.txt',
               Patient_3_Edge='patient3-edge.txt',
               Patient_4_Edge='patient4-edge.txt',
               Patient_5_Edge='patient5-edge.txt',
               Patient_6_Edge='patient6-edge.txt',
               Patient_7_Edge='patient7-edge.txt',
               Patient_8_Edge='patient8-edge.txt',
               Patient_9_Edge='patient9-edge.txt',
               Patient_10_Edge='patient10-edge.txt',
               Patient_11_Edge='patient11-edge.txt',
               Patient_12_Edge='patient12-edge.txt')

ipa_l <- lapply(ipa_files, importIPAenrichment)
ssdim(ipa_l)


#Analyze IPA enrichments from canonical enrichment test
enrichList_canonical <- lapply(ipa_l, function(i){
  i[["Upstream Regulators"]];
});
sdim(enrichList_canonical);

## Convert data.frame to enrichResult
## multienrichjam::enrichDF2enrichResult
er_canonical <- lapply(enrichList_canonical, function(i){
  enrichDF2enrichResult(i,
                        keyColname="Name",
                        pvalueColname="P-value",
                        geneColname="geneNames",
                        pathGenes="Bias Term",
                        pvalueCutoff=1)
});
sdim(er_canonical);
head(as.data.frame(er_canonical[[1]]));


?multiEnrichMap
#Running multiEnrichMap()
mem_canonical <- multiEnrichMap(er_canonical,
                                enrichBaseline=1,
                                cutoffRowMinP=0.05,
                                colorV=c("purple", "orange"),
                                topEnrichN=20)
sdim(mem_canonical)

nameColname = c("Name", "pathway", "Description", "itemsetID", 
                "ID")
nameColname <- find_colname(nameColname, er_canonical)
#
# for upstreeam regulator use Activation z-score
enrichIM <- enrichList2IM(er_canonical, valueColname = "Activation z-score", 
                          keyColname = nameColname, verbose = F, emptyValue = 0, 
                          GmtT = NULL)

enrichIM <- enrichIM[rownames(enrichIM) %in% rownames(mem_canonical$enrichIM),]

#change enrichIM if you want to change the whole plot

enrichIM <- enrichIM[rowSums(enrichIM[]) != 0 ,]

mem_canonical$enrichIM <- mem_canonical$enrichIM[rownames(mem_canonical$enrichIM) %in% rownames(mem_canonical$enrichIM),]

binarydf <- enrichIM
binarydf <- as.data.frame(ifelse(enrichIM != 0, 1, 0))
binarydf$total <- rowSums(binarydf)
keep <- rownames(binarydf[binarydf$total>10,])
enrichIM <- enrichIM[rownames(enrichIM) %in% keep,]

use_matrix <- enrichIM

name = "Activation z-score"
col_logp <- circlize::colorRamp2(breaks = c(min(enrichIM), 0, max(enrichIM)), colors = c(rev(jamba::getColorRamp("RdBu", 
                                                                                                                 lens = 2, n = 3, trimRamp = c(2, 2)))))

er_hc2 <- amap::hcluster(link = "ward", jamba::noiseFloor(enrichIM[rownames(enrichIM), 
                                                                   , drop = FALSE], minimum = min(enrichIM), 
                                                          newValue = 0, ceiling = max(enrichIM)), method = "euclidean")
er_hc2 <- as.dendrogram(er_hc2)
row_cex = 1
row_fontsize <- jamba::noiseFloor(row_cex * 60/(nrow(enrichIM))^(1/2), 
                                  minimum = 12, ceiling = 18)
column_cex = 1

column_fontsize <- jamba::noiseFloor(column_cex * 60/(ncol(enrichIM))^(1/2), 
                                     minimum = 1, ceiling = 20)

show_heatmap_legend <- TRUE
top_annotation <- NULL 
show = NULL
cell_fun_custom <- cell_fun_bivariate(list(use_matrix, 
                                           FALSE, mem_canonical$enrichIMgeneCount), pch = 21, 
                                      size_fun = ct_approxfun, size_by = 3, outline_style = "darker", 
                                      col_hm = col_logp, show = show, cex = 0, 
                                      type = "univariate", prefix = c("z-score: ", 
                                                                      "-log10P: ", "genes: ")[show])

hm <- call_fn_ellipsis(ComplexHeatmap::Heatmap, matrix = use_matrix, 
                       name = name, col = col_logp, cluster_rows = er_hc2, 
                       row_dend_reorder = TRUE, border = TRUE, 
                       row_names_gp = grid::gpar(fontsize = row_fontsize), 
                       row_names_max_width = grid::unit(12, "cm"), column_names_gp = grid::gpar(fontsize = column_fontsize), 
                       column_names_max_height = grid::unit(12, "cm"), 
                       cluster_columns = FALSE, row_dend_width = grid::unit(30, "cm"), 
                       column_title = NULL, heatmap_legend_param = list(border = "black", legend_height = grid::unit(6, "cm")), 
                       rect_gp = grid::gpar(type = "none"), cell_fun = cell_fun_custom, 
                       show_heatmap_legend = show_heatmap_legend, top_annotation = top_annotation)

gene_count_max <- NULL
min_count = 1
if (length(gene_count_max) == 0) {
  ctmax <- ceiling(max(mem_canonical$enrichIMgeneCount, na.rm = TRUE))
} else {
  ctmax <- gene_count_max
}
if (ctmax <= 1) {
  ct_ticks <- c(0, 1)
} else {
  n <- 8
  ct_ticks <- setdiff(unique(c(min_count, round(pretty(c(0, 
                                                         ctmax), n = n)))), 0)
  ct_step <- median(diff(ct_ticks))
  if (max(ct_ticks) > ctmax) {
    ct_ticks[which.max(ct_ticks)] <- ctmax
    if (tail(diff(ct_ticks), 1) < ceiling(ct_step/4)) {
      ct_ticks <- head(ct_ticks, -2)
    }
    else if (tail(diff(ct_ticks), 1) < ceiling(ct_step/2)) {
      ct_ticks <- c(head(ct_ticks, -2), ctmax)
    }
  }
}
point_size_factor = 1
point_size_max = 8
point_size_min = 1

ct_approxfun <- function(x, ...) {
  approxfun(x = sqrt(c(min_count, ctmax)), yleft = 0, 
            ties = "ordered", yright = point_size_max, y = c(point_size_min, 
                                                             point_size_max * point_size_factor))(sqrt(x), 
                                                                                                  ...)
}
ct_tick_sizes <- ct_approxfun(ct_ticks)
pt_legend_ncol <- 1
if (length(ct_ticks) >= 8) {
  pt_legend_ncol <- 2
}

ct_tick_sizes <- ct_approxfun(ct_ticks)
pch <- 21
pt_legend <- ComplexHeatmap::Legend(labels = ct_ticks, 
                                    title = "Gene Count", type = "points", pch = pch, 
                                    ncol = pt_legend_ncol, size = grid::unit(ct_tick_sizes, 
                                                                             "mm"), grid_height = grid::unit(max(ct_tick_sizes) * 
                                                                                                               0.95, "mm"), grid_width = grid::unit(max(ct_tick_sizes) * 
                                                                                                                                                      0.95, "mm"), background = "transparent", legend_gp = grid::gpar(col = "black", 
                                                                                                                                                                                                                      fill = "grey85"))
anno_legends <- list(pt_legend)
attr(hm, "annotation_legend_list") <- anno_legends
draw(hm, annotation_legend_list = anno_legends)


cor.exp2 = data.frame(group = c(rep("core", 12), 
                                rep("edge", 12)))


list1 = list(group = c("edge" = '#E64B35FF',
                       "core" = '#4DBBD5FF'))

list2 = list(setNames(list1$group, c("edge", "core")))
names(list2) = "group"

ha = HeatmapAnnotation(df = cor.exp2, col = list2)

hm <- call_fn_ellipsis(ComplexHeatmap::Heatmap, matrix = use_matrix, 
                       name = name, col = col_logp, cluster_rows = er_hc2, 
                       row_dend_reorder = TRUE, border = TRUE, 
                       row_names_gp = grid::gpar(fontsize = row_fontsize), 
                       row_names_max_width = grid::unit(12, "cm"), column_names_gp = grid::gpar(fontsize = 14), 
                       column_names_max_height = grid::unit(12, "cm"), 
                       cluster_columns = FALSE, row_dend_width = grid::unit(2, "cm"), 
                       column_title = NULL, heatmap_legend_param = list(border = "black", legend_height = grid::unit(6, "cm")),
                       bottom_annotation = ha)


png("IPA_Heatmap_UpReg.png", units = "in", res=300, width=9.5, height=16)
draw(hm)
dev.off()
