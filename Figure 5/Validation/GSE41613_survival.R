library(GEOquery)
library(tidyr)
GSE41613 <- GEOquery::getGEO(GEO = "GSE41613")

#process this dataset using limma package

GSE41613_meta <- GSE41613$GSE41613_series_matrix.txt.gz@phenoData@data

GSE41613_survival <- GSE41613_meta[,colnames(GSE41613_meta)%in% c("characteristics_ch1.4","characteristics_ch1.5")]
names(GSE41613_survival) <- c("survival","fu_time")

GSE41613_expr <- GSE41613$GSE41613_series_matrix.txt.gz@assayData$exprs
GSE41613_ex_info <- GSE41613$GSE41613_series_matrix.txt.gz@featureData@data
library(limma)

# perform background correction on fluorescent intensities
#GSE41613_expr <- backgroundCorrect(GSE41613_expr, method = 'normexp',offset=0)
GSE41613_expr <- avereps(GSE41613_expr,ID = GSE41613_ex_info$`Gene Symbol`)
rownames(GSE41613_expr) <- lapply(rownames(GSE41613_expr), function(x) stringr::str_split(x,"///")[[1]][1])
rownames(GSE41613_expr) <- lapply(rownames(GSE41613_expr), function(x) gsub(" ", "", x))


GSE41613_survival$time <- lapply(GSE41613_survival$fu_time, function(x) gsub("fu time:", "", x) )

GSE41613_survival$time <- as.numeric(GSE41613_survival$time)

GSE41613_survival$os <- ifelse(GSE41613_survival$survival == "vital: Alive", 0,1)
GSE41613_survival$dss <- ifelse(GSE41613_survival$survival == "vital: Alive", 0,ifelse(
  GSE41613_survival$survival == "vital: Dead-non OC",NA,1
))


###
library(singscore)
indication <- "HNSC"

core_genes <- readRDS("core_genes.RDS")
core_genes <- c(core_genes)

edge_genes <- readRDS("edge_genes.RDS")
edge_genes <- c(edge_genes)

geneset_of_interest <- list(core_genes,edge_genes)
names(geneset_of_interest) <- c("core","edge")

all_genes <- union(core_genes,edge_genes)

rankData <- rankGenes(as.data.frame(GSE41613_expr))

#add edge and croe scores
scoredf_core <- simpleScore(rankData, upSet = core_genes , knownDirection = F)
scoredf_edge <- simpleScore(rankData, upSet = edge_genes , knownDirection = F)

#calculate gene sets
names(scoredf_core) <- c("core","core_dispersion")
names(scoredf_edge) <- c("edge","edge_dispersion")

scores_df_singshot <- cbind(scoredf_core,scoredf_edge)

scores_df_singshot <- cbind(sample = rownames(scores_df_singshot), data.frame(scores_df_singshot, row.names=NULL))
rownames(scores_df_singshot) = scores_df_singshot$sample

merged_data = scores_df_singshot %>% 
  merge(GSE41613_survival, by = 0) %>% 
  dplyr::select(Row.names,names(geneset_of_interest), time, os,dss) 

#survival testing
library(survminer)
library(survival)
df <- data.frame (geneset  = rep(names(geneset_of_interest),2),
                  indication  = rep(c("os","dss"),each = 2),
                  p_val = 1,worse_outcome = 0, hazard_ratio = 0
)

for (survival_method in c("os","dss")) {
##switch stuff here more
for (gene in names(geneset_of_interest)) {
  
  cut = surv_cutpoint(data = merged_data, time = "time", event = survival_method, variables = names(geneset_of_interest),
                      minprop = 0.1)
  res.cat <- surv_categorize(cut)
  
  res.cat$binary = res.cat[[gene]]
  
  row = df[df$gene == gene & df$indication == survival_method,]
  rownumber = row.names(df[df$gene == gene & df$indication == survival_method,])
  
  res.cat$binary <- factor(res.cat$binary, levels = c("high","low"))
  res.cat$status <- res.cat[[survival_method]]
  fit.coxph <- coxph(Surv(time, status) ~ binary, data = res.cat)
  
  HR <- round(exp(coef(fit.coxph)), 2)[[1]]
  
  if (is.na(HR)) {
    row$hazard_ratio <- 1
    row$worse_outcome = "n/a"
  } else {
    if (coef(fit.coxph) > 0) {
      row$hazard_ratio <- round(exp(coef(fit.coxph)), 2)[[1]]
      row$worse_outcome = "Low"
    } else {
      row$hazard_ratio <- 1/round(exp(coef(fit.coxph)), 2)[[1]]
      row$worse_outcome = "High"
    }
  }
  
  #get p val
  
  p_value = summary(fit.coxph)$waldtest[[3]]
  
  row$p_val = p_value
  
  df[rownumber,] = row

}
}


df$sig = 1
df[is.na(df)] = 1


for (data_row in rownames(df)) {
  row = df[data_row,]
  if(row$p_val > 0.05) {
    row$sig = "ns"
    if(row$hazard_ratio == "Inf" || row$hazard_ratio > 100) {
      row$hazard_ratio = 1
    }
  } else {
    if(row$hazard_ratio == "Inf" || row$hazard_ratio > 100) {
      row$sig = row$worse_outcome 
      row$hazard_ratio = 1
    } else {
      row$sig = row$worse_outcome 
    }
  }
  df[data_row,] = row
}

dir.create("heatmap_GSE41613")
setwd("heatmap_GSE41613")

png("survival_dotplot.png", units="cm", width=9, height=7, res=300, type = "cairo")
ggplot(df, aes(x = geneset, y = indication, color = sig)) +
  geom_point(aes(size = hazard_ratio)) +
  scale_size(range = c(7, 10))+
  scale_colour_manual(values = c("High" = "#d73027",
                                 "Low" = "#4575b4",
                                 "ns" = "lightgray")) +
  labs(x="geneset", y="survival type",
       title="GSE41613") + 
  theme_minimal()
dev.off()
setwd("..")


