library(UCSCXenaTools)
library(dplyr)
library(survival)
library(survminer)
library(data.table)
library(ggplot2)
library(tibble)
library(stringr)
library(plyr)
library(singscore)
library(pROC)
indication <- "HNSC"
`%!in%` = Negate(`%in%`)

replacement_dictionary <- c(
  "ADIRF" = "C10orf116", "CTSV" = 'CTSL2',"ERO1A" = "ERO1L","IL36G" = "IL1F9",
  "IL36RN" = "IL1F5","ATP5F1A" = "ATP5A1","ATP5F1B" = "ATP5B", "ATP5MF" = "ATP5J2",
  'DDX39B' = "BAT1","HNRNPDL" = 'HNRPDL',"MZT2B" = "FAM128B","PKM" = 'PKM2', 
  'RACK1' = 'GNB2L1', 'SNHG29' = 'NCRNA00188','SRSF2' = "SFRS2" ,"TMA7" = 'CCDC72')

genes_to_remove <- c("DEFB4B","PRR9","RNF223","SLURP2","MIR205HG","TMSB4X")

core_genes <- readRDS("core_genes.RDS")
core_genes <- c(core_genes)
for (gene in core_genes) {
  if (gene %in% names(replacement_dictionary)) {
    core_genes <- core_genes[core_genes != gene]
    core_genes <- append(core_genes, replacement_dictionary[[gene]])
  }
} 

core_genes <- core_genes[core_genes %!in% genes_to_remove]

edge_genes <- readRDS("edge_genes.RDS")
edge_genes <- c(edge_genes)

for (gene in edge_genes) {
  if (gene %in% names(replacement_dictionary)) {
    edge_genes <- edge_genes[edge_genes != gene]
    edge_genes <- append(edge_genes, replacement_dictionary[[gene]])
  }
} 

edge_genes <- edge_genes[edge_genes %!in% genes_to_remove]

geneset_of_interest <- list(core_genes,edge_genes)
names(geneset_of_interest) <- c("core","edge")

all_genes <- union(core_genes,edge_genes)

#load in TCGA OSCC HPV negative codes
sample_ids <- read.csv(file = "TCGA_OSCC_HPV_neg_CODES.csv")
tcga_oscc_codes <- sample_ids$TCGA_codes 
tcga_oscc_sampleID <- lapply(tcga_oscc_codes, function(x) paste(x, "-01", sep = ""))

#initialize object
Xena_cohort = XenaData %>% 
  filter(XenaHostNames == "tcgaHub") %>% # select TCGA Hub
  XenaScan(indication)   # select cohort

#download clinical and survival datasets
cli_query = Xena_cohort %>% 
  filter(DataSubtype == "phenotype") %>%  # select clinical dataset
  XenaGenerate() %>%  # generate a XenaHub object
  XenaQuery() %>% 
  XenaDownload()

cli = XenaPrepare(cli_query)
survival_data = cli[[2]]
clinical_data = cli[[1]]
clinical_data$sample = clinical_data$sampleID

if ("xena_sample" %in% colnames(survival_data)) {
  survival_data$sample = survival_data$xena_sample 
  survival_data$xena_sample = NULL
}

clinical_data <- clinical_data[clinical_data$sampleID %in% tcga_oscc_sampleID, ]
survival_data <- survival_data[survival_data$sample %in% tcga_oscc_sampleID, ]  

ge <- XenaData %>%
  filter(XenaHostNames == "tcgaHub") %>% # select TCGA Hub
  XenaScan(indication) %>%
  filter(DataSubtype == "gene expression RNAseq", Label == "IlluminaHiSeq", Unit == "log2(norm_count+1)") 

gene_info_km <- data.frame()

gene_info_km <- fetch_dense_values(host = unique(ge$XenaHosts), 
                                   dataset = ge$XenaDatasets,
                                   samples = tcga_oscc_sampleID,
                                   use_probeMap = F,
                                   check = FALSE,
                                   time_limit = 10000)

gene_info_km <- as.data.frame(gene_info_km)


dir.create("tcga_plots_km")
setwd("tcga_plots_km")

###use multiple time methods
time_methods <- c("DSS.time","PFI.time","OS.time")


for(time_method in time_methods) { 
  
  rankData <- rankGenes(gene_info_km)
  
  scoredf_core <- simpleScore(rankData, upSet = core_genes , knownDirection = F)
  scoredf_edge <- simpleScore(rankData, upSet = edge_genes , knownDirection = F)
  
  #calculate gene sets
  names(scoredf_core) <- c("core","core_dispersion")
  names(scoredf_edge) <- c("edge","edge_dispersion")
  
  scores_df_singshot <- cbind(scoredf_core,scoredf_edge)
  
  scores_df_singshot <- cbind(sample = rownames(scores_df_singshot), data.frame(scores_df_singshot, row.names=NULL))
  
  merged_data = scores_df_singshot %>% 
    left_join(survival_data, by = "sample") %>% 
    dplyr::select(sample,names(geneset_of_interest), time_method, str_split(time_method,"[.]")[[1]][1])  %>% 
    left_join(clinical_data, by = "sample") %>% 
    dplyr::rename(time = time_method, 
                  status = str_split(time_method,"[.]")[[1]][1])
  
  #minprop 0.2
  cut = surv_cutpoint(data = merged_data, time = "time", event = "status", variables = c("core","edge"), minprop = 0.2)
  
  res.cat <- surv_categorize(cut)


coxph <- coxph(Surv(time, status) ~ edge, data = res.cat)
kableExtra::as_image(tab::tabcoxph(coxph), file = paste(indication, str_split(time_method,"[.]")[[1]][1],"edge.png",sep = "_"))

png(paste(indication, str_split(time_method,"[.]")[[1]][1],"edge_km.png",sep = "_"), units="in", width=7, height=5, res=300, type = "cairo")
print(ggsurvplot(survfit(Surv(time, status) ~ edge, data = res.cat), pval = T, conf.int = F,
                 palette=c("#EA738DFF","#89ABE3FF"),xscale = 365.25,surv.median.line = "v",break.time.by = 365.25,
                 xlab="Time in years",
                 ylab=paste("Probability of",str_split(time_method,"[.]")[[1]][1]),
                 surv.scale="percent"))
dev.off()


coxph <- coxph(Surv(time, status) ~ core, data = res.cat)
kableExtra::as_image(tab::tabcoxph(coxph), file = paste(indication, str_split(time_method,"[.]")[[1]][1],"core.png",sep = "_"))

png(paste(indication, str_split(time_method,"[.]")[[1]][1],"core_km.png",sep = "_"), units="in", width=7, height=5, res=300, type = "cairo")
print(ggsurvplot(survfit(Surv(time, status) ~ core, data = res.cat), pval = TRUE, conf.int = F,
      palette=c("#EA738DFF","#89ABE3FF"),xscale = 365.25,surv.median.line = "v",break.time.by = 365.25,
      xlab="Time in years",
      ylab=paste("Probability of",str_split(time_method,"[.]")[[1]][1]),
      surv.scale="percent"))
dev.off()


png(paste("figs/",indication, str_split(time_method,"[.]")[[1]][1],"edge_km.png",sep = "_"), units="in", width=6.5 , height=6, res=300, type = "cairo")
print(ggsurvplot(survfit(Surv(time, status) ~ edge, data = res.cat), pval = F, conf.int = F,risk.table = T, 
                 palette=c("#EA738DFF","#89ABE3FF"),xscale = 365.25,surv.median.line = "v",break.time.by = 365.25,
                 xlab="Time in years",
                 ylab=paste("Probability of",str_split(time_method,"[.]")[[1]][1]),
                 surv.scale="percent"))
dev.off()

png(paste("figs/",indication, str_split(time_method,"[.]")[[1]][1],"core_km.png",sep = "_"), units="in", width=6.5, height=6, res=300, type = "cairo")
print(ggsurvplot(survfit(Surv(time, status) ~ core, data = res.cat), pval = F, conf.int = F,risk.table = T, 
                 palette=c("#EA738DFF","#89ABE3FF"),xscale = 365.25,surv.median.line = "v",break.time.by = 365.25,
                 xlab="Time in years",
                 ylab=paste("Probability of",str_split(time_method,"[.]")[[1]][1]),
                 surv.scale="percent"))
dev.off()

}

setwd("..")

dir.create("other_tcga_plots")
setwd("other_tcga_plots")

library(immunedeconv)
png("rankplot_core.png", units="in", width=4, height=3, res=300, type = "cairo")
plotRankDensity(rankData[,2,drop = FALSE], upSet = core_genes, isInteractive = F) + scale_color_manual(values = c("#4DBBD5FF"))
dev.off()

png("rankplot_edge.png", units="in", width=4, height=3, res=300, type = "cairo")
plotRankDensity(rankData[,2,drop = FALSE], upSet = edge_genes, isInteractive = F) + scale_color_manual(values = c("#E64B35FF"))
dev.off()

res = deconvolute(gene_info_km, "epic")

t_res = as.data.frame(t(res))
names(t_res) <- lapply(t_res[1, ], as.character)
t_res <- t_res[-1, ] 

merged_data = scores_df_singshot %>% 
  left_join(clinical_data, by = "sample") 

rownames(merged_data) <- merged_data$sample

full_df_deconv  <- merge(t_res, merged_data, by = 0, all = TRUE)
png("core_edge_correlation.png", units="in", width=5, height=5, res=300, type = "cairo")

#correlation core edge
sp <- ggscatter(merged_data, x = "edge", y = "core",
                add = "reg.line",  # Add regression line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 0.87, label.y = 0.3)
dev.off()
library(rstatix)

###
full_df_deconv$CAF <- full_df_deconv$`Cancer associated fibroblast`
full_df_deconv$CAF <- as.numeric(full_df_deconv$CAF)
png("CAF_edge_correlation.png", units="in", width=6, height=5, res=300, type = "cairo")

#correlation edge CAF
sp <- ggscatter(full_df_deconv, x = "edge", y = "CAF",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 0.87, label.y = 0.01)
dev.off()

give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

merged_data_stage <- merged_data[merged_data$pathologic_stage %in% c("Stage I","Stage II","Stage III",
                                                                     "Stage IVA","Stage IVB"),]

png("stage_plot_core.png", units="in", width=7, height=5, res=300, type = "cairo")
ggplot(merged_data_stage, aes_string(x = "pathologic_stage", y = "core", fill = "pathologic_stage")) +
  geom_boxplot(alpha=0.7)+
  stat_summary(fun.data = give.n, geom = "text") +
  scale_y_continuous(name = "gene expression RNAseq") +
  scale_x_discrete(name = "Gene") +
  theme_bw() + 
  ggpubr::stat_compare_means()
dev.off()

png("stage_plot_edge.png", units="in", width=7, height=5, res=300, type = "cairo")
ggplot(merged_data_stage, aes_string(x = "pathologic_stage", y = "edge", fill = "pathologic_stage")) +
  geom_boxplot(alpha=0.7)+
  stat_summary(fun.data = give.n, geom = "text") +
  scale_y_continuous(name = "gene expression RNAseq") +
  scale_x_discrete(name = "Gene") +
  theme_bw() + 
  ggpubr::stat_compare_means()
dev.off()

library(RColorBrewer)

###comparison across clinical characteristics

merged_data = scores_df_singshot %>% 
  left_join(survival_data, by = "sample") %>% 
  dplyr::select(sample,names(geneset_of_interest), time_method, str_split(time_method,"[.]")[[1]][1])  %>% 
  left_join(clinical_data, by = "sample") %>% 
  dplyr::rename(time = time_method, 
                status = str_split(time_method,"[.]")[[1]][1])

edited_data <- merged_data %>%
              mutate(margin_status = replace(margin_status,margin_status %in% c("Positive"), "present"))  %>%
              mutate(pathologic_M = replace(pathologic_M,pathologic_M %in% c("MX"), "present"))  %>%
              mutate(pathologic_N = replace(pathologic_N,pathologic_N %in% c("N2a","N2b","N2","N2c"), "present")) %>%
              mutate(pathologic_T = replace(pathologic_T,pathologic_T %in% c("T3","T4","T4a"), "present")) %>%
              mutate(presence_of_pathological_nodal_extracapsular_spread = replace(presence_of_pathological_nodal_extracapsular_spread,
                                                                                   presence_of_pathological_nodal_extracapsular_spread %in% c("Gross Extension","Microscopic Extension"), "present")) %>%
              mutate(lymphovascular_invasion_present = replace(lymphovascular_invasion_present,lymphovascular_invasion_present %in% c("YES"), "present"))  %>% 
            mutate(neoplasm_histologic_grade = replace(neoplasm_histologic_grade,neoplasm_histologic_grade %in% c("G3"), "present"))  %>% 
              select(edge,core,pathologic_M,pathologic_T,pathologic_N,
                     lymphovascular_invasion_present,neoplasm_histologic_grade, margin_status,
                     presence_of_pathological_nodal_extracapsular_spread) %>%
                     melt(id.vars = c('core',"edge"), measure.vars = c("pathologic_M","pathologic_T","pathologic_N",
                                                            "lymphovascular_invasion_present","neoplasm_histologic_grade", "margin_status",
                                                            "presence_of_pathological_nodal_extracapsular_spread"))%>% na.omit
              
edited_data$value[!grepl('present', edited_data$value)] <- "absent"

###

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

wilcox_tests <- edited_data %>%
  group_by(variable) %>%
  do(tidy(wilcox.test(edge ~ value, data = .))) 

wilcox_tests$p.adj <- p.adjust(wilcox_tests$p.value, method = "BH")

wilcox_tests <- wilcox_tests %>%
  mutate(p.adj.signif = as.character(signif_abbreviation(p.adj)))

edited_data_with_p_adj <- left_join(edited_data, wilcox_tests, by = "variable")


png("comparision_plot_edge.png", units="in", width=6, height=6.5, res=300, type = "cairo")

ggplot(edited_data_with_p_adj, aes(x = variable, y = edge), fill = value)+
  geom_boxplot(outlier.shape = NA,aes(fill = factor(value)), alpha = 0.5)+
  scale_fill_manual(values = c('present' = '#Fa5d0f',
                               'absent' = '#2B1055'))+
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(title = "Sample")) +
  geom_text(data = edited_data_with_p_adj, aes(x = variable, y = Inf, label = p.adj.signif), vjust = 1.5, hjust = 0.5, fontface = "bold") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +   # Bold axis labels
  theme(plot.margin = unit(c(1, 5, 1, 1), "lines"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(size = 10, face = "bold"),
        axis.ticks.y = element_blank(),  panel.border = element_rect(linewidth = 0)) 
dev.off()

wilcox_tests <- edited_data %>%
  group_by(variable) %>%
  do(tidy(wilcox.test(core ~ value, data = .))) 

wilcox_tests$p.adj <- p.adjust(wilcox_tests$p.value, method = "BH")

wilcox_tests <- wilcox_tests %>%
  mutate(p.adj.signif = as.character(signif_abbreviation(p.adj)))

edited_data_with_p_adj <- left_join(edited_data, wilcox_tests, by = "variable")


png("comparision_plot_core.png", units="in", width=6, height=6.5, res=300, type = "cairo")
ggplot(edited_data, aes(x = variable, y = core), fill = value)+
  geom_boxplot(outlier.shape = NA,aes(fill = factor(value)), alpha = 0.5)+
  scale_fill_manual(values = c('present' = '#Fa5d0f',
                               'absent' = '#2B1055'))+
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(title = "Sample")) +
  geom_text(data = edited_data_with_p_adj, aes(x = variable, y = Inf, label = p.adj.signif), vjust = 1.5, hjust = 0.5, fontface = "bold") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +   # Bold axis labels
  theme(plot.margin = unit(c(1, 5, 1, 1), "lines"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(size = 10, face = "bold"),
        axis.ticks.y = element_blank(),  panel.border = element_rect(linewidth = 0)) 
dev.off()




