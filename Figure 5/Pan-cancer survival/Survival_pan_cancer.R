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

#replace gene names to gene names that exist in TCGA
replacement_dictionary <- c(
  "ADIRF" = "C10orf116", "CTSV" = 'CTSL2',"ERO1A" = "ERO1L","IL36G" = "IL1F9",
  "IL36RN" = "IL1F5","ATP5F1A" = "ATP5A1","ATP5F1B" = "ATP5B", "ATP5MF" = "ATP5J2",
  'DDX39B' = "BAT1","HNRNPDL" = 'HNRPDL',"MZT2B" = "FAM128B","PKM" = 'PKM2', 
  'RACK1' = 'GNB2L1', 'SNHG29' = 'NCRNA00188','SRSF2' = "SFRS2" ,"TMA7" = 'CCDC72')

#remove genes not in TCGA
genes_to_remove <- c("DEFB4B","PRR9","RNF223","SLURP2","MIR205HG","TMSB4X")


core_genes <- readRDS("core_genes.RDS")

for (gene in core_genes) {
  if (gene %in% names(replacement_dictionary)) {
    core_genes <- core_genes[core_genes != gene]
    core_genes <- append(core_genes, replacement_dictionary[[gene]])
  }
} 

core_genes <- core_genes[core_genes %!in% genes_to_remove]

edge_genes <- readRDS("edge_genes.RDS")

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


cancer_indications <- c("BRCA","LUAD","LUSC","LIHC","OSCC","COADREAD","PAAD","ACC",
                        "BLCA","STAD","THCA","PRAD",
                        "CESC","GBMLGG","KIRP","MESO",
                        "KIRC","OV","SARC","SKCM")

#make a df to keep track of geneset, indication, pval, worse outcome, and hazard ratio
df <- data.frame (geneset  = rep(names(geneset_of_interest),length(cancer_indications)),
                  indication  = rep(cancer_indications,each = length(geneset_of_interest)),
                  p_val = 1,worse_outcome = 0, hazard_ratio = 0
)

for (indication in cancer_indications) {
  #initialize object
  if (indication != "OSCC") {
  Xena_cohort = XenaData %>% 
    filter(XenaHostNames == "tcgaHub") %>% # select TCGA Hub
    XenaScan(indication)   # select cohort
  
  #download clinical and survival datasets
  cli_query = Xena_cohort %>% 
    filter(DataSubtype == "phenotype") %>%  # select clinical dataset
    XenaGenerate() %>%  # generate a XenaHub object
    XenaQuery() %>% 
    XenaDownload()
  
  if (indication != "OV") {
    cli = XenaPrepare(cli_query)
    survival_data = cli[[2]]
    clinical_data = cli[[1]]
    clinical_data$sample = clinical_data$sampleID
  } else {
    cli = XenaPrepare(cli_query)
    survival_data = cli[[3]]
    clinical_data = cli[[1]]
    clinical_data$sample = clinical_data$sampleID
  }
  
  if ("xena_sample" %in% colnames(survival_data)) {
    survival_data$sample = survival_data$xena_sample 
    survival_data$xena_sample = NULL
  }
  
  #process exceptions
  if (indication == "OV") {
    ge <- XenaData %>%
      filter(XenaHostNames == "tcgaHub") %>% # select TCGA Hub
      XenaScan(indication) %>%
      filter(DataSubtype == "gene expression RNAseq", Label == "IlluminaHiSeq UNC", Unit == "log2(norm_count+1)") 
  } else { 
    if (indication == "STAD") {
      ge <- XenaData %>%
        filter(XenaHostNames == "tcgaHub") %>% # select TCGA Hub
        XenaScan(indication) %>%
        filter(DataSubtype == "gene expression RNAseq", Label == "IlluminaHiSeq UNC", Unit == "log2(norm_count+1)") 
      ge <- ge[2,]
    } else {
      ge <- XenaData %>%
        filter(XenaHostNames == "tcgaHub") %>% # select TCGA Hub
        XenaScan(indication) %>%
        filter(DataSubtype == "gene expression RNAseq", Label == "IlluminaHiSeq", Unit == "log2(norm_count+1)") 
    }
  }
  
  if (indication == "SARC") {
    ge = ge[1,]
  } } else {
    sample_ids <- read.csv(file = "TCGA_OSCC_HPV_neg_CODES.csv")
    tcga_oscc_codes <- sample_ids$TCGA_codes 
    tcga_oscc_sampleID <- lapply(tcga_oscc_codes, function(x) paste(x, "-01", sep = ""))
    
    #initialize object
    Xena_cohort = XenaData %>% 
      filter(XenaHostNames == "tcgaHub") %>% # select TCGA Hub
      XenaScan("HNSC")   # select cohort
    
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
      XenaScan("HNSC") %>%
      filter(DataSubtype == "gene expression RNAseq", Label == "IlluminaHiSeq", Unit == "log2(norm_count+1)") 
  }
  
  gene_info_km <- data.frame()
  
  samples <- fetch_dataset_samples( unique(ge$XenaHosts),  ge$XenaDatasets, limit = NULL)
  #fetch gene expression data
  gene_info_km <- fetch_dense_values(host = unique(ge$XenaHosts), 
                                     dataset = ge$XenaDatasets,
                                     samples = samples,
                                     use_probeMap = F,
                                     check = FALSE,
                                     time_limit = 10000)
  
  gene_info_km <- as.data.frame((gene_info_km))
  
  
  #calculate gene sets
  
  rankData <- rankGenes(gene_info_km)
  
  scoredf_core <- simpleScore(rankData, upSet = core_genes , knownDirection = FALSE)
  scoredf_edge <- simpleScore(rankData, upSet = edge_genes , knownDirection = FALSE)
  
  #calculate gene sets
  names(scoredf_core) <- c("core","core_dispersion")
  names(scoredf_edge) <- c("edge","edge_dispersion")
  
  scores_df_singshot <- cbind(scoredf_core,scoredf_edge)
  
  scores_df_singshot <- cbind(sample = rownames(scores_df_singshot), data.frame(scores_df_singshot, row.names=NULL))
  
  merged_data = survival_data %>% 
    left_join(scores_df_singshot, by = "sample") %>% 
    dplyr::select(sample,names(geneset_of_interest), OS.time, OS)  %>% 
    left_join(clinical_data, by = "sample") %>% 
    dplyr::select(sample,names(geneset_of_interest), OS.time, OS,sample_type)  %>% 
    dplyr::rename(time = OS.time, 
                  status = OS)
  
  if (indication != "LAML") {
    merged_data = merged_data[merged_data$sample_type == "Primary Tumor",]
  }
  
  merged_data <- merged_data[complete.cases(merged_data), ]
  
  #identify cutpoint using minprop of 0.1
  cut = surv_cutpoint(data = merged_data, time = "time", event = "status", variables = names(geneset_of_interest), minprop = 0.1)
  res.cat <- surv_categorize(cut)
  
  for (gene in names(geneset_of_interest)) {
    
    res.cat$binary = res.cat[[gene]]
    
    row = df[df$gene == gene & df$indication == indication,]
    rownumber = row.names(df[df$gene == gene & df$indication == indication,])
    
    res.cat$binary <- factor(res.cat$binary, levels = c("high","low"))
    
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


dir.create("heatmap_tcga")
setwd("heatmap_tcga")


png("oss_survival_dotplot.png", units="in", width=3, height=5, res=300, type = "cairo")
ggplot(df, aes(x = geneset, y = indication, color = sig)) +
  geom_point(aes(size = hazard_ratio, fill = sig)) +
  scale_colour_manual(values = c("High" = "#d73027",
                                 "Low" = "#4575b4","ns" ="grey")) +
  labs(x="geneset", y="indication",
       title="Decreased OS") + 
  theme_minimal()
dev.off()
setwd("..")

###

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
`%!in%` = Negate(`%in%`)


replacement_dictionary <- c(
  "ADIRF" = "C10orf116", "CTSV" = 'CTSL2',"ERO1A" = "ERO1L","IL36G" = "IL1F9",
  "IL36RN" = "IL1F5","ATP5F1A" = "ATP5A1","ATP5F1B" = "ATP5B", "ATP5MF" = "ATP5J2",
  'DDX39B' = "BAT1","HNRNPDL" = 'HNRPDL',"MZT2B" = "FAM128B","PKM" = 'PKM2', 
  'RACK1' = 'GNB2L1', 'SNHG29' = 'NCRNA00188','SRSF2' = "SFRS2" ,"TMA7" = 'CCDC72')

genes_to_remove <- c("DEFB4B","PRR9","RNF223","SLURP2","MIR205HG","TMSB4X")


core_genes <- readRDS("core_genes.RDS")

for (gene in core_genes) {
  if (gene %in% names(replacement_dictionary)) {
    core_genes <- core_genes[core_genes != gene]
    core_genes <- append(core_genes, replacement_dictionary[[gene]])
  }
} 

core_genes <- core_genes[core_genes %!in% genes_to_remove]

edge_genes <- readRDS("edge_genes.RDS")

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


cancer_indications <- c("BRCA","LUAD","LUSC","LIHC","OSCC","COADREAD","PAAD","ACC",
                        "BLCA","STAD","THCA","PRAD",
                        "CESC","GBMLGG","KIRP","MESO",
                        "KIRC","OV","SARC","SKCM")

#make a df row 1 is gene, row 2 indication, row 3 p value, row 4 low or high
df <- data.frame (geneset  = rep(names(geneset_of_interest),length(cancer_indications)),
                  indication  = rep(cancer_indications,each = length(geneset_of_interest)),
                  p_val = 1,worse_outcome = 0, hazard_ratio = 0
)

for (indication in cancer_indications) {
  #initialize object
  if (indication != "OSCC") {
    Xena_cohort = XenaData %>% 
      filter(XenaHostNames == "tcgaHub") %>% # select TCGA Hub
      XenaScan(indication)   # select cohort
    
    #download clinical and survival datasets
    cli_query = Xena_cohort %>% 
      filter(DataSubtype == "phenotype") %>%  # select clinical dataset
      XenaGenerate() %>%  # generate a XenaHub object
      XenaQuery() %>% 
      XenaDownload()
    
    if (indication != "OV") {
      cli = XenaPrepare(cli_query)
      survival_data = cli[[2]]
      clinical_data = cli[[1]]
      clinical_data$sample = clinical_data$sampleID
    } else {
      cli = XenaPrepare(cli_query)
      survival_data = cli[[3]]
      clinical_data = cli[[1]]
      clinical_data$sample = clinical_data$sampleID
    }
    
    if ("xena_sample" %in% colnames(survival_data)) {
      survival_data$sample = survival_data$xena_sample 
      survival_data$xena_sample = NULL
    }
    
    if (indication == "OV") {
      ge <- XenaData %>%
        filter(XenaHostNames == "tcgaHub") %>% # select TCGA Hub
        XenaScan(indication) %>%
        filter(DataSubtype == "gene expression RNAseq", Label == "IlluminaHiSeq UNC", Unit == "log2(norm_count+1)") 
    } else { 
      if (indication == "STAD") {
        ge <- XenaData %>%
          filter(XenaHostNames == "tcgaHub") %>% # select TCGA Hub
          XenaScan(indication) %>%
          filter(DataSubtype == "gene expression RNAseq", Label == "IlluminaHiSeq UNC", Unit == "log2(norm_count+1)") 
        ge <- ge[2,]
      } else {
        ge <- XenaData %>%
          filter(XenaHostNames == "tcgaHub") %>% # select TCGA Hub
          XenaScan(indication) %>%
          filter(DataSubtype == "gene expression RNAseq", Label == "IlluminaHiSeq", Unit == "log2(norm_count+1)") 
      }
    }
    
    if (indication == "SARC") {
      ge = ge[1,]
    } } else {
      sample_ids <- read.csv(file = "TCGA_OSCC_HPV_neg_CODES.csv")
      tcga_oscc_codes <- sample_ids$TCGA_codes 
      tcga_oscc_sampleID <- lapply(tcga_oscc_codes, function(x) paste(x, "-01", sep = ""))
      
      #initialize object
      Xena_cohort = XenaData %>% 
        filter(XenaHostNames == "tcgaHub") %>% # select TCGA Hub
        XenaScan("HNSC")   # select cohort
      
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
        XenaScan("HNSC") %>%
        filter(DataSubtype == "gene expression RNAseq", Label == "IlluminaHiSeq", Unit == "log2(norm_count+1)") 
    }
  
  gene_info_km <- data.frame()
  
  samples <- fetch_dataset_samples( unique(ge$XenaHosts),  ge$XenaDatasets, limit = NULL)
  
  gene_info_km <- fetch_dense_values(host = unique(ge$XenaHosts), 
                                     dataset = ge$XenaDatasets,
                                     samples = samples,
                                     use_probeMap = F,
                                     check = FALSE,
                                     time_limit = 10000)
  
  gene_info_km <- as.data.frame((gene_info_km))
  
  
  #calculate gene sets
  
  rankData <- rankGenes(gene_info_km)
  
  scoredf_core <- simpleScore(rankData, upSet = core_genes , knownDirection = FALSE)
  scoredf_edge <- simpleScore(rankData, upSet = edge_genes , knownDirection = FALSE)
  
  #calculate gene sets
  names(scoredf_core) <- c("core","core_dispersion")
  names(scoredf_edge) <- c("edge","edge_dispersion")
  
  scores_df_singshot <- cbind(scoredf_core,scoredf_edge)
  
  scores_df_singshot <- cbind(sample = rownames(scores_df_singshot), data.frame(scores_df_singshot, row.names=NULL))
  
  merged_data = survival_data %>% 
    left_join(scores_df_singshot, by = "sample") %>% 
    dplyr::select(sample,names(geneset_of_interest), DSS.time, DSS)  %>% 
    left_join(clinical_data, by = "sample") %>% 
    dplyr::select(sample,names(geneset_of_interest), DSS.time, DSS,sample_type)  %>% 
    dplyr::rename(time = DSS.time, 
                  status = DSS)
  
  if (indication != "LAML") {
    merged_data = merged_data[merged_data$sample_type == "Primary Tumor",]
  }
  
  merged_data <- merged_data[complete.cases(merged_data), ]
  
  cut = surv_cutpoint(data = merged_data, time = "time", event = "status", variables = names(geneset_of_interest), minprop = 0.1)
  res.cat <- surv_categorize(cut)
  
  ##switch stuff here more
  for (gene in names(geneset_of_interest)) {
    
    res.cat$binary = res.cat[[gene]]
    
    row = df[df$gene == gene & df$indication == indication,]
    rownumber = row.names(df[df$gene == gene & df$indication == indication,])
    
    #try high ref
    res.cat$binary <- factor(res.cat$binary, levels = c("high","low"))
    
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


dir.create("heatmap_tcga")
setwd("heatmap_tcga")


png("survival_dotplot_DSS.png", units="in", width=3, height=5, res=300, type = "cairo")
ggplot(df, aes(x = geneset, y = indication, color = sig)) +
  geom_point(aes(size = hazard_ratio, fill = sig)) +
  scale_colour_manual(values = c("High" = "#d73027",
                                 "Low" = "#4575b4","ns" ="grey")) +
  labs(x="geneset", y="indication",
       title="Decreased DSS") + 
  theme_minimal()
dev.off()
setwd("..")
###

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
`%!in%` = Negate(`%in%`)


replacement_dictionary <- c(
  "ADIRF" = "C10orf116", "CTSV" = 'CTSL2',"ERO1A" = "ERO1L","IL36G" = "IL1F9",
  "IL36RN" = "IL1F5","ATP5F1A" = "ATP5A1","ATP5F1B" = "ATP5B", "ATP5MF" = "ATP5J2",
  'DDX39B' = "BAT1","HNRNPDL" = 'HNRPDL',"MZT2B" = "FAM128B","PKM" = 'PKM2', 
  'RACK1' = 'GNB2L1', 'SNHG29' = 'NCRNA00188','SRSF2' = "SFRS2" ,"TMA7" = 'CCDC72')

genes_to_remove <- c("DEFB4B","PRR9","RNF223","SLURP2","MIR205HG","TMSB4X")


core_genes <- readRDS("core_genes.RDS")

for (gene in core_genes) {
  if (gene %in% names(replacement_dictionary)) {
    core_genes <- core_genes[core_genes != gene]
    core_genes <- append(core_genes, replacement_dictionary[[gene]])
  }
} 

core_genes <- core_genes[core_genes %!in% genes_to_remove]

edge_genes <- readRDS("edge_genes.RDS")

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


cancer_indications <- c("BRCA","LUAD","LUSC","LIHC","OSCC","COADREAD","PAAD","ACC",
                        "BLCA","STAD","THCA","PRAD",
                        "CESC","GBMLGG","KIRP","MESO",
                        "KIRC","OV","SARC","SKCM")

#make a df row 1 is gene, row 2 indication, row 3 p value, row 4 low or high
df <- data.frame (geneset  = rep(names(geneset_of_interest),length(cancer_indications)),
                  indication  = rep(cancer_indications,each = length(geneset_of_interest)),
                  p_val = 1,worse_outcome = 0, hazard_ratio = 0
)

for (indication in cancer_indications) {
  #initialize object
  if (indication != "OSCC") {
    Xena_cohort = XenaData %>% 
      filter(XenaHostNames == "tcgaHub") %>% # select TCGA Hub
      XenaScan(indication)   # select cohort
    
    #download clinical and survival datasets
    cli_query = Xena_cohort %>% 
      filter(DataSubtype == "phenotype") %>%  # select clinical dataset
      XenaGenerate() %>%  # generate a XenaHub object
      XenaQuery() %>% 
      XenaDownload()
    
    if (indication != "OV") {
      cli = XenaPrepare(cli_query)
      survival_data = cli[[2]]
      clinical_data = cli[[1]]
      clinical_data$sample = clinical_data$sampleID
    } else {
      cli = XenaPrepare(cli_query)
      survival_data = cli[[3]]
      clinical_data = cli[[1]]
      clinical_data$sample = clinical_data$sampleID
    }
    
    if ("xena_sample" %in% colnames(survival_data)) {
      survival_data$sample = survival_data$xena_sample 
      survival_data$xena_sample = NULL
    }
    
    if (indication == "OV") {
      ge <- XenaData %>%
        filter(XenaHostNames == "tcgaHub") %>% # select TCGA Hub
        XenaScan(indication) %>%
        filter(DataSubtype == "gene expression RNAseq", Label == "IlluminaHiSeq UNC", Unit == "log2(norm_count+1)") 
    } else { 
      if (indication == "STAD") {
        ge <- XenaData %>%
          filter(XenaHostNames == "tcgaHub") %>% # select TCGA Hub
          XenaScan(indication) %>%
          filter(DataSubtype == "gene expression RNAseq", Label == "IlluminaHiSeq UNC", Unit == "log2(norm_count+1)") 
        ge <- ge[2,]
      } else {
        ge <- XenaData %>%
          filter(XenaHostNames == "tcgaHub") %>% # select TCGA Hub
          XenaScan(indication) %>%
          filter(DataSubtype == "gene expression RNAseq", Label == "IlluminaHiSeq", Unit == "log2(norm_count+1)") 
      }
    }
    
    if (indication == "SARC") {
      ge = ge[1,]
    } } else {
      sample_ids <- read.csv(file = "TCGA_OSCC_HPV_neg_CODES.csv")
      tcga_oscc_codes <- sample_ids$TCGA_codes 
      tcga_oscc_sampleID <- lapply(tcga_oscc_codes, function(x) paste(x, "-01", sep = ""))
      
      #initialize object
      Xena_cohort = XenaData %>% 
        filter(XenaHostNames == "tcgaHub") %>% # select TCGA Hub
        XenaScan("HNSC")   # select cohort
      
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
        XenaScan("HNSC") %>%
        filter(DataSubtype == "gene expression RNAseq", Label == "IlluminaHiSeq", Unit == "log2(norm_count+1)") 
    }
  
  gene_info_km <- data.frame()
  
  samples <- fetch_dataset_samples( unique(ge$XenaHosts),  ge$XenaDatasets, limit = NULL)
  
  gene_info_km <- fetch_dense_values(host = unique(ge$XenaHosts), 
                                     dataset = ge$XenaDatasets,
                                     samples = samples,
                                     use_probeMap = F,
                                     check = FALSE,
                                     time_limit = 10000)
  
  gene_info_km <- as.data.frame((gene_info_km))
  
  
  #calculate gene sets
  
  rankData <- rankGenes(gene_info_km)
  
  scoredf_core <- simpleScore(rankData, upSet = core_genes , knownDirection = FALSE)
  scoredf_edge <- simpleScore(rankData, upSet = edge_genes , knownDirection = FALSE)
  
  #calculate gene sets
  names(scoredf_core) <- c("core","core_dispersion")
  names(scoredf_edge) <- c("edge","edge_dispersion")
  
  scores_df_singshot <- cbind(scoredf_core,scoredf_edge)
  
  scores_df_singshot <- cbind(sample = rownames(scores_df_singshot), data.frame(scores_df_singshot, row.names=NULL))
  
  merged_data = survival_data %>% 
    left_join(scores_df_singshot, by = "sample") %>% 
    dplyr::select(sample,names(geneset_of_interest), PFI.time, PFI)  %>% 
    left_join(clinical_data, by = "sample") %>% 
    dplyr::select(sample,names(geneset_of_interest), PFI.time, PFI,sample_type)  %>% 
    dplyr::rename(time = PFI.time, 
                  status = PFI)
  
  if (indication != "LAML") {
    merged_data = merged_data[merged_data$sample_type == "Primary Tumor",]
  }
  
  merged_data <- merged_data[complete.cases(merged_data), ]
  
  cut = surv_cutpoint(data = merged_data, time = "time", event = "status", variables = names(geneset_of_interest), minprop = 0.1)
  res.cat <- surv_categorize(cut)
  
  ##switch stuff here more
  for (gene in names(geneset_of_interest)) {
    
    res.cat$binary = res.cat[[gene]]
    
    row = df[df$gene == gene & df$indication == indication,]
    rownumber = row.names(df[df$gene == gene & df$indication == indication,])
    
    #try high ref
    res.cat$binary <- factor(res.cat$binary, levels = c("high","low"))
    
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


dir.create("heatmap_tcga")
setwd("heatmap_tcga")


png("survival_dotplot_PFI.png", units="in", width=3, height=5, res=300, type = "cairo")
ggplot(df, aes(x = geneset, y = indication, color = sig)) +
  geom_point(aes(size = hazard_ratio, fill = sig)) +
  scale_colour_manual(values = c("High" = "#d73027",
                                 "Low" = "#4575b4","ns" ="grey")) +
  labs(x="geneset", y="indication",
       title="Decreased PFI") + 
  theme_minimal()
dev.off()
setwd("..")
