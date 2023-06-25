library(PharmacoGx)
'%!in%' <- function(x,y)!('%in%'(x,y))

CCLE <- readRDS("files/CCLE.rds")
CCLE <- updateObject(CCLE)

CTRPv2 <- readRDS("files/PSet_CTRPv2.rds")
CTRPv2 <- updateObject(CTRPv2)

PRISM <- readRDS("files/PSet_PRISM.rds")
PRISM <- updateObject(PRISM)

GDSC2 <- readRDS("files/PSet_GDSC2020.rds")
GDSC2 <- updateObject(GDSC2)

gCSI <- readRDS("files/PSet_gCSI2019.rds")
gCSI <- updateObject(gCSI)

##ccle
CCLE_cell_line_id <- CCLE@sample[CCLE@sample$tissueid %in% c("Head and Neck") ,]$sampleid

CTRPv2_cell_line_id <- CTRPv2@sample[CTRPv2@sample$tissueid %in% c("Head and Neck") ,]$sampleid

PRISM_cell_line_id <- PRISM@sample[PRISM@sample$tissueid %in% c("Head and Neck") ,]$sampleid

gCSI_cell_line_id <- gCSI@sample[gCSI@sample$tissueid %in% c("Head and Neck") ,]$sampleid

GDSC2_cell_line_id <- GDSC2@sample[GDSC2@sample$tissueid %in% c("Head and Neck") ,]$sampleid


#64 HNSC cell lines currently
all_cell_lines <- union(union(union(union(CCLE_cell_line_id,gCSI_cell_line_id),GDSC2_cell_line_id),CTRPv2_cell_line_id),PRISM_cell_line_id)

HNSC_cl_class <- read.csv(file = "HNSC_cl_class.csv")
HNSC_cl_class_pos <- HNSC_cl_class[HNSC_cl_class$HPV.Status == "HPV-negative",]$HNSC.Cell.Line
  
all_cell_lines = HNSC_cl_class_pos


i = 1

treat_response_studies <- c(CTRPv2,PRISM,GDSC2,CCLE,gCSI)

CCLE_drugs<- unique(CCLE@treatmentResponse$info$treatmentid)
CTRPv2_drugs<- unique(CTRPv2@treatmentResponse$info$treatmentid)
PRISM_drugs<- unique(PRISM@treatmentResponse$info$treatmentid)
gCSI_drugs <- unique(gCSI@treatmentResponse$info$treatmentid)
GDSC2_drugs <- unique(GDSC2@treatmentResponse$info$treatmentid)

all_drugs <- union(union(union(union(gCSI_drugs,GDSC2_drugs),CTRPv2_drugs),gCSI_drugs),PRISM_drugs)

auc_df = data.frame(matrix(nrow = length(all_drugs), ncol = length(all_cell_lines))) 
colnames(auc_df) <- all_cell_lines
rownames(auc_df) <- all_drugs

for (treat_response in treat_response_studies){
  
  response_to_treatment <- PharmacoGx::summarizeSensitivityProfiles(treat_response,
                                                                    sensitivity.measure='aac_recomputed',
                                                                    summary.stat="mean",
                                                                    verbose=FALSE)
  response_to_treatment <- as.data.frame(response_to_treatment)
  response_to_treatment <- response_to_treatment[,colnames(response_to_treatment) %in% all_cell_lines]
  response_to_treatment <- response_to_treatment[,colSums(is.na(response_to_treatment))<nrow(response_to_treatment)]
  for (rowenum in 1:nrow(response_to_treatment)){
    row = response_to_treatment[rowenum,]
    for (colenum in 1:ncol(row)) {
      value = response_to_treatment[rowenum,colenum]
      if (!is.na(value)){
        drug_name <- rownames(response_to_treatment)[rowenum]
        sample_name <- colnames(response_to_treatment)[colenum]
        if (is.na(auc_df[drug_name,sample_name])) {
          auc_df[drug_name,sample_name] = value
        } else {
          auc_df[drug_name,sample_name]  = paste0(auc_df[drug_name,sample_name], "&",value)
        }
      }
    }
  }
}

mean_applied = function(value) {
  if(!is.na(value)) {
  if(stringr:::str_detect(value, "&")){
    avg <- mean(as.numeric(unlist((stringr::str_split(value, "&")[[1]]))))
    return(avg)
  }
    return(value)
  }
    return(value)
}

auc_df <-apply(auc_df, FUN = mean_applied, MARGIN = c(1,2)) 

#rank drugs, only include drugs in more than half of cell lines
library(dplyr)
auc_df_reduced <- auc_df[,colSums(is.na(auc_df))<nrow(auc_df)]
auc_df_reduced <- auc_df_reduced[rowSums(is.na(auc_df_reduced))<ncol(auc_df_reduced)/2,]
auc_df_reduced <- as.data.frame(apply(auc_df_reduced,FUN = as.numeric, MARGIN = c(1,2)))
auc_df_reduced$mean = 0

#trimmed mean calculation
for (row in 1:nrow(auc_df_reduced)) {
  mean <- base::mean(t(auc_df_reduced[row,])[,1], trim=0.1, na.rm=TRUE)
  auc_df_reduced[row,]$mean <- mean
}

auc_df_reduced <- auc_df_reduced %>% mutate(rank=dense_rank(mean))

write.csv(auc_df_reduced, file = "HNSC_hpvneg_drug_AAC.csv")

#
drug_ic50summarized <- read.csv("HNSC_hpvneg_drug_AAC.csv")
drug_names <- drug_ic50summarized$X

library(rDGIdb)
library(httr)
library(jsonlite)

body <- list(drugs = paste(drug_names, collapse = ","))
            
url <- "https://dgidb.org/api/v2/interactions.json"
body <- body[!sapply(body, is.null)]
httr::verbose()
postRequest <- POST(url = url, body = body, encode = 'multipart')

text <- content(postRequest, as = "text")
result <- fromJSON(text, simplifyVector = TRUE)

#create a dataframe of #drugname, #IC50 value, #genes up- or down-regulated

dgidb_ic50 <- data.frame(drugs = drug_ic50summarized$X, IC50_value <- drug_ic50summarized$mean, dgidb_drug_name = 0 ,
                         genes_up = "",
                         genes_down = "")

names(dgidb_ic50) <- c("drugs","AUC_mean","dgidb_drug_name","genes_up","genes_down")
dgidb_ic50$genes_up <- as.list(dgidb_ic50$genes_up)
dgidb_ic50$genes_down <- as.list(dgidb_ic50$genes_down)

matched_terms <- result$matchedTerms
matched_terms <- as.data.frame(matched_terms)
typeof(matched_terms)

#upregulated gene terms
is_enriched <- c("activator", "agonist","inducer","partial agonist","positive modulator","potentiator","stimulator")
  
#downregulated gene terms
is_depeleted <- c("inhibitor", "antagonist","partial antagonist","blocker","inverse agonist","negative modulator","suppressor")
  
truefunc <- function(value) {
  if(length(value > 0)) {
    return(all(value))
  } else {
    return(FALSE)
  }
}


for (row in rownames(matched_terms)) {
  row_vals <- matched_terms[row,]
  
  dgidb_name <- row_vals$drugName
  interactions <- as.data.frame(row_vals$interactions)
  
  interactions_deplete <- lapply(interactions$interactionTypes, function(x) any(x %in% is_depeleted))
  interactions_deplete <- lapply(interactions_deplete, truefunc)
  genes_deplete <- interactions[unlist(interactions_deplete),]$geneName
  
  interactions_enrich <- lapply(interactions$interactionTypes, function(x) any(x %in% is_enriched))
  interactions_enrich <- lapply(interactions_enrich, truefunc)
  genes_enrich <- interactions[unlist(interactions_enrich),]$geneName

  dgidb_ic50[toupper(dgidb_ic50$drugs) == row_vals$searchTerm,]$dgidb_drug_name <- dgidb_name
  if (length(genes_enrich) != 0 ) {
    dgidb_ic50[toupper(dgidb_ic50$drugs) == row_vals$searchTerm,]$genes_up<- purrr::reduce(genes_enrich, paste)
  }
  if (length(genes_deplete) != 0 ) {
    dgidb_ic50[toupper(dgidb_ic50$drugs) == row_vals$searchTerm,]$genes_down <- purrr::reduce(genes_deplete, paste)
  }
}

dgidb_AUC_values <- dgidb_ic50

dgidb_AUC_values$genes_up <- lapply(dgidb_AUC_values$genes_up, plyr::ldply)

dgidb_AUC_values$genes_up <- lapply(dgidb_AUC_values$genes_up, `[[`, 1)
  
dgidb_AUC_values$genes_down <-  lapply(dgidb_AUC_values$genes_down, plyr::ldply)

dgidb_AUC_values$genes_down <- lapply(dgidb_AUC_values$genes_down, `[[`, 1)

dgidb_AUC_values$genes_up <- as.character(dgidb_AUC_values$genes_up)

dgidb_AUC_values$genes_down <- as.character(dgidb_AUC_values$genes_down)


write.csv(dgidb_AUC_values, file = "dgidb_AUC_values.csv")

