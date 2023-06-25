library("hardhat")
library("scPred")
library("Seurat")
library("magrittr")
library(tibble)

load("annotated_objects.Robj")

#reduce all objects
SCC_all_samples <- purrr::reduce(annotated_objects,merge)

#load in cancer cell annotations
load(file = "comb.Robj")
comb@meta.data$cluster_annotations = comb@active.ident
cancer_cell_meta <-comb@meta.data[,colnames(comb@meta.data) %in% c("new_annotations_trial","cluster_annotations")]

#transfer cancer cell state information, core, transitory, and edge
full_meta <- list(SCC_all_samples@meta.data, cancer_cell_meta) %>% 
  purrr::map(~ .x %>% 
        as.data.frame %>%
        rownames_to_column('rn')) %>% 
        purrr::reduce(left_join, by = 'rn') %>%
        column_to_rownames('rn')

full_meta$cluster_annotations <- as.character(full_meta$cluster_annotations)
#label noncancer cells
full_meta$cluster_annotations[is.na(full_meta$cluster_annotations)] <- "noncancer"

full_meta$cluster_annotations <- factor(full_meta$cluster_annotations,levels = unique(full_meta$cluster_annotations))

#transfer annotations to object
SCC_all_samples@meta.data <- full_meta

#generate reference
reference <- SCC_all_samples
DefaultAssay(reference) <- "Spatial"
reference <- reference %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunHarmony("sample_id")  %>% 
  RunPCA(pc.genes = reference@var.genes, npcs = 20, verbose = FALSE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

reference <- getFeatureSpace(reference, "cluster_annotations")
rm(list=setdiff(ls(), "reference"))
#train model on avNNet, svmRadial, and nb
reference_avNNet <- trainModel(reference,allowParallel = T,preProcess = c("center", "scale","YeoJohnson"),seed = 21,model = "avNNet",number = 10)
reference_svmRadial <- trainModel(reference,allowParallel = T,preProcess = c("center", "scale","YeoJohnson"),seed = 21,model = "svmRadial",number = 10)
reference_nb <- trainModel(reference,allowParallel = T,preProcess = c("center", "scale","YeoJohnson"),seed = 21,model = "nb",number = 10)

save(reference_avNNet, file = "reference_avNNet.Robj")
save(reference_svmRadial, file = "reference_svmRadial.Robj")
save(reference_nb, file = "reference_nb.Robj")
load(file = "reference_avNNet.Robj")
load(file = "reference_svmRadial.Robj")
load(file = "reference_nb.Robj")

get_scpred(reference_avNNet)
get_scpred(reference_svmRadial)
get_scpred(reference_nb)

#based on the best models for each subtype, we should do:
#core: svm
#edge: svm
#transitory: avNNet
#noncancer: svm

#reclassify svmRadial with avNNet for transitory 
reference <- trainModel(reference_svmRadial,allowParallel = T,preProcess = c("center", "scale","YeoJohnson"),seed = 21,model = "avNNet",reclassify = c("transitory"),number = 10)
save(reference, file = "reference.Robj")
rm(list=setdiff(ls(), "reference"))
### On server

library(Seurat)
library("scPred")
library("magrittr")
load(file = "reference.Robj")

get_scpred(reference)
model <- get_scpred(reference)
save(model,file = "model.Robj")
#plot model statistics
png(file = "plot_probabilities_reference_1.png", 
     width = 6, height = 6, units = "in", res = 300,type = "cairo")
plot_probabilities(reference)
dev.off()

#load all cancers and predict cancer cell states on them

## CHC
load(file = "objects/CHC.Robj")
query = CHC
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/CHC.Robj")

##CHC1_T 
load(file = "objects/CHC1_T.Robj")
query = CHC1_T
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/CHC1_T.Robj")


#crc 
load(file = "objects/crc.Robj")
query = crc
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/crc.Robj")

#gbm
load(file = "objects/gbm.Robj")
query = gbm
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/gbm.Robj")


#HBC_IDC
load(file = "objects/HBC_IDC.Robj")
query = HBC_IDC
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/HBC_IDC.Robj")

#HBC_ILC
load(file = "objects/HBC_ILC.Robj")
query = HBC_ILC
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/HBC_ILC.Robj")


#HCC1_T
load(file = "objects/HCC1_T.Robj")
query = HCC1_T
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/HCC1_T.Robj")

#HCC1
load(file = "objects/HCC1.Robj")
query = HCC1
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/HCC1.Robj")

#HCC2
load(file = "objects/HCC2_T.Robj")
query = HCC2_T
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/HCC2_T.Robj")

#ICC
load(file = "objects/ICC.Robj")
query = ICC
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/ICC.Robj")

#OV
load(file = "objects/OV.Robj")
query = OV
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/OV.Robj")

#p4_scc
load(file = "objects/p4_scc.Robj")
query = p4_scc
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/p4_cscc.Robj")

#p6_scc
load(file = "objects/p6_scc.Robj")
query = p6_scc
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/p6_cscc.Robj")

rm(list=setdiff(ls(), "reference"))

load(file = "objects/PDAC_A.Robj")
query = PDAC_A
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/PDAC_A.Robj")

load(file = "objects/PDAC_B.Robj")
query = PDAC_B
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/PDAC_B.Robj")

load(file = "objects/PDAC_C.Robj")
query = PDAC_C
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/PDAC_C.Robj")

load(file = "objects/cscc1.Robj")
query = cscc1
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/cscc1.Robj")

load(file = "objects/cscc2.Robj")
query = cscc2
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/cscc2.Robj")

load(file = "objects/cscc3.Robj")
query = cscc3
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/cscc3.Robj")

load(file = "objects/cscc4.Robj")
query = cscc4
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/cscc4.Robj")

rm(list=setdiff(ls(), "reference"))

load(file = "objects/cervical_scc.Robj")
query = cervical_scc
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/cervical_scc.Robj")

load(file = "objects/intestinal.Robj")
query = intestinal
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/intestinals.Robj")

load(file = "objects/prostate.Robj")
query = prostate
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/prostate.Robj")

load(file = "objects/prostate_acinar.Robj")
query = prostate_acinar
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/prostate_acinar.Robj")

load(file = "objects/lung_scc.Robj")
query = lung_scc
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/lung_scc.Robj")

load(file = "objects/melanoma.Robj")
query = melanoma
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/melanoma.Robj")

rm(list=setdiff(ls(), "reference"))

load(file = "objects/s1_paediatric_medulloblastoma.Robj")
query = s1_paediatric_medulloblastoma
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/s1_paediatric_medulloblastoma.Robj")

load(file = "objects/s2_paediatric_medulloblastoma.Robj")
query = s2_paediatric_medulloblastoma
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/s2_paediatric_medulloblastoma.Robj")

load(file = "objects/s3_paediatric_CNS_embryonal.Robj")
query = s3_paediatric_CNS_embryonal
query <- NormalizeData(query)
query <- scPredict(query, reference, threshold = 0.2)
save(query, file = "processed_data/s3_paediatric_CNS_embryonal.Robj")


