module load singularity/3.8.1

#singularity pull docker://pkharchenkolab/numbat-rbase:latest

salloc --time 4:59:00 --ntasks=40 --mem-per-cpu=40G

singularity exec --home /work/bose_lab/Rohit/numbat numbat-rbase_latest.sif Rscript /numbat/inst/bin/pileup_and_phase.R \
--label "sample_1" \
--samples "sample_1" \
--bams sample_data/sample_1/possorted_genome_bam.bam \
--barcodes sample_data/sample_1/barcodes.tsv \
--outdir "pap/sample_1" \
--gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
--snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
--paneldir /data/1000G_hg38 \
--ncores 30

singularity exec --home /work/bose_lab/Rohit/numbat numbat-rbase_latest.sif Rscript /numbat/inst/bin/pileup_and_phase.R \
--label "sample_2" \
--samples "sample_2" \
--bams sample_data/sample_2/possorted_genome_bam.bam \
--barcodes sample_data/sample_2/barcodes.tsv \
--outdir "pap/sample_2" \
--gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
--snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
--paneldir /data/1000G_hg38 \
--ncores 30

singularity exec --home /work/bose_lab/Rohit/numbat numbat-rbase_latest.sif Rscript /numbat/inst/bin/pileup_and_phase.R \
--label "sample_3" \
--samples "sample_3" \
--bams sample_data/sample_3/possorted_genome_bam.bam \
--barcodes sample_data/sample_3/barcodes.tsv \
--outdir "pap/sample_3" \
--gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
--snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
--paneldir /data/1000G_hg38 \
--ncores 30


singularity exec --home /work/bose_lab/Rohit/numbat numbat-rbase_latest.sif Rscript /numbat/inst/bin/pileup_and_phase.R \
--label "sample_4" \
--samples "sample_4" \
--bams sample_data/sample_4/possorted_genome_bam.bam \
--barcodes sample_data/sample_4/barcodes.tsv \
--outdir "pap/sample_4" \
--gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
--snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
--paneldir /data/1000G_hg38 \
--ncores 30


singularity exec --home /work/bose_lab/Rohit/numbat numbat-rbase_latest.sif Rscript /numbat/inst/bin/pileup_and_phase.R \
--label "sample_5" \
--samples "sample_5" \
--bams sample_data/sample_5/possorted_genome_bam.bam \
--barcodes sample_data/sample_5/barcodes.tsv \
--outdir "pap/sample_5" \
--gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
--snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
--paneldir /data/1000G_hg38 \
--ncores 30


singularity exec --home /work/bose_lab/Rohit/numbat numbat-rbase_latest.sif Rscript /numbat/inst/bin/pileup_and_phase.R \
--label "sample_6" \
--samples "sample_6" \
--bams sample_data/sample_6/possorted_genome_bam.bam \
--barcodes sample_data/sample_6/barcodes.tsv \
--outdir "pap/sample_6" \
--gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
--snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
--paneldir /data/1000G_hg38 \
--ncores 30


singularity exec --home /work/bose_lab/Rohit/numbat numbat-rbase_latest.sif Rscript /numbat/inst/bin/pileup_and_phase.R \
--label "sample_7" \
--samples "sample_7" \
--bams sample_data/sample_7/possorted_genome_bam.bam \
--barcodes sample_data/sample_7/barcodes.tsv \
--outdir "pap/sample_7" \
--gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
--snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
--paneldir /data/1000G_hg38 \
--ncores 30

singularity exec --home /work/bose_lab/Rohit/numbat numbat-rbase_latest.sif Rscript /numbat/inst/bin/pileup_and_phase.R \
--label "sample_8" \
--samples "sample_8" \
--bams sample_data/sample_8/possorted_genome_bam.bam \
--barcodes sample_data/sample_8/barcodes.tsv \
--outdir "pap/sample_8" \
--gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
--snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
--paneldir /data/1000G_hg38 \
--ncores 30


singularity exec --home /work/bose_lab/Rohit/numbat numbat-rbase_latest.sif Rscript /numbat/inst/bin/pileup_and_phase.R \
--label "sample_9" \
--samples "sample_9" \
--bams sample_data/sample_9/possorted_genome_bam.bam \
--barcodes sample_data/sample_9/barcodes.tsv \
--outdir "pap/sample_9" \
--gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
--snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
--paneldir /data/1000G_hg38 \
--ncores 30


singularity exec --home /work/bose_lab/Rohit/numbat numbat-rbase_latest.sif Rscript /numbat/inst/bin/pileup_and_phase.R \
--label "sample_10" \
--samples "sample_10" \
--bams sample_data/sample_10/possorted_genome_bam.bam \
--barcodes sample_data/sample_10/barcodes.tsv \
--outdir "pap/sample_10" \
--gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
--snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
--paneldir /data/1000G_hg38 \
--ncores 30


singularity exec --home /work/bose_lab/Rohit/numbat numbat-rbase_latest.sif Rscript /numbat/inst/bin/pileup_and_phase.R \
--label "sample_11" \
--samples "sample_11" \
--bams sample_data/sample_11/possorted_genome_bam.bam \
--barcodes sample_data/sample_11/barcodes.tsv \
--outdir "pap/sample_11" \
--gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
--snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
--paneldir /data/1000G_hg38 \
--ncores 30


singularity exec --home /work/bose_lab/Rohit/numbat numbat-rbase_latest.sif Rscript /numbat/inst/bin/pileup_and_phase.R \
--label "sample_12" \
--samples "sample_12" \
--bams sample_data/sample_12/possorted_genome_bam.bam \
--barcodes sample_data/sample_12/barcodes.tsv \
--outdir "pap/sample_12" \
--gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
--snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
--paneldir /data/1000G_hg38 \
--ncores 30



#prep h5 files for each sample
library(Seurat)
mtx1 <- Read10X_h5("mtx_files/s1_filtered_feature_bc_matrix.h5")
saveRDS(mtx1, file ="mtx_processed/sample_1/mtx.RDS")

mtx2 <- Read10X_h5("mtx_files/s2_filtered_feature_bc_matrix.h5")
saveRDS(mtx2, file ="mtx_processed/sample_2/mtx.RDS")

mtx3 <- Read10X_h5("mtx_files/s3_filtered_feature_bc_matrix.h5")
saveRDS(mtx3, file ="mtx_processed/sample_3/mtx.RDS")

mtx4 <- Read10X_h5("mtx_files/s4_filtered_feature_bc_matrix.h5")
saveRDS(mtx4, file ="mtx_processed/sample_4/mtx.RDS")

mtx5 <- Read10X_h5("mtx_files/s5_filtered_feature_bc_matrix.h5")
saveRDS(mtx5, file ="mtx_processed/sample_5/mtx.RDS")

mtx6 <- Read10X_h5("mtx_files/s6_filtered_feature_bc_matrix.h5")
saveRDS(mtx6, file ="mtx_processed/sample_6/mtx.RDS")

mtx7 <- Read10X_h5("mtx_files/s7_filtered_feature_bc_matrix.h5")
saveRDS(mtx7, file ="mtx_processed/sample_7/mtx.RDS")

mtx8 <- Read10X_h5("mtx_files/s8_filtered_feature_bc_matrix.h5")
saveRDS(mtx8, file ="mtx_processed/sample_8/mtx.RDS")

mtx9 <- Read10X_h5("mtx_files/s9_filtered_feature_bc_matrix.h5")
saveRDS(mtx9, file ="mtx_processed/sample_9/mtx.RDS")

mtx10 <- Read10X_h5("mtx_files/s10_filtered_feature_bc_matrix.h5")
saveRDS(mtx10, file ="mtx_processed/sample_10/mtx.RDS")

mtx11 <- Read10X_h5("mtx_files/s11_filtered_feature_bc_matrix.h5")
saveRDS(mtx11, file ="mtx_processed/sample_11/mtx.RDS")

mtx12 <- Read10X_h5("mtx_files/s12_filtered_feature_bc_matrix.h5")
saveRDS(mtx12, file ="mtx_processed/sample_12/mtx.RDS")

#gunzip allele counts file
gunzip pap/sample_1/sample_1_allele_counts.tsv.gz 
gunzip pap/sample_2/sample_2_allele_counts.tsv.gz 
gunzip pap/sample_3/sample_3_allele_counts.tsv.gz 
gunzip pap/sample_4/sample_4_allele_counts.tsv.gz 
gunzip pap/sample_5/sample_5_allele_counts.tsv.gz 
gunzip pap/sample_6/sample_6_allele_counts.tsv.gz 
gunzip pap/sample_7/sample_7_allele_counts.tsv.gz 
gunzip pap/sample_8/sample_8_allele_counts.tsv.gz 
gunzip pap/sample_9/sample_9_allele_counts.tsv.gz 
gunzip pap/sample_10/sample_10_allele_counts.tsv.gz 
gunzip pap/sample_11/sample_11_allele_counts.tsv.gz 
gunzip pap/sample_12/sample_12_allele_counts.tsv.gz 


salloc --time 4:59:00 --ntasks=40 --mem-per-cpu=40G
#run once for each sample
singularity exec --home /work/bose_lab/Rohit/numbat numbat-rbase_latest.sif R

library(numbat)

count_mat <- readRDS("sample_data/sample_1/mtx.RDS")
allele_df <- read.delim("pap/sample_1/sample_1_allele_counts.tsv") 
out = run_numbat(
  count_mat, # gene x cell integer UMI count matrix 
  ref_hca, # reference expression profile, a gene x cell type normalized expression level matrix
  allele_df, # allele dataframe generated by pileup_and_phase script
  genome = "hg38",
  t = 1e-5,
  ncores = 4,
  plot = TRUE,
  out_dir = 'outs/sample_1',
  max_entropy = 0.8
)

q() 

n 

singularity exec --home /work/bose_lab/Rohit/numbat numbat-rbase_latest.sif R

library(numbat)

count_mat <- readRDS("sample_data/sample_2/mtx.RDS")
allele_df <- read.delim("pap/sample_2/sample_2_allele_counts.tsv") 
out = run_numbat(
  count_mat, # gene x cell integer UMI count matrix 
  ref_hca, # reference expression profile, a gene x cell type normalized expression level matrix
  allele_df, # allele dataframe generated by pileup_and_phase script
  genome = "hg38",
  t = 1e-5,
  ncores = 4,
  plot = TRUE,
  out_dir = 'outs/sample_2',
  max_entropy = 0.8
)

q() 

n 

singularity exec --home /work/bose_lab/Rohit/numbat numbat-rbase_latest.sif R

library(numbat)

count_mat <- readRDS("sample_data/sample_3/mtx.RDS")
allele_df <- read.delim("pap/sample_3/sample_3_allele_counts.tsv") 
out = run_numbat(
  count_mat, # gene x cell integer UMI count matrix 
  ref_hca, # reference expression profile, a gene x cell type normalized expression level matrix
  allele_df, # allele dataframe generated by pileup_and_phase script
  genome = "hg38",
  t = 1e-5,
  ncores = 4,
  plot = TRUE,
  out_dir = 'outs/sample_3',
  max_entropy = 0.8
)

q() 

n 

singularity exec --home /work/bose_lab/Rohit/numbat numbat-rbase_latest.sif R

library(numbat)

count_mat <- readRDS("sample_data/sample_4/mtx.RDS")
allele_df <- read.delim("pap/sample_4/sample_4_allele_counts.tsv") 
out = run_numbat(
  count_mat, # gene x cell integer UMI count matrix 
  ref_hca, # reference expression profile, a gene x cell type normalized expression level matrix
  allele_df, # allele dataframe generated by pileup_and_phase script
  genome = "hg38",
  t = 1e-5,
  ncores = 4,
  plot = TRUE,
  out_dir = 'outs/sample_4',
  max_entropy = 0.8
)

q() 

n 

singularity exec --home /work/bose_lab/Rohit/numbat numbat-rbase_latest.sif R

library(numbat)

count_mat <- readRDS("sample_data/sample_5/mtx.RDS")
allele_df <- read.delim("pap/sample_5/sample_5_allele_counts.tsv") 
out = run_numbat(
  count_mat, # gene x cell integer UMI count matrix 
  ref_hca, # reference expression profile, a gene x cell type normalized expression level matrix
  allele_df, # allele dataframe generated by pileup_and_phase script
  genome = "hg38",
  t = 1e-5,
  ncores = 4,
  plot = TRUE,
  out_dir = 'outs/sample_5',
  max_entropy = 0.8
)

q() 

n 

singularity exec --home /work/bose_lab/Rohit/numbat numbat-rbase_latest.sif R

library(numbat)

count_mat <- readRDS("sample_data/sample_6/mtx.RDS")
allele_df <- read.delim("pap/sample_6/sample_6_allele_counts.tsv") 
out = run_numbat(
  count_mat, # gene x cell integer UMI count matrix 
  ref_hca, # reference expression profile, a gene x cell type normalized expression level matrix
  allele_df, # allele dataframe generated by pileup_and_phase script
  genome = "hg38",
  t = 1e-5,
  ncores = 4,
  plot = TRUE,
  out_dir = 'outs/sample_6',
  max_entropy = 0.8
)

q() 

n 

singularity exec --home /work/bose_lab/Rohit/numbat numbat-rbase_latest.sif R

library(numbat)

count_mat <- readRDS("sample_data/sample_7/mtx.RDS")
allele_df <- read.delim("pap/sample_7/sample_7_allele_counts.tsv") 
out = run_numbat(
  count_mat, # gene x cell integer UMI count matrix 
  ref_hca, # reference expression profile, a gene x cell type normalized expression level matrix
  allele_df, # allele dataframe generated by pileup_and_phase script
  genome = "hg38",
  t = 1e-5,
  ncores = 4,
  plot = TRUE,
  out_dir = 'outs/sample_7',
  max_entropy = 0.8
)

q() 

n 


singularity exec --home /work/bose_lab/Rohit/numbat numbat-rbase_latest.sif R

library(numbat)

count_mat <- readRDS("sample_data/sample_8/mtx.RDS")
allele_df <- read.delim("pap/sample_8/sample_8_allele_counts.tsv") 
out = run_numbat(
  count_mat, # gene x cell integer UMI count matrix 
  ref_hca, # reference expression profile, a gene x cell type normalized expression level matrix
  allele_df, # allele dataframe generated by pileup_and_phase script
  genome = "hg38",
  t = 1e-5,
  ncores = 4,
  plot = TRUE,
  out_dir = 'outs/sample_8',
  max_entropy = 0.8
)

q() 

n 

singularity exec --home /work/bose_lab/Rohit/numbat numbat-rbase_latest.sif R

library(numbat)

count_mat <- readRDS("sample_data/sample_9/mtx.RDS")
allele_df <- read.delim("pap/sample_9/sample_9_allele_counts.tsv") 
out = run_numbat(
  count_mat, # gene x cell integer UMI count matrix 
  ref_hca, # reference expression profile, a gene x cell type normalized expression level matrix
  allele_df, # allele dataframe generated by pileup_and_phase script
  genome = "hg38",
  t = 1e-5,
  ncores = 4,
  plot = TRUE,
  out_dir = 'outs/sample_9',
  max_entropy = 0.8
)

q() 

n 

singularity exec --home /work/bose_lab/Rohit/numbat numbat-rbase_latest.sif R

library(numbat)

count_mat <- readRDS("sample_data/sample_10/mtx.RDS")
allele_df <- read.delim("pap/sample_10/sample_10_allele_counts.tsv") 
out = run_numbat(
  count_mat, # gene x cell integer UMI count matrix 
  ref_hca, # reference expression profile, a gene x cell type normalized expression level matrix
  allele_df, # allele dataframe generated by pileup_and_phase script
  genome = "hg38",
  t = 1e-5,
  ncores = 4,
  plot = TRUE,
  out_dir = 'outs/sample_10',
  max_entropy = 0.8
)

q() 

n 

singularity exec --home /work/bose_lab/Rohit/numbat numbat-rbase_latest.sif R

library(numbat)

count_mat <- readRDS("sample_data/sample_11/mtx.RDS")
allele_df <- read.delim("pap/sample_11/sample_11_allele_counts.tsv") 
out = run_numbat(
  count_mat, # gene x cell integer UMI count matrix 
  ref_hca, # reference expression profile, a gene x cell type normalized expression level matrix
  allele_df, # allele dataframe generated by pileup_and_phase script
  genome = "hg38",
  t = 1e-5,
  ncores = 4,
  plot = TRUE,
  out_dir = 'outs/sample_11',
  max_entropy = 0.8
)

q() 

n 

singularity exec --home /work/bose_lab/Rohit/numbat numbat-rbase_latest.sif R

library(numbat)

count_mat <- readRDS("sample_data/sample_12/mtx.RDS")
allele_df <- read.delim("pap/sample_12/sample_12_allele_counts.tsv") 
out = run_numbat(
  count_mat, # gene x cell integer UMI count matrix 
  ref_hca, # reference expression profile, a gene x cell type normalized expression level matrix
  allele_df, # allele dataframe generated by pileup_and_phase script
  genome = "hg38",
  t = 1e-5,
  ncores = 4,
  plot = TRUE,
  out_dir = 'outs/sample_12',
  max_entropy = 0.8
)

q() 

n 

