# ---------------------------------------------------------------------------
# This script runs MAGIC imputation for selected genes
# ---------------------------------------------------------------------------
library(dplyr)         
library(Seurat)       
library(Rmagic)

input_dir <- "/gladstone/corces/lab/users/cpetersen/cluster_benchmarking/Wang_2022_multiome"
count_matrix <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/RNA_matrix_raw.rds"))

seurat_object <- CreateSeuratObject(counts = count_matrix)
cell_metadata <- read.csv(paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata.csv"))
seurat_object <- subset(seurat_object, cells = cell_metadata$CellID)
seurat_object <- NormalizeData(seurat_object)

seurat_magic <- magic(seurat_object, knn=5, genes = c("ENSG00000182674", "ENSG00000049759",
                                                      "ENSG00000144278"))
saveRDS(seurat_magic, paste0(input_dir, "/seurat_magic.rds"))