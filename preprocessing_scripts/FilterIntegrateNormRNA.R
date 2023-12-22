# ---------------------------------------------------------------------------
# Filter, normalize & integrate/concatenate samples for sn/scRNA-seq data
# ---------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(Seurat)
library(reticulate)
library(sceasy)
library(stringr)

# Set up ---------------------------

args <- commandArgs(trailingOnly = TRUE)

# Define args
input_dir <- args[1]
strategy <- args[2]
print(paste0("Strategy: ", strategy))
sample_ids <- c(args[3:length(args)])
print(sample_ids)

# Import sample data
n_samples = length(sample_ids)
print(paste0("Number of samples: ", n_samples))

file_list <- c()

# Run filtering, normalization & merging ---------------------------

seurat_obj_list <- list()
for (j in 1:n_samples) {
  sample_j <- sample_ids[j]
  sce_j <- readRDS(paste0(input_dir, "/intermediate_files/qc/", sample_j, "/qc_matrix.rds"))
  matrix_j <- sce_j@assays@data$counts
  sample_j <- sample_ids[j]
  cell_metadata_j <- readRDS(paste0(input_dir, "/intermediate_files/qc/", sample_j, "/qc_cell_metadata.rds")) %>%
    mutate(Sample_ID = sample_j)
  if (strategy == "ground_truth") {
    no_ground_truth <- dplyr::filter(cell_metadata_j, is.na(Ground_truth))$CellID
    # Exclude cells based on ground truth availability
    matrix_j <- matrix_j[ , !(colnames(matrix_j) %in% no_ground_truth)]
    cell_metadata_j <- cell_metadata_j %>% dplyr::filter(!is.na(Ground_truth))
  } else {
    low_quality <- cell_metadata_j$CellID[cell_metadata_j[,"quality"] == "low quality"]
    doublets <- cell_metadata_j$CellID[cell_metadata_j[,"doublet_class"] == "Doublet"]
    # Exclude cells based on QC
    matrix_j <- matrix_j[ , !(colnames(matrix_j) %in% low_quality)]
    matrix_j <- matrix_j[ , !(colnames(matrix_j) %in% doublets)]
    cell_metadata_j <- cell_metadata_j %>% dplyr::filter(!(CellID %in% low_quality), !(CellID %in% doublets))
  }
  rownames(cell_metadata_j) <- cell_metadata_j$CellID
  seurat_obj_j <- CreateSeuratObject(matrix_j, meta.data = cell_metadata_j)
  seurat_obj_list <- append(seurat_obj_list, seurat_obj_j)
}

# Merge/normalize ---------------------------
if (n_samples > 1) {
  print("Merging..")
  
  # Simple merge ---------------------------
  seurat_obj_vector <- unlist(seurat_obj_list[2:length(seurat_obj_list)])
  seurat_obj <- merge(seurat_obj_list[[1]], y = seurat_obj_vector)
  
  # Metadata
  write.csv(seurat_obj@meta.data, paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata_RNA.csv"), row.names = FALSE)
  
  # Raw matrix
  saveRDS(seurat_obj@assays$RNA$counts, paste0(input_dir, "/intermediate_files/preprocessed/RNA_matrix_raw.rds"))
  write.csv(seurat_obj@assays$RNA$counts, paste0(input_dir, "/intermediate_files/preprocessed/RNA_matrix_raw.csv"))
  file_list <- c(file_list, paste0("RNA_matrix_raw"))
  
  # Normalized matrix
  seurat_norm <- seurat_obj %>%
    NormalizeData()
  norm_data <- seurat_norm@assays$RNA$data
  rownames(norm_data) <- rownames(seurat_norm)
  colnames(norm_data) <- colnames(seurat_norm)
  saveRDS(norm_data, paste0(input_dir, "/intermediate_files/preprocessed/RNA_matrix_norm_log.rds"))
  write.csv(norm_data, paste0(input_dir, "/intermediate_files/preprocessed/RNA_matrix_norm_log.csv"))
  file_list <- c(file_list, paste0("RNA_matrix_norm_log"))
  
  # Preview UMAP colored by sample
  seurat_norm %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30) %>%
    DimPlot(group.by = "Sample_ID", shuffle = TRUE) + NoLegend()
  ggsave(paste0(input_dir, "/quality_checks/", "sample_umap_norm_log.jpg"))
  
  # Preview UMAP colored by batch
  seurat_norm %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30) %>%
    DimPlot(group.by = "Batch", shuffle = TRUE)
  ggsave(paste0(input_dir, "/quality_checks/", "batch_umap_norm_log.jpg"))
  
} else {
  # 1 sample ---------------------------
  seurat_obj <- seurat_obj_list[[1]]
  
  # Metadata
  write.csv(seurat_obj@meta.data, paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata_RNA.csv"), row.names = FALSE)
  
  # Raw matrix
  saveRDS(seurat_obj@assays$RNA$counts, paste0(input_dir, "/intermediate_files/preprocessed/RNA_matrix_raw.rds"))
  write.csv(seurat_obj@assays$RNA$counts, paste0(input_dir, "/intermediate_files/preprocessed/RNA_matrix_raw.csv"))
  file_list <- c(file_list, paste0("RNA_matrix_raw"))
  
  # Normalized matrix
  seurat_obj <- seurat_obj %>%
    NormalizeData()
  norm_data <- seurat_obj@assays$RNA$data
  rownames(norm_data) <- rownames(seurat_obj)
  colnames(norm_data) <- colnames(seurat_obj)
  saveRDS(norm_data, paste0(input_dir, "/intermediate_files/preprocessed/RNA_matrix_norm_log.rds"))
  write.csv(norm_data, paste0(input_dir, "/intermediate_files/preprocessed/RNA_matrix_norm_log.csv"))
  file_list <- c(file_list, paste0("RNA_matrix_", assay, "_norm_log"))
  
  # Preview UMAP
  seurat_obj %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30) %>%
    DimPlot(group.by = "Sample_ID", shuffle = TRUE)
  ggsave(paste0(input_dir, "/quality_checks/", "sample_umap_norm_log.jpg"))
  
}

# Output file
write.csv(data.frame(File = file_list),
          file = paste0(input_dir, "/quality_checks/RNA_preprocessed_file_list.csv"))
