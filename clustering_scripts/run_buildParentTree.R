# ---------------------------------------------------------------------------
# Run CHOIR buildParentTree on Siletti et al. 2023 human brain atlas dataset
# ---------------------------------------------------------------------------
library(Seurat)
library(stringr)
library(CHOIR)
library(BPCells)
library(SeuratObject)
library(dplyr)

start_time <- Sys.time()
input_dir <- "/gladstone/mucke/sequencing/Cathrine/cluster_benchmarking/datasets/Siletti_2023"
input_data1 <- "bpcells_raw"

# Import data
cell_metadata <- read.csv(paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata.csv"))
options(Seurat.object.assay.version = 'v5')
modality1 <- open_matrix_dir(dir = paste0(input_dir, "/intermediate_files/preprocessed/", input_data1))

# Create Seurat object
seurat_object <- CreateSeuratObject(counts = modality1, assay = "RNA")
# Normalize
seurat_object <- NormalizeData(seurat_object)
  
# Cell metadata
rownames(cell_metadata) <- cell_metadata$CellID
cell_metadata <- cell_metadata[colnames(seurat_object),]
seurat_object@meta.data$Batch <- as.character(cell_metadata$Batch)
batch_labels <- "Batch"

# Build parent tree
seurat_object <- buildParentTree(seurat_object, batch_correction_method = "Harmony", batch_labels = "Batch", n_cores = 16)
end_time <- Sys.time()
diff_time <- difftime(end_time, start_time, units = "secs")

saveRDS(seurat_object@misc$CHOIR, paste0(input_dir, "/intermediate_files/preprocessed/CHOIR_object_parent_clusters.rds"))
saveRDS(seurat_object, paste0(input_dir, "/intermediate_files/preprocessed/object_parent_clusters.rds"))
saveRDS(diff_time, paste0(input_dir, "/buildParentTree_difftime.rds"))
