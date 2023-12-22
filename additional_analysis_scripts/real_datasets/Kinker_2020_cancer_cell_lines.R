# ---------------------------------------------------------------------------
# This script contains the benchmarking analysis of the Kinker et al. 2020
# scRNA-seq cancer cell line dataset
# ---------------------------------------------------------------------------
library(BPCells)
library(Seurat)
library(CHOIR)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)

# Set color palette
color_palette <- CHOIR::CHOIRpalette(100)
alt_color_palette <- CHOIR::CHOIRpalette(500)

# Import data
input_dir <- "/gladstone/corces/lab/users/cpetersen/cluster_benchmarking/Kinker_2020_all"
count_matrix <- open_matrix_dir(dir = paste0(input_dir, "/bpcells"))

# Create Seurat Object
options(Seurat.object.assay.version = 'v5')
seurat_object <- CreateSeuratObject(counts = count_matrix)

# Add metadata
cell_metadata <- read.csv(paste0(input_dir, "/cell_metadata.csv"))
identical(colnames(seurat_object), cell_metadata$CellID)
seurat_object@meta.data <- cbind(seurat_object@meta.data, cell_metadata)

# Normalize & find variable features
seurat_object <- NormalizeData(seurat_object)

# Import CHOIR UMAP coordinates
choir_data <- readRDS(paste0(input_dir, "/CHOIR_output_1.rds"))
umap_coords <- RunUMAP(choir_data$reduction$P0_reduction)
# Add UMAP coordinates to Seurat object
seurat_object[["P0_umap"]] <- CreateDimReducObject(embeddings = umap_coords@cell.embeddings,
                                                   assay = "RNA")

# Import clusters
# CHOIR
CHOIR_clusters <- read.csv(paste0(input_dir, 
                                  "/compiled_clusters_CHOIR_parameter_list_RNA_real.csv"))
rownames(CHOIR_clusters) <- CHOIR_clusters$CellID
colnames(CHOIR_clusters) <- sub("Parameters", "CHOIR", colnames(CHOIR_clusters))
CHOIR_clusters <- CHOIR_clusters[seurat_object@meta.data$CellID,]
# Other methods
Other_methods_clusters <- read.csv(paste0(input_dir, 
                                          "/compiled_clusters_cluster_parameter_list_RNA_real.csv"))
rownames(Other_methods_clusters) <- Other_methods_clusters$CellID
Other_methods_clusters <- Other_methods_clusters[seurat_object@meta.data$CellID,]

# Add clustering results to Seurat object
seurat_object@meta.data <- cbind(seurat_object@meta.data, CHOIR_clusters[,-1])
seurat_object@meta.data <- cbind(seurat_object@meta.data, Other_methods_clusters[,-1])

# UMAPs for clusters
current_label <- "CHOIR_1"
plt <- DimPlot(seurat_object, group.by = current_label, reduction = "P0_umap",
               shuffle = TRUE, label = FALSE, raster = FALSE) +
  theme_void() + NoLegend() + ggtitle(current_label) +
  scale_color_manual(values = alt_color_palette) +
  coord_fixed(ratio = 1) +
  xlim(-20,20) + ylim(-20,20)
plt[[1]]$layers[[1]]$aes_params$alpha = 0.2
plt[[1]]$layers[[1]]$aes_params$size = 0.2
plt[[1]]$layers[[1]]$aes_params$shape = 16
plt

# Specific cell lines
current_label <- "T47D"
seurat_object$Cell_line <- grepl(current_label, seurat_object$Ground_truth)
plt <- DimPlot(seurat_object, group.by = "Cell_line", reduction = "P0_umap",
               shuffle = TRUE, label = FALSE, raster = FALSE) +
  theme_void() + NoLegend() + ggtitle("") +
  scale_color_manual(values = c("grey", "#9EB7FA")) +
  coord_fixed(ratio = 1) +
  xlim(-20,20) + ylim(-20,20)
plt[[1]]$layers[[1]]$aes_params$alpha = 0.2
plt[[1]]$layers[[1]]$aes_params$size = 0.2
plt[[1]]$layers[[1]]$aes_params$shape = 16
plt

# Entropy of cluster accuracy (H accuracy)
H_accuracies <- data.frame(label = NULL,
                           n_clusters = NULL,
                           H_acc = NULL)

for (current_label in colnames(seurat_object@meta.data)[c(1,2,13:132)]) {
  print(current_label)
  M_clusters <- unique(seurat_object@meta.data[, current_label])
  M <- dplyr::n_distinct(seurat_object@meta.data[, current_label])
  sum_M = 0
  for (j in 1:M) {
    cluster_j <- M_clusters[j]
    cluster_j_cells <- seurat_object@meta.data$CellID[seurat_object@meta.data[, current_label] == cluster_j]
    ground_truth_clusters_j <- unique(seurat_object@meta.data$Ground_truth[seurat_object@meta.data$CellID %in% cluster_j_cells])
    N <- length(ground_truth_clusters_j) # of true clusters w/in cluster j
    n_j <- length(cluster_j_cells) # of cells in cluster j
    sum_N = 0
    for (k in 1:N) {
      cluster_k <- ground_truth_clusters_j[k]
      cluster_k_cells <- seurat_object@meta.data$CellID[seurat_object@meta.data$Ground_truth == cluster_k]
      n_k = length(intersect(cluster_j_cells, cluster_k_cells)) # of cells in cluster k w/in cluster j
      p_k = n_k/n_j
      stat = p_k*log(p_k)
      sum_N = sum_N + stat
    }
    sum_M = sum_M + sum_N
  }
  H_accuracy <- (sum_M/M)*(-1)
  print(H_accuracy)
  H_accuracies <- rbind(H_accuracies, data.frame(label = current_label,
                                                 n_clusters = M,
                                                 H_acc = H_accuracy))
}
write.csv(H_accuracies, paste0(input_dir, "/H_accuracies.csv"))

# Plot entropy of cluster accuracy
H_accuracies %>%
  mutate(Default = ifelse(label %in% c("CHOIR_1", "Cytocipher_1", "GiniClust3_1", 
                                       "SCCAF_1", "scSHC_4", "Seurat_1"), TRUE, FALSE)) %>%
  mutate(Method = substr(label, 1, 5)) %>%
  ggplot(aes(x = Method, y = H_acc, color = Method, fill = Default)) +
  theme_classic() +
  geom_beeswarm(shape = 21, cex = 0.5) +
  scale_fill_manual(values = c("white", "black"))

# Projection onto independent dataset
cell_line <- "T47D"
cells <- dplyr::filter(seurat_object@meta.data, grepl(cell_line, Ground_truth))$CellID
clust <- seurat_object@meta.data[cells, "CHOIR_1"]
kinker_pool_obj <- subset(seurat_object, cells = cells)

if (cell_line == "A375") {
  projection_data <- readRDS(paste0(input_dir, "/Yang_A375_RNA_matrix_raw.rds"))
} else if (cell_line == "T47D") {
  projection_data <- readRDS(paste0(input_dir, "/T47D_RNA_matrix_raw.rds"))
}
projection_data <- Azimuth:::ConvertEnsembleToSymbol(mat = projection_data, species = "human")
shared_features <- intersect(rownames(kinker_pool_obj)[rowSums(kinker_pool_obj@assays$RNA$counts) > 0], rownames(projection_data)[rowSums(projection_data) > 0])
kinker_pool_obj <- subset(kinker_pool_obj, features = shared_features)
projection_data <- projection_data[shared_features,]
projection <- CreateSeuratObject(projection_data)
projection <- projection %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>% RunPCA() %>%
  RunUMAP(dims = 1:30, return.model = TRUE) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters()
projection.anchors <- FindTransferAnchors(reference = projection, query = kinker_pool_obj,
                                          dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = projection.anchors, 
                            refdata = as.character(projection$seurat_clusters),
                            dims = 1:30)
kinker_pool_obj <- AddMetaData(kinker_pool_obj, metadata = predictions)
kinker_pool_obj <- MapQuery(anchorset = projection.anchors, reference = projection, query = kinker_pool_obj,
                            refdata = list(celltype = "seurat_clusters"), reference.reduction = "pca", 
                            reduction.model = "umap")
# Plot projection
projection@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  ggplot(aes(x = umap_1, y = umap_2)) +
  theme_void() +
  geom_point(color = "grey", alpha = 0.2) +
  geom_point(data = data.frame(kinker_pool_obj@reductions$ref.umap@cell.embeddings), 
             aes(x = refUMAP_1, y = refUMAP_2,
                 color = factor(kinker_pool_obj@meta.data$CHOIR_1))) +
  scale_color_manual(values = alt_color_palette[c(75, 249, 346)]) +
  NoLegend() +
  coord_fixed(ratio = 1)

# Import MAGIC-imputed data
if (cell_line == "A375") {
  seurat_magic <- readRDS(paste0(input_dir, "/Yang_A375_RNA_matrix_raw_magic.rds"))
} else if (cell_line == "T47D") {
  seurat_magic <- readRDS(paste0(input_dir, "/T47D_RNA_matrix_raw_magic.rds"))
}
seurat_magic[["umap"]] <- CreateDimReducObject(embeddings = projection@reductions$umap@cell.embeddings,
                                               assay = "RNA")
kinker_magic <- readRDS(paste0(input_dir, "/Kinker_all_raw_magic.rds"))
kinker_magic <- subset(kinker_magic, cells = colnames(kinker_pool_obj))
kinker_magic[["umap"]] <- CreateDimReducObject(embeddings = kinker_pool_obj@reductions$ref.umap@cell.embeddings,
                                               assay = "RNA")

# Feature plots
color_gradient <- c("#DDDDDD", "#DDDDDD","#FFFC5A",
                    "#FFD32B","#FF9965",
                    "#FD619D","#C732D5","#6E34FC",
                    "#1632FB","#021EA9", "#000436")
color_values <- c(0, 0.08, 0.081, 
                  0.195, 0.31,
                  0.425, 0.54, 0.655,
                  0.77, 0.885, 1)
df <- NULL
gene <- "DEPTOR"

# Independent dataset
ggplot(df) +
  theme_void() +
  geom_point(aes(x = seurat_magic@reductions$umap@cell.embeddings[,1],
                 y = seurat_magic@reductions$umap@cell.embeddings[,2],
                 color = seurat_magic@assays$MAGIC_RNA$data[gene,]), shape = 16) +
  coord_fixed(ratio = 1) +
  scale_color_gradientn(colors = color_gradient,
                        values = color_values) +
  NoLegend()

# Projected Kinker et al. data
ggplot(df) +
  theme_void() +
  geom_point(aes(x = projection@reductions$umap@cell.embeddings[,1],
                 y = projection@reductions$umap@cell.embeddings[,2]), 
             color = "white", shape = 16) +
  geom_point(aes(x = kinker_magic@reductions$umap@cell.embeddings[,1],
                 y = kinker_magic@reductions$umap@cell.embeddings[,2],
                 color = c(as.matrix(kinker_magic@assays$MAGIC_RNA$data[gene,]))), shape = 16) +
  coord_fixed(ratio = 1) +
  scale_color_gradientn(colors = color_gradient,
                        values = color_values) +
  NoLegend()

# Computational time
CHOIR_time <- read.csv(paste0(input_dir, "/compiled_time_CHOIR_parameter_list_RNA_real.csv"))
CHOIR_time$method <- "CHOIR"
other_methods_time <- read.csv(paste0(input_dir, "/compiled_time_cluster_parameter_list_RNA_real.csv"))
other_methods_time <- other_methods_time
time_df <- rbind(CHOIR_time, other_time)
time_df <- time_df %>%
  mutate(default = ifelse((method == "scSHC" & parameter_index == 4) |
                            (method != "scSHC" & parameter_index == 1), TRUE, FALSE))
time_df %>%
  ggplot(aes(x = method, y = time_secs/3600, color = default)) +
  geom_beeswarm()