# ---------------------------------------------------------------------------
# This script contains the benchmarking analysis of the Hao et al. 2021
# CITE-seq human PBMC dataset
# ---------------------------------------------------------------------------
library(BPCells)
library(Seurat)
library(CHOIR)
library(dplyr)
library(ggplot2)
library(progress)     
library(Libra)       
library(edgeR)       
library(Matrix)
library(pbmcapply)

# Set color palette
color_palette <- CHOIR::CHOIRpalette(100)
alt_color_palette <- CHOIR::CHOIRpalette(500)

# Import data
input_dir <- "/gladstone/corces/lab/users/cpetersen/cluster_benchmarking/Hao_2021_human_PBMCs_CITEseq"
count_matrix <- open_matrix_dir(dir = paste0(input_dir, "/bpcells"))

# Create Seurat Object
options(Seurat.object.assay.version = 'v5')
seurat_object <- CreateSeuratObject(counts = count_matrix)

# Add metadata
cell_metadata <- read.csv(paste0(input_dir, "/cell_metadata.csv"))
seurat_object@meta.data <- cbind(seurat_object@meta.data, cell_metadata)

# Normalize
seurat_object <- NormalizeData(seurat_object)

# Import CHOIR UMAP coordinates
choir_data <- readRDS(paste0(input_dir, "/CHOIR_output_1.rds"))
umap_coords <- RunUMAP(choir_data$reduction$P0_reduction)
# Add UMAP coordinates to Seurat object
seurat_object[["P0_umap"]] <- CreateDimReducObject(embeddings = umap_coords@cell.embeddings,
                                                   assay = "RNA")

# Feature plot
# Import MAGIC imputed features
magic_object <- readRDS(paste0(input_dir, "/seurat_magic.rds"))
magic_object[["P0_umap"]] <- CreateDimReducObject(embeddings = umap_coords@cell.embeddings,
                                                   assay = "RNA")
color_gradient <- c("#DDDDDD", "#DDDDDD","#FFFC5A",
                    "#FFD32B","#FF9965",
                    "#FD619D","#C732D5","#6E34FC",
                    "#1632FB","#021EA9", "#000436")
color_values <- c(0, 0.08, 0.081, 
                  0.195, 0.31,
                  0.425, 0.54, 0.655,
                  0.77, 0.885, 1)
df <- NULL
gene <- "CCR7"
ggplot(df) +
  theme_void() +
  geom_point(aes(x = magic_object@reductions$P0_umap@cell.embeddings[,1],
                 y = magic_object@reductions$P0_umap@cell.embeddings[,2],
                 color = magic_object@assays$MAGIC_RNA$data[gene,]), shape = 16,
             size = 0.2, alpha = 0.2) +
  coord_fixed(ratio = 1) +
  scale_color_gradientn(colors = color_gradient,
                        values = color_values) +
  NoLegend()
min(magic_object@assays$MAGIC_RNA$data[gene,])
max(magic_object@assays$MAGIC_RNA$data[gene,])

# Import clusters
# CHOIR
CHOIR_clusters <- read.csv(paste0(input_dir, 
                                  "/compiled_clusters_CHOIR_parameter_list_RNA_real_batch.csv"))
rownames(CHOIR_clusters) <- CHOIR_clusters$CellID
colnames(CHOIR_clusters) <- sub("Parameters", "CHOIR", colnames(CHOIR_clusters))
CHOIR_clusters <- CHOIR_clusters[seurat_object@meta.data$CellID,]
# Other methods
Other_methods_clusters <- read.csv(paste0(input_dir, 
                                          "/compiled_clusters_cluster_parameter_list_RNA_real_batch.csv"))
rownames(Other_methods_clusters) <- Other_methods_clusters$CellID
Other_methods_clusters <- Other_methods_clusters[seurat_object@meta.data$CellID,]

# Add clustering results to Seurat object
seurat_object@meta.data <- cbind(seurat_object@meta.data, CHOIR_clusters[,-1])
seurat_object@meta.data <- cbind(seurat_object@meta.data, Other_methods_clusters[,-1])

# UMAPs for clusters
current_label <- "CHOIR_1"
plt <- DimPlot(seurat_object, group.by = current_label, reduction = "P0_umap",
               shuffle = TRUE, label = TRUE, raster = FALSE) +
  theme_void() +
  NoLegend() +
  ggtitle("") +
  scale_color_manual(values = color_palette) +
  coord_fixed(ratio = 1) +
  xlim(-16,16) +
  ylim(-16,16)
plt[[1]]$layers[[1]]$aes_params$alpha = 0.2
plt[[1]]$layers[[1]]$aes_params$size = 0.05
plt[[1]]$layers[[1]]$aes_params$shape = 16
plt

# Assess pseudobulk differential expression between pairwise cluster comparisons
# Import read count matrix for surface marker protein data
matrix_data <- readMM(paste0(input_dir, "/ADT_matrix.mtx"))
barcodes <- read.table(file = paste0(input_dir, "/ADT_barcodes.tsv"),
                       sep = "\t", header = FALSE)
barcodes[,1] <- paste0("CITEseq_3P", "#", barcodes[,1])
features <- read.table(file = paste0(input_dir, "/ADT_features.tsv"),
                       sep = "\t", header = FALSE)
rownames(matrix_data) <- features$V1
colnames(matrix_data) <- barcodes$V1

# Import cell metadata
cell_metadata <- read.csv(paste0(input_dir, "/cell_metadata.csv"))
matrix_data <- matrix_data[, cell_metadata$CellID]
# Create Seurat Object
adt_object <- CreateSeuratObject(counts = matrix_data)
adt_object@meta.data <- cbind(adt_object@meta.data, cell_metadata)
# Normalize
norm_data1 <- sweep(matrix_data,2,colSums(matrix_data),`/`)
norm_data2 <- sweep(norm_data1,2,10000,`*`)
norm_data3 <- log1p(norm_data2)
adt_object@assays$RNA[["data"]] <- norm_data3
# Add PCA dimensions
adt_object[["P0_pca"]] <- CreateDimReducObject(embeddings = choir_data$reduction$P0_reduction,
                                                  assay = "RNA")

# For each set of clusters
current_clusters <- cbind(CHOIR_clusters[,-1], Other_methods_clusters[,-1])
rownames(current_clusters) <- rownames(CHOIR_clusters)
output <- data.frame(index = NULL,
                     n_clusters = NULL,
                     min_DEGs = NULL)
for (i in 1:ncol(current_clusters)) {
  # Add to object
  adt_object@meta.data$Clusters <- current_clusters[adt_object@meta.data$CellID, i]
  
  if (is.na(adt_object@meta.data$Clusters[1])) {
    output_i <- data.frame(index = i,
                         n_clusters = NA,
                         min_DEGs = NA)
  } else if (dplyr::n_distinct(adt_object@meta.data$Clusters) == 1) {
    output_i <- data.frame(index = i,
                         n_clusters = 1,
                         min_DEGs = NA)
  } else {
    # Get centroid distances
    centroid_distances <- CHOIR:::.getCentroidDistance(reduction = adt_object@reductions$P0_pca@cell.embeddings,
                                                                    clusters = adt_object@meta.data$Clusters)
    centroid_distances[!lower.tri(centroid_distances)] <- NA
    centroid_distances <- data.frame(cluster1 = rep(rownames(centroid_distances), ncol(centroid_distances)),
                                     cluster2 = rep(colnames(centroid_distances), each = nrow(centroid_distances)),
                                     centroid_distance = c(centroid_distances)) %>%
      dplyr::filter(!is.na(centroid_distance)) %>%
      arrange(centroid_distance)
    
    print(paste0(nrow(centroid_distances), " pairs"))
    # Get pairwise DEGs in ADT matrix
    degs <- pbmclapply(1:min(50,nrow(centroid_distances)), function(j) {
      adt_object@meta.data <- adt_object@meta.data %>%
        mutate(subset = ifelse(Clusters %in% c(centroid_distances$cluster1[j],
                                               centroid_distances$cluster2[j]), "yes", "no"),
               cluster_j = ifelse(Clusters == centroid_distances$cluster1[j], 1, 0))
      # cluster DE
      try(cluster_DE_j <- run_de(adt_object, min_reps = 0, min_cells = 0,
                                 replicate_col = 'donor',
                                 cell_type_col = "subset",
                                 label_col = 'cluster_j',
                                 de_method = 'edgeR'))
      if (exists("cluster_DE_j")) {
        cluster_DE_j %>% arrange(p_val) %>% head()
        n_degs <- cluster_DE_j %>% dplyr::filter(p_val_adj < 0.05, cell_type == "yes") %>% nrow()
      } else {
        n_degs <- 0
      }
      return(n_degs)
    }, mc.cores = 8)
    degs <- unlist(degs)
    output_i <- data.frame(index = i,
                           n_clusters = dplyr::n_distinct(adt_object@meta.data$Clusters),
                           min_DEGs = min(degs))
  }
  output <- rbind(output, output_i)
}
write.csv(output, paste0(input_dir, "/pairwise_DEPs.csv"))

# Plot
overclustered_counts <- output %>% 
  mutate(min_DEGs = ifelse(n_clusters == 1, 1, min_DEGs)) %>%
  mutate(overclustered = ifelse(min_DEGs == 0, TRUE, FALSE)) %>%
  group_by(Method, overclustered) %>% summarise(n = n())

overclustered_counts %>%
  group_by(Method) %>%
  mutate(n_total = sum(n)) %>%
  ungroup() %>%
  data.frame() %>%
  ggplot(aes(x = Method, y = n/n_total, fill = overclustered)) +
  theme_classic() +
  geom_bar(stat = "identity") +
  ylab("Percent of parameters tested")

DEGs_adj %>% 
  mutate(min_DEGs = ifelse(n_clusters == 1, 1, min_DEGs)) %>%
  mutate(overclustered = ifelse(min_DEGs == 0, TRUE, FALSE)) %>%
  mutate(method_index = paste0(Method, index)) %>%
  mutate(default = ifelse(method_index %in% c("CHOIR_1", "Cytocipher_1", "GiniClust3_1",
                                              "SCCAF_1", "scSHC_4", "Seurat_1"), TRUE, FALSE)) %>%
  dplyr::filter(overclustered == FALSE) %>%
  ggplot(aes(x = Method, y = n_clusters, color = default)) +
  theme_classic() +
  geom_beeswarm() +
  ylim(0,25)

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