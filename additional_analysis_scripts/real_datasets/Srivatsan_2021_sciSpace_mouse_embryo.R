# ---------------------------------------------------------------------------
# This script contains the benchmarking analysis of the Srivatsan et al. 2021
# sci-Space mouse embryo dataset
# ---------------------------------------------------------------------------
library(BPCells)
library(Seurat)
library(CHOIR)
library(dplyr)
library(ggplot2)
library(scales)
library(scCustomize)
library(ggforce)
library(ggbeeswarm)

# Set color palette
color_palette <- ClusteringBetaTest::CHOIRpalette(100)
alt_color_palette <- ClusteringBetaTest::CHOIRpalette(500)

# Import data
input_dir <- "/gladstone/corces/lab/users/cpetersen/cluster_benchmarking/Srivatsan_2021_mouse_embryo_spatial"
count_matrix <- open_matrix_dir(dir = paste0(input_dir, "/bpcells"))

# Create Seurat Object
options(Seurat.object.assay.version = 'v5')
seurat_object <- CreateSeuratObject(counts = count_matrix)

# Add metadata
cell_metadata <- read.csv(paste0(input_dir, "/cell_metadata.csv"))
seurat_object@meta.data <- cbind(seurat_object@meta.data, cell_metadata)

# Normalize & find variable features
seurat_object <- NormalizeData(seurat_object)

# Import CHOIR UMAP coordinates
choir_data <- readRDS(paste0(input_dir, "/CHOIR_output_1.rds"))
umap_coords <- RunUMAP(choir_data$reduction$P0_reduction)
# Add UMAP coordinates to Seurat object
seurat_object[["P0_umap"]] <- CreateDimReducObject(embeddings = umap_coords@cell.embeddings,
                                                   assay = "RNA")

# Import imaging data
imaging_data <- readRDS(paste0(input_dir, "/misc/GSE166692_sciSpace_imaging.RDS"))

# Import MAGIC imputed features
magic_object <- readRDS(paste0(input_dir, "/seurat_magic.rds"))

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
seurat_object@meta.data$CHOIR_parent_clusters <- choir_data$clusters$full_tree[seurat_object@meta.data$CellID, "L5"]

# UMAPs for clusters
current_label <- "CHOIR_1"
plt <- DimPlot(seurat_object, group.by = current_label, reduction = "P0_umap",
               shuffle = TRUE, label = FALSE, raster = FALSE) +
  theme_void() + NoLegend() + ggtitle("") +
  scale_color_manual(values = alt_color_palette) +
  coord_fixed(ratio = 1) +
  xlim(-16,16) + ylim(-16,16)
plt[[1]]$layers[[1]]$aes_params$alpha = 0.2
plt[[1]]$layers[[1]]$aes_params$size = 0.05
plt[[1]]$layers[[1]]$aes_params$shape = 16
plt

# UMAPs for each parent cluster
current_label <- 4
seurat_object$CHOIR_parent <- grepl(paste0("P", current_label), seurat_object$CHOIR_parent_clusters)
plt <- DimPlot(seurat_object, group.by = "CHOIR_parent", reduction = "P0_umap",
               order = TRUE, shuffle = FALSE, label = FALSE, raster = FALSE) +
  theme_void() + NoLegend() + ggtitle("") +
  scale_color_manual(values = c("#DDDDDD", color_palette[1])) +
  coord_fixed(ratio = 1) +
  xlim(-16,16) + ylim(-16,16)
plt[[1]]$layers[[1]]$aes_params$alpha = 0.2
plt[[1]]$layers[[1]]$aes_params$size = 0.05
plt[[1]]$layers[[1]]$aes_params$shape = 16
plt

# UMAP for specific clusters
current_cluster <- 39
current_label <- "CHOIR_1"
seurat_object$Cluster <- seurat_object@meta.data[,current_label] == current_cluster
plt <- DimPlot(seurat_object, group.by = "Cluster", reduction = "P0_umap", order = TRUE,
               shuffle = FALSE, label = FALSE, raster = FALSE) +
  theme_void() + NoLegend() + ggtitle("") +
  scale_color_manual(values = c("#DDDDDD", color_palette[current_cluster])) +
  coord_fixed(ratio = 1) +
  xlim(-16,16) + ylim(-16,16)
plt[[1]]$layers[[1]]$aes_params$alpha = 0.2
plt[[1]]$layers[[1]]$aes_params$size = 0.05
plt[[1]]$layers[[1]]$aes_params$shape = 16
plt

# P3 subtree
umap_coords_P3 <- RunUMAP(choir_data$reduction$P3_reduction)
P3_object <- subset(seurat_object, subset = CellID %in% rownames(umap_coords_P3@cell.embeddings))
P3_object[["P3_umap"]] <- CreateDimReducObject(embeddings = umap_coords_P3@cell.embeddings[colnames(P3_object),],
                                               assay = "RNA")
# P3 UMAP
current_label <- "CHOIR_1"
selected_colors <- seurat_object@meta.data %>% 
  dplyr::filter(CHOIR_parent_clusters == "P3_L5_1") %>%
  group_by(!!sym(current_label)) %>% 
  summarise(n = n())
plt <- DimPlot(P3_object, group.by = current_label, reduction = "P3_umap",
               shuffle = TRUE, label = FALSE, raster = FALSE) +
  theme_void() + NoLegend() + ggtitle("") +
  scale_color_manual(values = color_palette[unlist(selected_colors[, current_label])]) +
  coord_fixed(ratio = 1) + xlim(-16,16) + ylim(-16,16)
plt[[1]]$layers[[1]]$aes_params$alpha = 0.2
plt[[1]]$layers[[1]]$aes_params$size = 0.2
plt[[1]]$layers[[1]]$aes_params$shape = 16
plt

# P4 subtree
umap_coords_P4 <- RunUMAP(choir_data$reduction$P4_reduction)
P4_object <- subset(seurat_object, subset = CellID %in% rownames(umap_coords_P4@cell.embeddings))
P4_object[["P4_umap"]] <- CreateDimReducObject(embeddings = umap_coords_P4@cell.embeddings[colnames(P4_object),],
                                               assay = "RNA")
# P4 UMAP
current_label <- "CHOIR_1"
selected_colors <- seurat_object@meta.data %>% 
  dplyr::filter(CHOIR_parent_clusters == "P4_L5_1") %>%
  group_by(!!sym(current_label)) %>% 
  summarise(n = n())
plt <- DimPlot(P4_object, group.by = current_label, reduction = "P4_umap",
               shuffle = TRUE, label = FALSE, raster = FALSE) +
  theme_void() + NoLegend() + ggtitle("") +
  scale_color_manual(values = color_palette[unlist(selected_colors[, current_label])]) +
  coord_fixed(ratio = 1) + xlim(-16,16) + ylim(-16,16)
plt[[1]]$layers[[1]]$aes_params$alpha = 0.2
plt[[1]]$layers[[1]]$aes_params$size = 0.6
plt[[1]]$layers[[1]]$aes_params$shape = 16
plt

# Spatial plots of clusters
# Compile spatial coordinates with cluster IDs
current_label <- "CHOIR_1"
spatial_data <- data.frame(CellID = cell_metadata$CellID,
                           Cluster = seurat_object@meta.data[, current_label],
                           Slide = cell_metadata$slide_id,
                           X_coord = cell_metadata$coords.x1,
                           Y_coord = cell_metadata$coords.x2)
# Slide by slide
for (current_slide in 1:14) {
  print(current_slide)
  slide_outline <- data.frame(imaging_data$hull_polygon[[current_slide]][[1]])
  slide_data <- spatial_data %>% 
    dplyr::filter(Slide == paste0("Slide ", current_slide))
  slide_data <- slide_data %>% 
    group_by(Cluster) %>%
    mutate(n_cells_cluster = n()) %>%
    ungroup()
  slide_data <- slide_data %>%
    group_by(X_coord, Y_coord, Cluster, n_cells_cluster) %>%
    summarise(n_spot = n()) %>%
    mutate(percent_clust = n_spot/n_cells_cluster)
  slide_data <- rbind(slide_data, data.frame(X_coord = 0,
                                             Y_coord = 0, 
                                             Cluster = unique(seurat_object@meta.data[, current_label]), 
                                             n_cells_cluster = 0,
                                             n_spot = 0,
                                             percent_clust = 0))
  # Mode cluster at each spot
  slide_data %>%
    group_by(X_coord, Y_coord) %>%
    slice_max(n_spot, with_ties = TRUE) %>%
    arrange(X_coord, Y_coord, Cluster) %>%
    ggplot(aes(x = X_coord, y = Y_coord, color = as.factor(Cluster))) +
    theme_classic() +
    geom_polygon(data = slide_outline, 
                 aes(x = x_scaled, y = y_scaled), 
                 size = 0.2, fill = "black", color = "black") +
    geom_point(size = 0.5) +
    scale_color_manual(values = color_palette) +
    labs(color = paste0(current_label, " clusters\n(Mode per spot)")) +
    coord_fixed(ratio = 1)
  ggsave(paste0(input_dir, "/slide_", current_slide, "_", current_label, "_mode.pdf"),
         height = 5, width = 5)
}

# By cluster
current_label <- "CHOIR_1"
spatial_data <- data.frame(CellID = cell_metadata$CellID,
                           Cluster = seurat_object@meta.data[, current_label],
                           Slide = cell_metadata$slide_id,
                           X_coord = cell_metadata$coords.x1,
                           Y_coord = cell_metadata$coords.x2)
for (cluster in as.numeric(unique(seurat_object@meta.data[, current_label]))) {
  color_1 = muted(color_3, l = 20, c = 50)
  colfunc <- colorRampPalette(c(color_1, color_3))
  color_2 = colfunc(3)[2]
  color_5 = muted(color_3, l = 90, c = 100)
  colfunc <- colorRampPalette(c(color_3, color_5))
  color_4 = colfunc(3)[2]
  for (current_slide in 1:14) {
    slide_outline <- data.frame(imaging_data$hull_polygon[[current_slide]][[1]])
    slide_data <- spatial_data %>% 
      dplyr::filter(Slide == paste0("Slide ", current_slide))
    slide_data <- slide_data %>% dplyr::filter(Cluster == cluster) %>%
      group_by(X_coord, Y_coord, Cluster) %>%
      summarise(n_cells = n()) %>%
      mutate(n_cells = ifelse(n_cells > 5, 5, n_cells)) %>%
      data.frame()
    color_gradient <- c("black", color_1, color_2, color_3, color_4, color_5)
    slide_data <- rbind(slide_data, data.frame(X_coord = 0,
                                               Y_coord = 0,
                                               Cluster = cluster,
                                               n_cells = c(0,1,2,3,4,5)))
    slide_data %>%
      ggplot(aes(x = X_coord, y = Y_coord, color = factor(n_cells))) +
      theme_classic() +
      geom_polygon(data = slide_outline, 
                   aes(x = x_scaled, y = y_scaled), 
                   size = 0.2, fill = "black", color = "black") +
      geom_point(size = 1, shape = 16) +
      scale_color_manual(values = color_gradient) +
      labs(color = "Number of cells") +
      coord_fixed(ratio = 1) +
      ggtitle(paste0("Cluster ", cluster))
    ggsave(paste0(input_dir, "/Plots/Spatial/", current_label, "/slide_", current_slide, "_", 
                  current_label, "_", cluster, ".pdf"),
           height = 5, width = 5)
  }
}

# Spatial feature plots
gene <- "Lhx6"
current_slide <- 14
spatial_data <- data.frame(CellID = cell_metadata$CellID,
                           Gene = magic_object@assays$MAGIC_RNA@data[gene, ],
                           Slide = cell_metadata$slide_id,
                           X_coord = cell_metadata$coords.x1,
                           Y_coord = cell_metadata$coords.x2)
slide_outline <- data.frame(imaging_data$hull_polygon[[current_slide]][[1]])
slide_data <- spatial_data %>% 
  dplyr::filter(Slide == paste0("Slide ", current_slide)) %>%
  group_by(X_coord, Y_coord) %>%
  summarise(sum_reads = sum(Gene),
            mean_reads = mean(Gene),
            n_cells = n()) %>%
  mutate(reads_per_cell = sum_reads/n_cells) %>%
  dplyr::filter(reads_per_cell != 0)
color_gradient <- c("#DDDDDD", "#DDDDDD","#FFFC5A",
                    "#FFD32B","#FF9965",
                    "#FD619D","#C732D5","#6E34FC",
                    "#1632FB","#021EA9", "#000436")
color_values <- c(0, 0.08, 0.08001, 
                  0.195, 0.31,
                  0.425, 0.54, 0.655,
                  0.77, 0.885, 1)
slide_data %>%
  ggplot(aes(x = X_coord, y = Y_coord, color = mean_reads)) +
  theme_void() +
  geom_polygon(data = slide_outline, 
               aes(x = x_scaled, y = y_scaled), 
               size = 0.2, fill = "white", color = "black") +
  geom_point(size = 1, shape = 16) +
  scale_color_gradientn(colors = color_gradient,
                        values = color_values) +
  labs(color = "Expression") +
  ggtitle(gene) +
  coord_fixed(ratio = 1)
ggsave(paste0(input_dir, "/Plots/Spatial/slide_", current_slide, "_", 
              gene, ".pdf"),
       height = 5, width = 5)

# Spatial distance from mean by cluster
spatial_distance <- data.frame(Label = NULL, CellID = NULL, Cluster = NULL,
                               Slide = NULL, X_coord = NULL, Y_coord = NULL,
                               min_X = NULL, max_X = NULL, min_Y = NULL, max_Y = NULL,
                               slide_area = NULL, mean_X = NULL, mean_Y = NULL,
                               n_cells = NULL, dist_from_centroid = NULL)
for (current_label in c("CHOIR", "Cytocipher", "GiniClust3",
                        "SCCAF", "scSHC", "Seurat")) {
  if (current_label == "CHOIR") {
    params <- 1:7
  } else if (current_label == "Cytocipher") {
    params <- 1:16
  } else if (current_label == "GiniClust3") {
    params <- 1:29
  } else if (current_label == "SCCAF") {
    params <- 1:13
  } else if (current_label == "scSHC") {
    params <- 1:13
  } else if (current_label == "Seurat") {
    params <- 1:28
  } 
  for (p in 1:length(params)) {
    current_param <- paste0(current_label, "_", params[p])
    print(current_param)
    spatial_data <- data.frame(Label = current_label,
                               Label_Param = current_param,
                               CellID = cell_metadata$CellID,
                               Cluster = seurat_object@meta.data[, current_param],
                               Slide = cell_metadata$slide_id,
                               X_coord = cell_metadata$coords.x1,
                               Y_coord = cell_metadata$coords.x2)
    spatial_data <- spatial_data %>% group_by(Slide) %>%
      mutate(min_X = min(X_coord), max_X = max(X_coord),
             min_Y = min(Y_coord), max_Y = max(Y_coord),
             slide_area = (max_X - min_X)*(max_Y - min_Y))
    spatial_data <- spatial_data %>% group_by(Slide, Cluster) %>%
      mutate(mean_X = mean(X_coord), mean_Y = mean(Y_coord), n_cells = n())
    spatial_data <- spatial_data %>%
      mutate(dist_from_centroid = sqrt((X_coord - mean_X)^2 + (Y_coord - mean_Y)^2))
    spatial_distance <- rbind(spatial_distance, spatial_data)
  }
}
# Slide-aware label summary
label_summary <- spatial_distance %>% group_by(Label, Label_Param, Cluster, Slide, slide_area) %>%
  summarise(mean_distance = mean(dist_from_centroid),
            cluster_n_cells = sum(n_cells)) %>%
  mutate(focal = ifelse(pi*(mean_distance^2) < 0.05*slide_area, TRUE, FALSE))
label_summary %>% group_by(Label, Label_Param, Slide) %>%
  summarise(n_clusters = n(),
            sum_focal = sum(focal)) %>%
  mutate(percent_focal = sum_focal/n_clusters,
         default = ifelse(Label_Param %in% c("CHOIR_1", "Cytocipher_1", "GiniClust3_1",
                                             "SCCAF_1", "scSHC_4", "Seurat_1"), 
                          TRUE, FALSE)) %>%
  group_by(Label, Label_Param, default) %>%
  summarise(mean_percent_focal = mean(percent_focal)) %>%
  arrange(default) %>%
  ggplot(aes(x = Label, y = mean_percent_focal, color = default)) +
  theme_classic() +
  geom_beeswarm()
ggsave(paste0(input_dir, "/percent_focal_5percent_slide_area_median.pdf"), width = 5, height = 5)

percent_focal_df <- label_summary %>% group_by(Label, Label_Param, Slide) %>%
  summarise(n_clusters = n(),
            sum_focal = sum(focal)) %>%
  mutate(percent_focal = sum_focal/n_clusters,
         default = ifelse(Label_Param %in% c("CHOIR_1", "Cytocipher_1", "GiniClust3_1",
                                             "SCCAF_1", "scSHC_4", "Seurat_1"), 
                          TRUE, FALSE)) %>%
  group_by(Label, Label_Param, default) %>%
  summarise(mean_percent_focal = mean(percent_focal))

# Which are the focal clusters?
focal_clusters <- label_summary %>% dplyr::filter(Label_Param == "CHOIR_1",
                                                  Slide == "Slide 14",
                                                  focal == TRUE)

# Feature plots for parent clusters
P4_magic <- subset(magic_object, subset = CellID %in% rownames(umap_coords_P4@cell.embeddings))
P4_magic[["P4_umap"]] <- CreateDimReducObject(embeddings = umap_coords_P4@cell.embeddings[colnames(P4_magic),],
                                               assay = "RNA")

P3_magic <- subset(magic_object, subset = CellID %in% rownames(umap_coords_P3@cell.embeddings))
P3_magic[["P3_umap"]] <- CreateDimReducObject(embeddings = umap_coords_P3@cell.embeddings[colnames(P3_magic),],
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
gene <- "Col1a1"
ggplot(df) +
  theme_void() +
  geom_point(aes(x = P4_magic@reductions$P4_umap@cell.embeddings[,1],
                 y = P4_magic@reductions$P4_umap@cell.embeddings[,2],
                 color = P4_magic3@assays$MAGIC_RNA$data[gene,]), shape = 16) +
  coord_fixed(ratio = 1) +
  scale_color_gradientn(colors = color_gradient,
                        values = color_values) +
  NoLegend()
ggsave(paste0(input_dir, "/P4_", gene, ".pdf"))

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