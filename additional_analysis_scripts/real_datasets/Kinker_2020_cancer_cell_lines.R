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

# ---------------------------------------------------------------------------
# Addition during revision: proliferation module score comparisons
# ---------------------------------------------------------------------------

# Import feature sets and create key to rename Ensembl IDs
feature_files <- list.files(paste0(input_dir, "/features"))
feature_key <- data.frame(Ensembl = rownames(count_matrix),
                          Gene = NA)
for (i in 1:length(feature_files)) {
  print(i)
  feature_file_i <- read.table(paste0(input_dir, "/features/", feature_files[i]))
  rownames(feature_file_i) <- feature_file_i$V1
  for (j in 1:nrow(feature_key)) {
    if (j%%1000 == 0) {
      print(paste0(j, "/", nrow(feature_key)))
    }
    if (feature_key$Ensembl[j] %in% rownames(feature_file_i) & is.na(feature_key$Gene[j])) {
      feature_key$Gene[j] <- feature_file_i[feature_key$Ensembl[j], "V2"]
    }
  }
}

feature_key$Gene <- make.unique(feature_key$Gene)
rownames(count_matrix) <- feature_key$Gene

# Instances where CHOIR identifies multiple clusters within a cell line
CHOIR_divisions_by_cell_line <- seurat_object@meta.data %>% dplyr::group_by(Ground_truth2, CHOIR_1) %>% summarise(n = n()) %>% data.frame()
#write.csv(CHOIR_divisions_by_cell_line, paste0(input_dir, "/CHOIR_divisions_by_cell_line.csv"), row.names = FALSE)

# Cell proliferation module score
# Gene set: CELL_PROLIFERATION_GO_0008283
# Source: https://www.gsea-msigdb.org/gsea/msigdb/cards/CELL_PROLIFERATION_GO_0008283
module_genes <- read.csv(paste0(input_dir, 
                                "/module_gene_sets/CELL_PROLIFERATION_GO_0008283.v2023.2.Hs.tsv"))
module_genes <- names(module_genes)

# Which module genes are not found
module_genes[!(module_genes %in% rownames(seurat_object))]
# Add gene aliases when appropriate
# Add WARS instead of WARS1
module_genes <- c(module_genes, "WARS")
# Check whether duplicate gene names need to be included (they don't)
for (i in 1:length(module_genes)) {
  if (module_genes[i] %in% rownames(seurat_object) &
      length(rownames(seurat_object)[grepl(paste0("^", module_genes[i]), rownames(seurat_object))]) > 1) {
    print(rownames(seurat_object)[grepl(paste0("^", module_genes[i]), rownames(seurat_object))])
  }
}
# Subset module genes to genes in data
module_genes <- intersect(module_genes, rownames(seurat_object))

# Create module score
seurat_object <- AddModuleScore(
  object = seurat_object,
  features = list(module_genes),
  name = 'Proliferation_Module')

# Plot module score for different clusters within a cell line
df = NULL
ggplot(df) +
  theme_classic() +
  geom_density(aes(x = seurat_object$Proliferation_Module1[rownames(dplyr::filter(seurat_object@meta.data, 
                                                                                  CHOIR_1 == 13))]), color = alt_color_palette[13]) +
  geom_density(aes(x = seurat_object$Proliferation_Module1[rownames(dplyr::filter(seurat_object@meta.data, 
                                                                                  CHOIR_1 == 29))]), color = alt_color_palette[29]) +
  xlab("Module Score") +
  ylab("Density") +
  xlim(0, 0.15) +
  ylim(0, 38)
ggsave(paste0(input_dir, "/CHOIR_13_29_proliferation_module.pdf"), height = 5, width = 5)

ggplot(df) +
  theme_classic() +
  geom_density(aes(x = seurat_object$Proliferation_Module1[rownames(dplyr::filter(seurat_object@meta.data, 
                                                                                  CHOIR_1 == 346))]), color = "#8E8AE4") +
  geom_density(aes(x = seurat_object$Proliferation_Module1[rownames(dplyr::filter(seurat_object@meta.data, 
                                                                                  CHOIR_1 == 75))]), color = alt_color_palette[75]) +
  xlab("Module Score") +
  ylab("Density") +
  xlim(0, 0.15) +
  ylim(0, 38)
ggsave(paste0(input_dir, "/CHOIR_346_75_proliferation_module.pdf"), height = 5, width = 5)

# Compare proliferation scores across intra-cell line clusters

# Mann-Whitney test
CHOIR_divisions_by_cell_line_rm_dups <- CHOIR_divisions_by_cell_line %>% group_by(CHOIR_1) %>% slice_max(n, n = 1) %>% ungroup() %>%
  arrange(Ground_truth2)

comparisons <- data.frame(cluster1 = NULL,
                          cluster2 = NULL,
                          Ground_truth2 = NULL,
                          p = NULL)
for (i in 1:190) {
  print(i)
  n_rows <- nrow(dplyr::filter(CHOIR_divisions_by_cell_line_rm_dups, Ground_truth2 == i))
  if (n_rows > 1) {
    for (j in 1:(n_rows - 1)) {
      cluster_j <- dplyr::filter(CHOIR_divisions_by_cell_line_rm_dups, Ground_truth2 == i)$CHOIR_1[j]
      for (k in (j + 1):n_rows) {
        cluster_k <- dplyr::filter(CHOIR_divisions_by_cell_line_rm_dups, Ground_truth2 == i)$CHOIR_1[k]
        p_val <- wilcox.test(seurat_object$Proliferation_Module1[rownames(dplyr::filter(seurat_object@meta.data, 
                                                                                        CHOIR_1 == cluster_j))],
                             seurat_object$Proliferation_Module1[rownames(dplyr::filter(seurat_object@meta.data, 
                                                                                        CHOIR_1 == cluster_k))])$p.value
        comparisons <- rbind(comparisons,
                             data.frame(cluster1 = cluster_j,
                                        cluster2 = cluster_k,
                                        Ground_truth2 = i,
                                        p = p_val))
      }
    }
  }
}

comparisons$p_adj <- p.adjust(comparisons$p, method = "fdr")
comparison_results <- comparisons %>% 
  group_by(Ground_truth2) %>% summarise(sum_diff = sum(p < 0.05),
                                        sum_diff_adj = sum(p_adj < 0.05),
                                        sum_same = sum(p >= 0.05),
                                        sum_same_adj = sum(p_adj >= 0.05)) %>% 
  data.frame()


mean_prolif_df <- seurat_object@meta.data %>% data.frame() %>% dplyr::group_by(CHOIR_1, Ground_truth, Ground_truth2) %>%
  summarise(n = n(),
            mean_prolif = mean(Proliferation_Module1))
prolif_merged <- merge(mean_prolif_df, comparison_results, by = "Ground_truth2", all = TRUE)
n_clusters <- CHOIR_divisions_by_cell_line_rm_dups %>% group_by(Ground_truth2) %>% summarise(n_clusters = n())
prolif_merged <- merge(prolif_merged, n_clusters, by = "Ground_truth2", all = TRUE)
prolif_merged <-prolif_merged %>% arrange(Ground_truth) %>%
  mutate(sig = ifelse(is.na(sum_diff), NA,
                      ifelse(sum_diff_adj > 0, TRUE, FALSE)))
write.csv(prolif_merged, paste0(input_dir, "/CHOIR_cell_line_divisions.csv"))

# Percentages
sum(comparisons$p_adj < 0.05)/nrow(comparisons)

n_rows <- CHOIR_divisions_by_cell_line_rm_dups %>% group_by(Ground_truth2) %>% summarise(n_clusters = n()) %>%
  dplyr::filter(n_clusters > 1) %>% nrow()
pdf(file = paste0(input_dir, "/Kinker_pie_chart_CHOIR_1.pdf"), height = 5, width = 5)
pie(c(sum(comparison_results$sum_diff_adj > 0), 
      sum(comparison_results$sum_diff_adj == 0),
      190 - n_rows))
dev.off()
sum(comparison_results$sum_diff_adj > 0)/190
sum(comparison_results$sum_diff_adj == 0)/190
(190 - n_rows)/190

0.4421053/(0.4421053 + 0.1526316)

# Cell lines with multiple clusters but not significantly diff cell proliferation module score
current_label <- "CCFSTTG1_CENTRAL_NERVOUS_SYSTEM"
seurat_object$Cell_line <- grepl(current_label, seurat_object$Ground_truth)
plt <- DimPlot(seurat_object, group.by = "Cell_line", reduction = "P0_umap",
               shuffle = TRUE, label = FALSE, raster = FALSE) +
  theme_void() + NoLegend() + ggtitle("") +
  scale_color_manual(values = c("grey", "#FA7B8E")) +
  coord_fixed(ratio = 1) +
  xlim(-20,20) + ylim(-20,20)
plt[[1]]$layers[[1]]$aes_params$alpha = 0.2
plt[[1]]$layers[[1]]$aes_params$size = 0.2
plt[[1]]$layers[[1]]$aes_params$shape = 16
plt
ggsave(paste0(input_dir, "/Plots/Cell_line_", current_label, "_new.png"),
       height = 5, width = 5, dpi = 2400)

table(subset(seurat_object, subset = Ground_truth == "CCFSTTG1_CENTRAL_NERVOUS_SYSTEM")@meta.data$CHOIR_1)
DimPlot(subset(seurat_object, subset = Ground_truth == "CCFSTTG1_CENTRAL_NERVOUS_SYSTEM"),
        group.by = "CHOIR_1") +
  coord_fixed(ratio = 1)
ggsave(paste0(input_dir, "/CCFSTTG1.pdf"), width = 5, height = 5)

DimPlot(subset(seurat_object),
        group.by = "CHOIR_1", label = TRUE) +
  scale_y_continuous(limits = c(-6,-4)) +
  scale_x_continuous(limits = c(-14,-12.5)) +
  NoLegend()

VlnPlot(subset(seurat_object, subset = Ground_truth == "CCFSTTG1_CENTRAL_NERVOUS_SYSTEM"),
        features = "IL1RAPL1", group.by = "CHOIR_1")
ggsave(paste0(input_dir, "/IL1RAPL1.pdf"), width = 5, height = 5)

mks <- FindMarkers(seurat_object, ident.1 = 341, ident.2 = 340, group.by = "CHOIR_1")
head(mks, 20)

# A375 and T47D

FeaturePlot_scCustom(subset(seurat_object, subset = CHOIR_1 %in% c(13,29)), 
                     "Proliferation_Module1", split.by = "CHOIR_1", order = FALSE)

FeaturePlot_scCustom(subset(seurat_object, subset = CHOIR_1 %in% c(75,346)), 
                     "Proliferation_Module1", split.by = "CHOIR_1", order = FALSE)

# Projection

cell_line <- "A375"
cells <- dplyr::filter(seurat_object@meta.data, grepl(cell_line, Ground_truth))$CellID
length(cells)
clust <- seurat_object@meta.data[cells, "CHOIR_1"]
table(clust)
kinker_pool_obj_temp <- subset(seurat_object, cells = cells)

kinker_pool_obj_mat <- as(kinker_pool_obj_temp@assays$RNA@layers$counts, "dgCMatrix")
rownames(kinker_pool_obj_mat) <- rownames(kinker_pool_obj_temp)
colnames(kinker_pool_obj_mat) <- colnames(kinker_pool_obj_temp)
kinker_pool_obj <- CreateSeuratObject(kinker_pool_obj_mat, meta.data = kinker_pool_obj_temp@meta.data)

projection_data <- readRDS(paste0(input_dir, "/Yang_A375_RNA_matrix_noSoupX_raw.rds"))
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

kinker_pool_obj <- kinker_pool_obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>% RunPCA() %>%
  RunUMAP(dims = 1:30, return.model = TRUE) %>%
  FindNeighbors(dims = 1:30)

kinker_pool_obj$seurat_clusters <- kinker_pool_obj$CHOIR_1

DimPlot(projection, group.by = "seurat_clusters", label = TRUE)
DimPlot(kinker_pool_obj, group.by = "seurat_clusters", label = TRUE)

projection.anchors <- FindTransferAnchors(reference = kinker_pool_obj, query = projection,
                                          dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = projection.anchors, 
                            refdata = as.character(kinker_pool_obj$seurat_clusters),
                            dims = 1:30)
projection <- AddMetaData(projection, metadata = predictions)
projection <- MapQuery(anchorset = projection.anchors, reference = kinker_pool_obj, query = projection,
                       refdata = list(celltype = "seurat_clusters"), reference.reduction = "pca", 
                       reduction.model = "umap")

p1 <- DimPlot(kinker_pool_obj, reduction = "umap", group.by = "CHOIR_1", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations") + 
  xlim(-7,9) + ylim(-7,7)

p2 <- DimPlot(projection, reduction = "umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels") +
  xlim(-7,9) + ylim(-7,7)
p1 + p2

projection <- AddModuleScore(projection, features = list(module_genes), name = 'Proliferation_Module')

df = NULL
ggplot(df) +
  theme_classic() +
  geom_density(aes(x = projection$Proliferation_Module1[rownames(dplyr::filter(projection@meta.data, 
                                                                               predicted.celltype == 13))]), color = alt_color_palette[13]) +
  geom_density(aes(x = projection$Proliferation_Module1[rownames(dplyr::filter(projection@meta.data, 
                                                                               predicted.celltype == 29))]), color = alt_color_palette[29]) +
  xlab("Module Score") +
  ylab("Density") +
  xlim(0, 0.15) +
  ylim(0, 38)
ggsave(paste0(input_dir, "/CHOIR_13_29_proliferation_module_projection.pdf"), height = 5, width = 5)

wilcox.test(projection$Proliferation_Module1[rownames(dplyr::filter(projection@meta.data, 
                                                                    predicted.celltype == 13))],
            projection$Proliferation_Module1[rownames(dplyr::filter(projection@meta.data, 
                                                                    predicted.celltype == 29))])
# T47D

cell_line <- "T47D"
cells <- dplyr::filter(seurat_object@meta.data, grepl(cell_line, Ground_truth))$CellID
length(cells)
clust <- seurat_object@meta.data[cells, "CHOIR_1"]
table(clust)
kinker_pool_obj_temp <- subset(seurat_object, cells = cells)

kinker_pool_obj_mat <- as(kinker_pool_obj_temp@assays$RNA@layers$counts, "dgCMatrix")
rownames(kinker_pool_obj_mat) <- rownames(kinker_pool_obj_temp)
colnames(kinker_pool_obj_mat) <- colnames(kinker_pool_obj_temp)
kinker_pool_obj <- CreateSeuratObject(kinker_pool_obj_mat, meta.data = kinker_pool_obj_temp@meta.data)

projection_data <- readRDS(paste0(input_dir, "/T47D_RNA_matrix_noSoupX_raw.rds"))
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

kinker_pool_obj <- kinker_pool_obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>% RunPCA() %>%
  RunUMAP(dims = 1:30, return.model = TRUE) %>%
  FindNeighbors(dims = 1:30)

kinker_pool_obj$seurat_clusters <- kinker_pool_obj$CHOIR_1

DimPlot(projection, group.by = "seurat_clusters", label = TRUE)
DimPlot(kinker_pool_obj, group.by = "seurat_clusters", label = TRUE)

projection.anchors <- FindTransferAnchors(reference = kinker_pool_obj, query = projection,
                                          dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = projection.anchors, 
                            refdata = as.character(kinker_pool_obj$seurat_clusters),
                            dims = 1:30)
projection <- AddMetaData(projection, metadata = predictions)
projection <- MapQuery(anchorset = projection.anchors, reference = kinker_pool_obj, query = projection,
                       refdata = list(celltype = "seurat_clusters"), reference.reduction = "pca", 
                       reduction.model = "umap")

p1 <- DimPlot(kinker_pool_obj, reduction = "umap", group.by = "CHOIR_1", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations") + 
  xlim(-7,9) + ylim(-7,7)

p2 <- DimPlot(projection, reduction = "umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels") +
  xlim(-7,9) + ylim(-7,7)
p1 + p2

projection <- AddModuleScore(projection, features = list(module_genes), name = 'Proliferation_Module')

df = NULL
ggplot(df) +
  theme_classic() +
  geom_density(aes(x = projection$Proliferation_Module1[rownames(dplyr::filter(projection@meta.data, 
                                                                               predicted.celltype == 346))]), color = "#8E8AE4") +
  geom_density(aes(x = projection$Proliferation_Module1[rownames(dplyr::filter(projection@meta.data, 
                                                                               predicted.celltype == 75))]), color = alt_color_palette[75]) +
  xlab("Module Score") +
  ylab("Density") +
  xlim(0, 0.15) +
  ylim(0, 38)
ggsave(paste0(input_dir, "/CHOIR_75_346_proliferation_module_projection.pdf"), height = 5, width = 5)

wilcox.test(projection$Proliferation_Module1[rownames(dplyr::filter(projection@meta.data, 
                                                                    predicted.celltype == 346))],
            projection$Proliferation_Module1[rownames(dplyr::filter(projection@meta.data, 
                                                                    predicted.celltype == 75))])

# Other way projection to get nearest neighbors

cell_line <- "A375"
cells <- dplyr::filter(seurat_object@meta.data, grepl(cell_line, Ground_truth))$CellID
length(cells)
clust <- seurat_object@meta.data[cells, "CHOIR_1"]
table(clust)
kinker_pool_obj_temp <- subset(seurat_object, cells = cells)

kinker_pool_obj_mat <- as(kinker_pool_obj_temp@assays$RNA@layers$counts, "dgCMatrix")
rownames(kinker_pool_obj_mat) <- rownames(kinker_pool_obj_temp)
colnames(kinker_pool_obj_mat) <- colnames(kinker_pool_obj_temp)
kinker_pool_obj <- CreateSeuratObject(kinker_pool_obj_mat, meta.data = kinker_pool_obj_temp@meta.data)

projection_data <- readRDS(paste0(input_dir, "/Yang_A375_RNA_matrix_noSoupX_raw.rds"))
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

kinker_pool_obj <- kinker_pool_obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>% RunPCA() %>%
  RunUMAP(dims = 1:30, return.model = TRUE) %>%
  FindNeighbors(dims = 1:30)

kinker_pool_obj$seurat_clusters <- kinker_pool_obj$CHOIR_1

DimPlot(projection, group.by = "seurat_clusters", label = TRUE)
DimPlot(kinker_pool_obj, group.by = "seurat_clusters", label = TRUE)

projection.anchors <- FindTransferAnchors(reference = projection, query = kinker_pool_obj,
                                          dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = projection.anchors, 
                            refdata = as.character(projection$seurat_clusters),
                            dims = 1:30)
kinker_pool_obj <- AddMetaData(kinker_pool_obj, metadata = predictions)
kinker_pool_obj <- MapQuery(anchorset = projection.anchors, reference = projection, query = kinker_pool_obj,
                            refdata = list(celltype = "seurat_clusters"), reference.reduction = "pca", 
                            reduction.model = "umap")

p1 <- DimPlot(projection, reduction = "pca", group.by = "seurat_clusters", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")

p2 <- DimPlot(kinker_pool_obj, reduction = "ref.pca", group.by = "CHOIR_1", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

pca_projection <- data.frame(projection@reductions$pca@cell.embeddings)[,1:30]
pca_kinker <- data.frame(kinker_pool_obj@reductions$ref.pca@cell.embeddings)
colnames(pca_kinker) <- paste0("PC_", 1:30)
pca_combined <- rbind(pca_projection, pca_kinker)
combined_distance <- dist(pca_combined, method = "euclidean")
combined_distance_mat <- as.matrix(combined_distance)

neighbor_df <- data.frame(Kinker_cell = colnames(kinker_pool_obj),
                          Kinker_cluster = kinker_pool_obj$CHOIR_1,
                          Nearest_neighbor = rep("None", 593),
                          Distance = rep(0, 593))
for (i in 1:nrow(neighbor_df)) {
  print(i)
  cell_i <- neighbor_df$Kinker_cell[i]
  ranked_neighbors <- sort(combined_distance_mat[cell_i, 1:4794])
  neighbor_used <- TRUE
  j <- 1
  while (neighbor_used == TRUE) {
    distance_i <- ranked_neighbors[j]
    neighbor_i <- names(ranked_neighbors)[j]
    if (neighbor_i %in% neighbor_df$Nearest_neighbor) {
      print(neighbor_df$Distance[neighbor_df$Nearest_neighbor == neighbor_i])
      if (distance_i < neighbor_df$Distance[neighbor_df$Nearest_neighbor == neighbor_i]) {
        neighbor_used <- FALSE
        print(neighbor_df$Nearest_neighbor[neighbor_df$Nearest_neighbor == neighbor_i])
        neighbor_df$Nearest_neighbor[neighbor_df$Nearest_neighbor == neighbor_i] <- "Redo"
      } else {
        j <- j + 1
      }
    } else {
      neighbor_used <- FALSE
    }
  }
  neighbor_df$Distance[i] <- distance_i
  neighbor_df$Nearest_neighbor[i] <- neighbor_i
}

for (i in 1:nrow(neighbor_df)) {
  if (neighbor_df$Nearest_neighbor[i] == "Redo") {
    print(i)
    cell_i <- neighbor_df$Kinker_cell[i]
    ranked_neighbors <- sort(combined_distance_mat[cell_i, 1:4794])
    neighbor_used <- TRUE
    j <- 1
    while (neighbor_used == TRUE) {
      distance_i <- ranked_neighbors[j]
      neighbor_i <- names(ranked_neighbors)[j]
      if (neighbor_i %in% neighbor_df$Nearest_neighbor) {
        print(neighbor_df$Distance[neighbor_df$Nearest_neighbor == neighbor_i])
        if (distance_i < neighbor_df$Distance[neighbor_df$Nearest_neighbor == neighbor_i]) {
          neighbor_used <- FALSE
          print(neighbor_df$Nearest_neighbor[neighbor_df$Nearest_neighbor == neighbor_i])
          neighbor_df$Nearest_neighbor[neighbor_df$Nearest_neighbor == neighbor_i] <- "Redo"
        } else {
          j <- j + 1
        }
      } else {
        neighbor_used <- FALSE
      }
    }
    neighbor_df$Distance[i] <- distance_i
    neighbor_df$Nearest_neighbor[i] <- neighbor_i
  }
}
# Re-run for loop above until no instances of "Redo"
dplyr::filter(neighbor_df, Nearest_neighbor == "Redo")

projection$CellID <- colnames(projection)
projection_subset <- subset(projection, subset = CellID %in% neighbor_df$Nearest_neighbor)

rownames(neighbor_df) <- neighbor_df$Nearest_neighbor
projection_subset$CHOIR_1 <- neighbor_df[colnames(projection_subset),]$Kinker_cluster

projection_subset <- AddModuleScore(projection_subset, features = list(module_genes), name = 'Proliferation_Module')

df = NULL
ggplot(df) +
  theme_classic() +
  geom_density(aes(x = projection_subset$Proliferation_Module1[rownames(dplyr::filter(projection_subset@meta.data, 
                                                                                      CHOIR_1 == 13))]), color = alt_color_palette[13]) +
  geom_density(aes(x = projection_subset$Proliferation_Module1[rownames(dplyr::filter(projection_subset@meta.data, 
                                                                                      CHOIR_1 == 29))]), color = alt_color_palette[29]) +
  xlab("Module Score") +
  ylab("Density") +
  xlim(0, 0.15) +
  ylim(0, 50)
ggsave(paste0(input_dir, "/CHOIR_13_29_proliferation_module_projection_nearest_neighbors.pdf"), height = 5, width = 5)

x <- wilcox.test(projection_subset$Proliferation_Module1[rownames(dplyr::filter(projection_subset@meta.data, 
                                                                                CHOIR_1 == 13))],
                 projection_subset$Proliferation_Module1[rownames(dplyr::filter(projection_subset@meta.data, 
                                                                                CHOIR_1 == 29))])
x$p.value

p1 <- DimPlot(projection_subset, reduction = "umap", group.by = "CHOIR_1", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Yang") + coord_fixed(ratio = 1)

p2 <- DimPlot(kinker_pool_obj, reduction = "ref.umap", group.by = "CHOIR_1", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Kinker") + coord_fixed(ratio = 1)
p1 + p2
ggsave(paste0(input_dir, "/CHOIR_13_29_proliferation_module_projection_nearest_neighbors_UMAP.pdf"), height = 5, width = 10)

min(projection_subset$Proliferation_Module1[rownames(dplyr::filter(projection_subset@meta.data, 
                                                                   CHOIR_1 == 13))])
max(projection_subset$Proliferation_Module1[rownames(dplyr::filter(projection_subset@meta.data, 
                                                                   CHOIR_1 == 13))])
min(projection_subset$Proliferation_Module1[rownames(dplyr::filter(projection_subset@meta.data, 
                                                                   CHOIR_1 == 29))])
max(projection_subset$Proliferation_Module1[rownames(dplyr::filter(projection_subset@meta.data, 
                                                                   CHOIR_1 == 29))])
# T47D

cell_line <- "T47D"
cells <- dplyr::filter(seurat_object@meta.data, grepl(cell_line, Ground_truth))$CellID
length(cells)
clust <- seurat_object@meta.data[cells, "CHOIR_1"]
table(clust)
kinker_pool_obj_temp <- subset(seurat_object, cells = cells)

kinker_pool_obj_mat <- as(kinker_pool_obj_temp@assays$RNA@layers$counts, "dgCMatrix")
rownames(kinker_pool_obj_mat) <- rownames(kinker_pool_obj_temp)
colnames(kinker_pool_obj_mat) <- colnames(kinker_pool_obj_temp)
kinker_pool_obj <- CreateSeuratObject(kinker_pool_obj_mat, meta.data = kinker_pool_obj_temp@meta.data)

projection_data <- readRDS(paste0(input_dir, "/T47D_RNA_matrix_noSoupX_raw.rds"))
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

kinker_pool_obj <- kinker_pool_obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>% RunPCA() %>%
  RunUMAP(dims = 1:30, return.model = TRUE) %>%
  FindNeighbors(dims = 1:30)

kinker_pool_obj$seurat_clusters <- kinker_pool_obj$CHOIR_1

DimPlot(projection, group.by = "seurat_clusters", label = TRUE)
DimPlot(kinker_pool_obj, group.by = "seurat_clusters", label = TRUE)

projection.anchors <- FindTransferAnchors(reference = projection, query = kinker_pool_obj,
                                          dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = projection.anchors, 
                            refdata = as.character(projection$seurat_clusters),
                            dims = 1:30)
kinker_pool_obj <- AddMetaData(kinker_pool_obj, metadata = predictions)
kinker_pool_obj <- MapQuery(anchorset = projection.anchors, reference = projection, query = kinker_pool_obj,
                            refdata = list(celltype = "seurat_clusters"), reference.reduction = "pca", 
                            reduction.model = "umap")

p1 <- DimPlot(projection, reduction = "pca", group.by = "seurat_clusters", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")

p2 <- DimPlot(kinker_pool_obj, reduction = "ref.pca", group.by = "CHOIR_1", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

kinker_pool_obj <- subset(kinker_pool_obj, subset = CHOIR_1 != 249)

pca_projection <- data.frame(projection@reductions$pca@cell.embeddings)[,1:30]
pca_kinker <- data.frame(kinker_pool_obj@reductions$ref.pca@cell.embeddings)
colnames(pca_kinker) <- paste0("PC_", 1:30)
pca_combined <- rbind(pca_projection, pca_kinker)
combined_distance <- dist(pca_combined, method = "euclidean")
combined_distance_mat <- as.matrix(combined_distance)

neighbor_df <- data.frame(Kinker_cell = colnames(kinker_pool_obj),
                          Kinker_cluster = kinker_pool_obj$CHOIR_1,
                          Nearest_neighbor = rep("None", length(colnames(kinker_pool_obj))),
                          Distance = rep(0, length(colnames(kinker_pool_obj))))
for (i in 1:nrow(neighbor_df)) {
  print(i)
  cell_i <- neighbor_df$Kinker_cell[i]
  ranked_neighbors <- sort(combined_distance_mat[cell_i, 1:length(colnames(projection))])
  neighbor_used <- TRUE
  j <- 1
  while (neighbor_used == TRUE) {
    distance_i <- ranked_neighbors[j]
    neighbor_i <- names(ranked_neighbors)[j]
    if (neighbor_i %in% neighbor_df$Nearest_neighbor) {
      print(neighbor_df$Distance[neighbor_df$Nearest_neighbor == neighbor_i])
      if (distance_i < neighbor_df$Distance[neighbor_df$Nearest_neighbor == neighbor_i]) {
        neighbor_used <- FALSE
        print(neighbor_df$Nearest_neighbor[neighbor_df$Nearest_neighbor == neighbor_i])
        neighbor_df$Nearest_neighbor[neighbor_df$Nearest_neighbor == neighbor_i] <- "Redo"
      } else {
        j <- j + 1
      }
    } else {
      neighbor_used <- FALSE
    }
  }
  neighbor_df$Distance[i] <- distance_i
  neighbor_df$Nearest_neighbor[i] <- neighbor_i
}

for (i in 1:nrow(neighbor_df)) {
  if (neighbor_df$Nearest_neighbor[i] == "Redo") {
    print(i)
    cell_i <- neighbor_df$Kinker_cell[i]
    ranked_neighbors <- sort(combined_distance_mat[cell_i, 1:4794])
    neighbor_used <- TRUE
    j <- 1
    while (neighbor_used == TRUE) {
      distance_i <- ranked_neighbors[j]
      neighbor_i <- names(ranked_neighbors)[j]
      if (neighbor_i %in% neighbor_df$Nearest_neighbor) {
        print(neighbor_df$Distance[neighbor_df$Nearest_neighbor == neighbor_i])
        if (distance_i < neighbor_df$Distance[neighbor_df$Nearest_neighbor == neighbor_i]) {
          neighbor_used <- FALSE
          print(neighbor_df$Nearest_neighbor[neighbor_df$Nearest_neighbor == neighbor_i])
          neighbor_df$Nearest_neighbor[neighbor_df$Nearest_neighbor == neighbor_i] <- "Redo"
        } else {
          j <- j + 1
        }
      } else {
        neighbor_used <- FALSE
      }
    }
    neighbor_df$Distance[i] <- distance_i
    neighbor_df$Nearest_neighbor[i] <- neighbor_i
  }
}
# Re-run for loop above until no instances of "Redo"
dplyr::filter(neighbor_df, Nearest_neighbor == "Redo")

projection$CellID <- colnames(projection)
projection_subset <- subset(projection, subset = CellID %in% neighbor_df$Nearest_neighbor)

rownames(neighbor_df) <- neighbor_df$Nearest_neighbor
projection_subset$CHOIR_1 <- neighbor_df[colnames(projection_subset),]$Kinker_cluster

projection_subset <- AddModuleScore(projection_subset, features = list(module_genes), name = 'Proliferation_Module')

df = NULL
ggplot(df) +
  theme_classic() +
  geom_density(aes(x = projection_subset$Proliferation_Module1[rownames(dplyr::filter(projection_subset@meta.data, 
                                                                                      CHOIR_1 == 346))]), color = "#8E8AE4") +
  geom_density(aes(x = projection_subset$Proliferation_Module1[rownames(dplyr::filter(projection_subset@meta.data, 
                                                                                      CHOIR_1 == 75))]), color = alt_color_palette[75]) +
  xlab("Module Score") +
  ylab("Density") +
  xlim(0, 0.15) +
  ylim(0, 38)
ggsave(paste0(input_dir, "/CHOIR_75_346_proliferation_module_projection_nearest_neighbors.pdf"), height = 5, width = 5)

wilcox.test(projection_subset$Proliferation_Module1[rownames(dplyr::filter(projection_subset@meta.data, 
                                                                           CHOIR_1 == 75))],
            projection_subset$Proliferation_Module1[rownames(dplyr::filter(projection_subset@meta.data, 
                                                                           CHOIR_1 == 346))])

p1 <- DimPlot(projection_subset, reduction = "umap", group.by = "CHOIR_1", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Independent") + coord_fixed(ratio = 1)

p2 <- DimPlot(kinker_pool_obj, reduction = "ref.umap", group.by = "CHOIR_1", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Kinker") + coord_fixed(ratio = 1)
p1 + p2
ggsave(paste0(input_dir, "/CHOIR_75_346_proliferation_module_projection_nearest_neighbors_UMAP.pdf"), height = 5, width = 10)

min(projection_subset$Proliferation_Module1[rownames(dplyr::filter(projection_subset@meta.data, 
                                                                   CHOIR_1 == 75))])
max(projection_subset$Proliferation_Module1[rownames(dplyr::filter(projection_subset@meta.data, 
                                                                   CHOIR_1 == 75))])
min(projection_subset$Proliferation_Module1[rownames(dplyr::filter(projection_subset@meta.data, 
                                                                   CHOIR_1 == 346))])
max(projection_subset$Proliferation_Module1[rownames(dplyr::filter(projection_subset@meta.data, 
                                                                   CHOIR_1 == 346))])

# ---------------------------------------------------------------------------
# Addition during revision: downsampling analysis
# ---------------------------------------------------------------------------

## Create downsampled datasets (these are then clustered by each method just as the original dataset was)

cell_metadata_original <- read.csv(paste0(input_dir, "/cell_metadata.csv"))
cell_lines <- unique(cell_metadata_original$Ground_truth)
# Current cell line
for (j in 1:190) {
  current_cell_line <- cell_lines[j]
  other_ids <- dplyr::filter(cell_metadata_original, Ground_truth != current_cell_line)$CellID
  current_ids <- dplyr::filter(cell_metadata_original, Ground_truth == current_cell_line)$CellID
  set.seed(442266)
  current_ids <- sample(current_ids, 50)
  cell_metadata_current <- cell_metadata_original %>% dplyr::filter(CellID %in% c(other_ids, current_ids))
  write.csv(cell_metadata_current, 
            paste0(input_dir, "/new_metadata_files/cell_metadata_", 
                   current_cell_line, "_", i, ".csv"), 
            row.names = FALSE)
}

## Evaluate underclustering in original runs

input_dir <- "/gladstone/corces/lab/users/cpetersen/cluster_benchmarking/parameter_lists"
cell_lines <- read.csv("/gladstone/corces/lab/users/cpetersen/cluster_benchmarking/Kinker_2020_all/cell_line_names.csv")
input_dir <- "/gladstone/corces/lab/users/cpetersen/cluster_benchmarking/Kinker_2020_all"
parameter_list_directory <- "/gladstone/corces/lab/users/cpetersen/cluster_benchmarking/parameter_lists"
metadata_file <- paste0("cell_metadata.csv")

# Repeat across each method
cluster_file <- paste0("compiled_clusters_CHOIR_parameter_list_RNA_real.csv")
method <- "CHOIR"
parameters <- c(1:7)

id <- paste0(method, "_", parameters)
cell_id <- c()
for (i in id) {
  cell_id <- c(cell_id, paste0(i, "_", cell_lines$C))
}

df <- data.frame(method_name = method,
                 method_id = rep(paste0(method, "_", parameters), each = 190),
                 index = cell_id,
                 works = 1)
rownames(df) <- df$index

for (i in 1:length(cell_lines$C)) {
  print(i)
  cell_line_i <- cell_lines$C[i]
  
  # Add metadata
  cell_metadata <- read.csv(paste0(input_dir, "/intermediate_files/preprocessed/", metadata_file))
  
  # Import clusters
  try(clusters <- read.csv(paste0(input_dir, 
                                  "/intermediate_files/clusters/", cluster_file)))
  if (exists("clusters")) {
    rownames(clusters) <- clusters$CellID
    colnames(clusters) <- sub("Parameters", method, colnames(clusters))
    clusters <- clusters[cell_metadata$CellID,]
    
    # Add clustering results
    cell_metadata <- cbind(cell_metadata, clusters[,-1])
    
    parameter_names <- paste0(method, "_", parameters)
    run_parameters <- rep(1, length(parameter_names))
    # Assess whether cell line is underclustered
    for (p in 1:length(parameter_names)) {
      parameter <- parameter_names[p]
      if (is.na(cell_metadata[1,parameter])) {
        df[paste0(parameter, "_", cell_line_i), "works"] <- 0
      } else {
        clusters_i <- cell_metadata %>% 
          dplyr::filter(Ground_truth == cell_line_i) %>%
          group_by(!!sym(parameter)) %>%
          summarise(n = n()) %>%
          arrange(-n) %>% 
          data.frame()
        primary_cluster_i <- clusters_i[1, parameter]
        cell_line_i_n <- sum(clusters_i$n)
        primary_cluster_i_n <- clusters_i$n[1]
        
        # Get other cell lines in primary cluster
        primary_cluster_i_cell_lines <- cell_metadata %>% 
          dplyr::filter(!!sym(parameter) == primary_cluster_i, Ground_truth != cell_line_i) %>%
          group_by(Ground_truth) %>%
          summarise(n = n()) %>%
          arrange(-n) %>% 
          data.frame()
        
        if (nrow(primary_cluster_i_cell_lines) >= 1) {
          # Check what proportion of primary cluster is made up by cell line i
          prop_cell_line_i <- primary_cluster_i_n/(sum(primary_cluster_i_cell_lines$n) + primary_cluster_i_n)
          
          if (prop_cell_line_i < 0.5) {
            df[paste0(parameter, "_", cell_line_i), "works"] <- 0
          }
          # Check if primary cluster is primary cluster for any other cell lines
          for (j in 1:nrow(primary_cluster_i_cell_lines)) {
            cell_line_j <- primary_cluster_i_cell_lines$Ground_truth[j]
            cluster_j_n <- primary_cluster_i_cell_lines$n[j]
            if (cluster_j_n > 2) {
              clusters_j <- cell_metadata %>% 
                dplyr::filter(Ground_truth == cell_line_j) %>%
                group_by(!!sym(parameter)) %>%
                summarise(n = n()) %>%
                arrange(-n) %>% 
                data.frame()
              primary_cluster_j <- clusters_j[1, parameter]
              if (primary_cluster_j == primary_cluster_i) {
                df[paste0(parameter, "_", cell_line_i), "works"] <- 0
              }
            }
          }
        } 
      }
    }
    rm(clusters)
  }
}
# Switch object name by method
choir_orig <- df %>% transmute(method_name, method_id, index, cell_line, works_orig = works)

# Merge
orig <- rbind(choir_orig, cyto_orig)
orig <- rbind(orig, gini_orig)
orig <- rbind(orig, sccaf_orig)
orig <- rbind(orig, scshc_orig)
orig <- rbind(orig, seurat_orig)

## Evaluate underclustering across downsampled runs

# Repeat across each method
cluster_file_prefix <- "compiled_clusters_CHOIR_parameter_list_"
time_file_prefix <- "compiled_time_CHOIR_parameter_list_"
method <- "CHOIR"
parameters <- c(1:7)

id <- paste0(method, "_", parameters)
cell_id <- c()
for (i in id) {
  cell_id <- c(cell_id, paste0(i, "_", cell_lines$C))
}

df <- data.frame(method_name = method,
                 method_id = rep(paste0(method, "_", parameters), each = 190),
                 index = cell_id,
                 n_clusters = NA,
                 n_distinguished = NA,
                 H_accuracy = NA,
                 works = 1,
                 time_secs = NA,
                 cpu_time = NA,
                 maxvmem = NA,
                 mem = NA)
rownames(df) <- df$index

for (i in 1:length(cell_lines$C)) {
  print(i)
  cell_line_i <- cell_lines$C[i]
  
  metadata_file <- paste0("cell_metadata_", cell_line_i, "_50.csv")
  cluster_file <- paste0(cluster_file_prefix, cell_line_i, "_50.csv")
  
  # Add metadata
  cell_metadata <- read.csv(paste0(input_dir, "/intermediate_files/preprocessed/", metadata_file))
  
  # Import clusters
  try(clusters <- read.csv(paste0(input_dir, 
                                  "/intermediate_files/clusters/", cluster_file)))
  if (exists("clusters")) {
    rownames(clusters) <- clusters$CellID
    colnames(clusters) <- sub("Parameters", method, colnames(clusters))
    clusters <- clusters[cell_metadata$CellID,]
    
    # Add clustering results
    cell_metadata <- cbind(cell_metadata, clusters[,-1])
    
    parameter_names <- paste0(method, "_", parameters)
    run_parameters <- rep(1, length(parameter_names))
    
    n_clusters <- apply(clusters[-1], 2, n_distinct)
    df[paste0(parameter_names, "_", cell_line_i), "n_clusters"] <- n_clusters
    
    time <- read.csv(paste0(input_dir, 
                                "/intermediate_files/clusters/", 
                                time_file_prefix, cell_line_i, "_50.csv"))
    
    df[paste0(parameter_names, "_", cell_line_i), "time_secs"] <- time$time_secs
    df[paste0(parameter_names, "_", cell_line_i), "cpu_time"] <- time$total_cpu_time
    df[paste0(parameter_names, "_", cell_line_i), "maxvmem"] <- time$maxvmem
    df[paste0(parameter_names, "_", cell_line_i), "mem"] <- time$mem
    
    # Assess whether cell line is underclustered
    for (p in 1:length(parameter_names)) {
      parameter <- parameter_names[p]
      if (is.na(cell_metadata[1,parameter])) {
        df[paste0(parameter, "_", cell_line_i), "works"] <- 0
      } else {
        clusters_i <- cell_metadata %>% 
          dplyr::filter(Ground_truth == cell_line_i) %>%
          group_by(!!sym(parameter)) %>%
          summarise(n = n()) %>%
          arrange(-n) %>% 
          data.frame()
        primary_cluster_i <- clusters_i[1, parameter]
        cell_line_i_n <- sum(clusters_i$n)
        primary_cluster_i_n <- clusters_i$n[1]
        
        # Get other cell lines in primary cluster
        primary_cluster_i_cell_lines <- cell_metadata %>% 
          dplyr::filter(!!sym(parameter) == primary_cluster_i, Ground_truth != cell_line_i) %>%
          group_by(Ground_truth) %>%
          summarise(n = n()) %>%
          arrange(-n) %>% 
          data.frame()
        
        if (nrow(primary_cluster_i_cell_lines) >= 1) {
          # Check what proportion of primary cluster is made up by cell line i
          prop_cell_line_i <- primary_cluster_i_n/(sum(primary_cluster_i_cell_lines$n) + primary_cluster_i_n)
          
          if (prop_cell_line_i < 0.5) {
            df[paste0(parameter, "_", cell_line_i), "works"] <- 0
          }
          # Check if primary cluster is primary cluster for any other cell lines
          for (j in 1:nrow(primary_cluster_i_cell_lines)) {
            cell_line_j <- primary_cluster_i_cell_lines$Ground_truth[j]
            cluster_j_n <- primary_cluster_i_cell_lines$n[j]
            if (cluster_j_n > 2) {
              clusters_j <- cell_metadata %>% 
                dplyr::filter(Ground_truth == cell_line_j) %>%
                group_by(!!sym(parameter)) %>%
                summarise(n = n()) %>%
                arrange(-n) %>% 
                data.frame()
              primary_cluster_j <- clusters_j[1, parameter]
              if (primary_cluster_j == primary_cluster_i) {
                df[paste0(parameter, "_", cell_line_i), "works"] <- 0
              }
            }
          }
        } 
      }
    }
    rm(clusters)
    rm(time)
  } else {
    parameter_names <- paste0(method, "_", parameters)
    for (p in 1:length(parameter_names)) {
      parameter <- parameter_names[p]
      df[paste0(parameter, "_", cell_line_i), "works"] <- NA
    }
  }
}
# Switch object name by method
choir_50 <- df %>% mutate(n_enough = ifelse(n_clusters >= 190, TRUE, FALSE))
choir_50 <- choir_50 %>% transmute(method_name, method_id, index, cell_line, works_50 = works)

d50 <- rbind(choir_50, cyto_50)
d50 <- rbind(d50, gini_50)
d50 <- rbind(d50, sccaf_50)
d50 <- rbind(d50, scshc_50)
d50 <- rbind(d50, seurat_50)

## Compile results

# Merge
d_merge <- merge(orig, d50, by = c("index", "method_id", "method_name", "cell_line"), all = TRUE)
d_merge <- d_merge %>%
  mutate(works_50 = ifelse(works_orig == 0, NA, works_50))
d_merge %>% dplyr::filter(works_orig != 0) %>% nrow()
d_stats <- d_merge %>% 
  group_by(method_name, method_id) %>% 
  summarise(sum_works_orig = sum(works_orig, na.rm = TRUE),
            sum_works_50 = sum(works_50, na.rm = TRUE)) %>% 
  data.frame()
d_stats <- d_stats %>% 
  tidyr::pivot_longer(cols = c("sum_works_orig", "sum_works_50"), 
                      names_to = "group", 
                      values_to = "works") 
d_stats$group <- factor(d_stats$group, levels = c("sum_works_orig", "sum_works_50"))
d_stats <- d_stats %>%
  mutate(default = ifelse(method_id %in% c("CHOIR_1", "Cytocipher_1", "GiniClust3_1", "SCCAF_1", "scSHC_4", "Seurat_1"),
                          TRUE, FALSE))

d_stats %>%
  ggplot(aes(x = group, y = works, group = method_id, color = method_name, alpha = default, linewidth = default)) +
  theme_classic() +
  geom_line() +
  scale_color_manual(values = CHOIR::CHOIRpalette(6)) +
  scale_alpha_manual(values = c(0.3, 1)) +
  scale_linewidth_manual(values = c(0.3, 1)) +
  scale_y_continuous(breaks = seq(0,190,10))
ggsave("~/Desktop/Kinker_downsample_cell_lines.pdf", width = 5, height = 5)

d_stats %>%
  ggplot(aes(x = group, y = works/190, group = method_id, color = method_name, alpha = default, linewidth = default)) +
  theme_classic() +
  geom_line() +
  scale_color_manual(values = CHOIR::CHOIRpalette(6)) +
  scale_alpha_manual(values = c(0.3, 1)) +
  scale_linewidth_manual(values = c(0.3, 1)) +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
ggsave("~/Desktop/Kinker_downsample_cell_lines_percent.pdf", width = 5, height = 5)

for_table <- d_stats %>% dplyr::filter(group == "sum_works_orig") %>% data.frame() %>%
  mutate(index = ifelse(method_name == "CHOIR", as.numeric(sub("CHOIR_", "", method_id)),
                        ifelse(method_name == "Cytocipher", as.numeric(sub("Cytocipher_", "", method_id)),
                               ifelse(method_name == "GiniClust3", as.numeric(sub("GiniClust3_", "", method_id)),
                                      ifelse(method_name == "SCCAF", as.numeric(sub("SCCAF_", "", method_id)),
                                             ifelse(method_name == "scSHC", as.numeric(sub("scSHC_", "", method_id)),
                                                    as.numeric(sub("Seurat_", "", method_id)))))))) %>% arrange(method_name, index)

write.csv(for_table, "~/Desktop/Kinker_downsample_cell_lines_n_distinguished.csv")


