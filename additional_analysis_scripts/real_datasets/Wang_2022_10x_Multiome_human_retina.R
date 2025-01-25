# ---------------------------------------------------------------------------
# This script contains the analysis and plots for the Wang et al. 2022
# 10x Multiome human retina dataset
# ---------------------------------------------------------------------------
library(BPCells)
library(Seurat)
library(dplyr)
library(ggplot2)
library(Azimuth)
library(scales)
library(scCustomize)
library(ggforce)
library(ggbeeswarm)
library(Matrix)
library(pbmcapply)
library(clustree)
library(ArchR)
library(CHOIR)
library(countsplit)

# Set color palette
color_palette <- ClusteringBetaTest::CHOIRpalette(100)

input_dir <- "/gladstone/corces/lab/users/cpetersen/cluster_benchmarking/Wang_2022_multiome"

raw_matrix <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/RNA_matrix_raw.rds"))

# Create Seurat Object
options(Seurat.object.assay.version = 'v5')
seurat_object <- CreateSeuratObject(counts = raw_matrix)

# Add metadata
cell_metadata <- read.csv(paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata.csv"))
seurat_object <- subset(seurat_object, cells = cell_metadata$CellID)
rownames(cell_metadata) <- cell_metadata$CellID
cell_metadata <- cell_metadata[colnames(seurat_object),]
identical(colnames(seurat_object), cell_metadata$CellID)
seurat_object@meta.data <- cbind(seurat_object@meta.data, cell_metadata)
seurat_object@meta.data <- seurat_object@meta.data[,-c(2:3)]

# Normalize & find variable features
seurat_object <- NormalizeData(seurat_object)

# Import CHOIR UMAP coordinates
parameter_list <- "CHOIR_parameter_list_multiome_batch_revision_actual"
choir_data <- readRDS(paste0(input_dir, "/intermediate_files/clusters/", parameter_list, "/CHOIR_output_1.rds"))
umap_coords <- RunUMAP(choir_data$reduction$P0_reduction)
# Add UMAP coordinates to Seurat object
umap_coords@cell.embeddings <- umap_coords@cell.embeddings[colnames(seurat_object),]
seurat_object[["P0_umap"]] <- CreateDimReducObject(embeddings = umap_coords@cell.embeddings,
                                                   assay = "RNA")

# Clustering tree
clustree(choir_data$clusters$full_tree, prefix = "L")
ggsave(paste0(input_dir, "/clustering_tree_new.pdf"), width = 45, height = 10)

# N_clusters for sup table
n_clusters_per_param <- data.frame(param_id = NULL,
                                   n_clusters = NULL)

# Import clusters
# CHOIR
parameter_list <- "CHOIR_parameter_list_multiome_batch"
CHOIR_clusters <- read.csv(paste0(input_dir, 
                                  "/intermediate_files/clusters/compiled_clusters_", parameter_list, ".csv"))
rownames(CHOIR_clusters) <- CHOIR_clusters$CellID
colnames(CHOIR_clusters) <- sub("Parameters", "CHOIR", colnames(CHOIR_clusters))
CHOIR_clusters <- CHOIR_clusters[seurat_object@meta.data$CellID,]
current_n_clusters <- apply(CHOIR_clusters, 2, n_distinct)
n_clusters_per_param <- rbind(n_clusters_per_param, data.frame(param_id = names(current_n_clusters)[-1],
                                                               n_clusters = current_n_clusters[-1]))
# CHOIR - ATAC
parameter_list <- "CHOIR_parameter_list_ATAC_batch"
CHOIR_clusters_ATAC <- read.csv(paste0(input_dir, 
                                       "/intermediate_files/clusters/compiled_clusters_", parameter_list, ".csv"))
rownames(CHOIR_clusters_ATAC) <- CHOIR_clusters_ATAC$CellID
colnames(CHOIR_clusters_ATAC) <- sub("Parameters", "CHOIR_ATAC", colnames(CHOIR_clusters_ATAC))
CHOIR_clusters_ATAC <- CHOIR_clusters_ATAC[seurat_object@meta.data$CellID,]
current_n_clusters <- apply(CHOIR_clusters_ATAC, 2, n_distinct)
n_clusters_per_param <- rbind(n_clusters_per_param, data.frame(param_id = names(current_n_clusters)[-1],
                                                               n_clusters = current_n_clusters[-1]))
# CHOIR - RNA
parameter_list <- "CHOIR_parameter_list_RNA_batch_final_new"
CHOIR_clusters_RNA <- read.csv(paste0(input_dir, 
                                      "/intermediate_files/clusters/compiled_clusters_", parameter_list, ".csv"))
rownames(CHOIR_clusters_RNA) <- CHOIR_clusters_RNA$CellID
colnames(CHOIR_clusters_RNA) <- sub("Parameters", "CHOIR_RNA", colnames(CHOIR_clusters_RNA))
CHOIR_clusters_RNA <- CHOIR_clusters_RNA[seurat_object@meta.data$CellID,]
current_n_clusters <- apply(CHOIR_clusters_RNA, 2, n_distinct)
n_clusters_per_param <- rbind(n_clusters_per_param, data.frame(param_id = names(current_n_clusters)[-1],
                                                               n_clusters = current_n_clusters[-1]))
# CHOIR - RNA - countsplit
parameter_list <- "CHOIR_parameter_list_RNA_batch_countsplit"
CHOIR_clusters_RNA_countsplit <- read.csv(paste0(input_dir, 
                                                 "/intermediate_files/clusters/compiled_clusters_", parameter_list, ".csv"))
rownames(CHOIR_clusters_RNA_countsplit) <- CHOIR_clusters_RNA_countsplit$CellID
colnames(CHOIR_clusters_RNA_countsplit) <- sub("Parameters", "CHOIR_RNA_countsplit", colnames(CHOIR_clusters_RNA_countsplit))
CHOIR_clusters_RNA_countsplit <- CHOIR_clusters_RNA_countsplit[seurat_object@meta.data$CellID,]
current_n_clusters <- apply(CHOIR_clusters_RNA_countsplit, 2, n_distinct)
n_clusters_per_param <- rbind(n_clusters_per_param, data.frame(param_id = names(current_n_clusters)[-1],
                                                               n_clusters = current_n_clusters[-1]))
# ArchR and Signac - ATAC
Other_clusters_ATAC <- read.csv(paste0(input_dir, 
                                       "/outputs/compiled_clusters_cluster_parameter_list_short_ATAC_batch.csv"))
rownames(Other_clusters_ATAC) <- Other_clusters_ATAC$CellID
colnames(Other_clusters_ATAC) <- sub("Parameters", "Other", colnames(Other_clusters_ATAC))
colnames(Other_clusters_ATAC)[c(8,9,10,11,12,14,15,16,17,18,19,20,21,22,23)] <- sub("Other", "ArchR_ATAC", colnames(Other_clusters_ATAC)[c(8,9,10,11,12,14,15,16,17,18,19,20,21,22,23)])
colnames(Other_clusters_ATAC) <- sub("Other", "Signac_ATAC", colnames(Other_clusters_ATAC))
Other_clusters_ATAC <- Other_clusters_ATAC[seurat_object@meta.data$CellID,]
current_n_clusters <- apply(Other_clusters_ATAC[,grepl("ArchR", colnames(Other_clusters_ATAC))], 2, n_distinct)
n_clusters_per_param <- rbind(n_clusters_per_param, data.frame(param_id = names(current_n_clusters),
                                                               n_clusters = current_n_clusters))
current_n_clusters <- apply(Other_clusters_ATAC[,grepl("Signac", colnames(Other_clusters))], 2, n_distinct)
current_n_clusters <- current_n_clusters[order(as.numeric(sub("Signac_ATAC_", "", names(current_n_clusters))))]
n_clusters_per_param <- rbind(n_clusters_per_param, data.frame(param_id = names(current_n_clusters),
                                                               n_clusters = current_n_clusters))
# Seurat - RNA
parameter_list <- "cluster_parameter_list_Seurat_RNA_batch"
Seurat_clusters_RNA <- read.csv(paste0(input_dir, 
                                       "/intermediate_files/clusters/compiled_clusters_", parameter_list, ".csv"))
rownames(Seurat_clusters_RNA) <- Seurat_clusters_RNA$CellID
colnames(Seurat_clusters_RNA) <- sub("Parameters", "Seurat_RNA", colnames(Seurat_clusters_RNA))
Seurat_clusters_RNA <- Seurat_clusters_RNA[seurat_object@meta.data$CellID,]
current_n_clusters <- apply(Seurat_clusters_RNA, 2, n_distinct)
current_n_clusters <- current_n_clusters[-1]
current_n_clusters <- current_n_clusters[order(as.numeric(sub("Seurat_RNA_", "", names(current_n_clusters))))]
n_clusters_per_param <- rbind(n_clusters_per_param, data.frame(param_id = names(current_n_clusters),
                                                               n_clusters = current_n_clusters))
# ArchR and Signac - multiome
Other_clusters <- read.csv(paste0(input_dir, 
                                  "/intermediate_files/clusters/compiled_clusters_cluster_parameter_list_short_multiome_batch.csv"))
rownames(Other_clusters) <- Other_clusters$CellID
colnames(Other_clusters) <- sub("Parameters", "Other", colnames(Other_clusters))
colnames(Other_clusters)[c(8,9,10,11,12,14,15,16,17,18,19,20,21,22,23)] <- sub("Other", "ArchR", colnames(Other_clusters)[c(8,9,10,11,12,14,15,16,17,18,19,20,21,22,23)])
colnames(Other_clusters) <- sub("Other", "Signac", colnames(Other_clusters))
Other_clusters <- Other_clusters[seurat_object@meta.data$CellID,]
current_n_clusters <- apply(Other_clusters[,grepl("ArchR", colnames(Other_clusters))], 2, n_distinct)
n_clusters_per_param <- rbind(n_clusters_per_param, data.frame(param_id = names(current_n_clusters),
                                                               n_clusters = current_n_clusters))
current_n_clusters <- apply(Other_clusters[,grepl("Signac", colnames(Other_clusters))], 2, n_distinct)
current_n_clusters <- current_n_clusters[order(as.numeric(sub("Signac_", "", names(current_n_clusters))))]
n_clusters_per_param <- rbind(n_clusters_per_param, data.frame(param_id = names(current_n_clusters),
                                                               n_clusters = current_n_clusters))
write.csv(n_clusters_per_param, paste0(input_dir, "/n_clusters.csv"))

# Add clustering results to Seurat object
seurat_object@meta.data <- cbind(seurat_object@meta.data, CHOIR_clusters[,-1])
seurat_object@meta.data <- cbind(seurat_object@meta.data, CHOIR_clusters_ATAC[,-1])
seurat_object@meta.data <- cbind(seurat_object@meta.data, CHOIR_clusters_RNA[,-1])
seurat_object@meta.data <- cbind(seurat_object@meta.data, CHOIR_clusters_RNA_countsplit[,-1])
seurat_object@meta.data <- cbind(seurat_object@meta.data, Other_clusters[,-1])
seurat_object@meta.data <- cbind(seurat_object@meta.data, Other_clusters_ATAC[,-1])
seurat_object@meta.data <- cbind(seurat_object@meta.data, Seurat_clusters_RNA[,-1])

cluster_labels_for_table <- seurat_object@meta.data[,c("CellID", "Batch", n_clusters_per_param$param_id)]
write.csv(cluster_labels_for_table, paste0(input_dir, "/cluster_labels.csv"))

seurat_object@meta.data$CHOIR_parent_clusters <- choir_data$clusters$full_tree[seurat_object@meta.data$CellID, "L3"]
seurat_object@meta.data$CHOIR_1_labels <- choir_data$clusters$CHOIR_clusters_0.05[seurat_object@meta.data$CellID, "Record_cluster_label"]

# UMAPs for clusters
current_label <- "CHOIR_1"
plt <- DimPlot(seurat_object, group.by = current_label, reduction = "P0_umap",
               shuffle = TRUE, label = FALSE, raster = FALSE) +
  theme_void() + NoLegend() + ggtitle("") +
  scale_color_manual(values = color_palette) +
  coord_fixed(ratio = 1) +
  xlim(-18,18) + ylim(-18,18)
plt[[1]]$layers[[1]]$aes_params$alpha = 0.2
plt[[1]]$layers[[1]]$aes_params$size = 0.2
plt[[1]]$layers[[1]]$aes_params$shape = 16
plt
ggsave(paste0(input_dir, 
              "/", current_label, ".png"),
       height = 5, width = 5, dpi = 2400)

seurat_object@misc$CHOIR <- choir_data
plotCHOIR(seurat_object, 
          accuracy_scores = TRUE,
          plot_nearest = FALSE)

# Features
parameter_list <- "CHOIR_parameter_list_multiome_batch_revision"
choir_data <- readRDS(paste0(input_dir, "/intermediate_files/clusters/", parameter_list, "/CHOIR_output_1.rds"))
features <- colnames(choir_data$records$feature_importance_records)
vals <- dplyr::filter(choir_data$records$feature_importance_records, 
                      cluster1 == "P3_L3_1")[,3:(ncol(choir_data$records$feature_importance_records))] %>% c() %>% unlist()
rna_vals <- vals[grepl(" _", names(vals))]
genes <- Azimuth:::ConvertEnsembleToSymbol(choir_data$var_features$P0_var_features[[2]]$name, species = "human")
mks <- FindMarkers(seurat_object, ident.1 = 1, ident.2 = 4, group.by = "CHOIR_1", min.pct = 0, logfc.threshold = 0)
gene_key <- choir_data$var_features$P0_var_features[[2]] %>% data.frame()
gene_key$gene_seq_id <- paste0(gene_key$seqnames, " ", gene_key$start, " _")
vals_rna_df <- data.frame(gene_seq_id = names(rna_vals),
                          feature_importance = rna_vals)
vals_rna_df_merged <- merge(vals_rna_df, gene_key, by = "gene_seq_id", all.x = TRUE)
mks$name <- rownames(mks)
vals_rna_df_merged_mks <- merge(vals_rna_df_merged, mks, by = "name", all = TRUE)

vals_rna_df_merged_mks %>%
  ggplot(aes(x = feature_importance, y = abs(avg_log2FC))) +
  geom_point(shape = 21, fill = "white") +
  theme_classic() +
  scale_x_continuous(limits = c(0,0.8), breaks = c(0, 0.2, 0.4, 0.6, 0.8)) +
  scale_y_continuous(limits = c(0,12), breaks = c(0, 4, 8, 12))
ggsave(paste0(input_dir, "/feature_importance_find_markers.pdf"), width = 5, height = 5)

# ArchR Project
archr_proj <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/ATAC_obj_ArchR.rds"))
setwd(paste0(input_dir, "/intermediate_files/arrow_files"))
archr_proj_sub <- subsetArchRProject(archr_proj, cells = colnames(seurat_object), force = TRUE)
archr_proj_sub@cellColData$CHOIR_1 <- paste0("C", seurat_object@meta.data[rownames(archr_proj_sub@cellColData), "CHOIR_1"])
gene_annotation <- getGeneAnnotation(archr_proj)
atac_vals <- vals[!grepl(" _", names(vals))]

chrom <- "chr2"
start_val <- 153871913
end_val <- 154453849
tss <- getTSS(ArchRProj = archr_proj_sub)
gr <- makeGRangesFromDataFrame(data.frame(chr=chrom, start=start_val - 100000, end=end_val + 100000))
gr_highlight <- makeGRangesFromDataFrame(data.frame(chr=chrom, start=start_val, end=end_val))
p <- plotBrowserTrack(archr_proj_sub, region = gr, groupBy = "CHOIR_1",
                      highlight = gr_highlight, maxCells = 50000000000000000000)
grid::grid.newpage()
grid::grid.draw(p)
ggsave(paste0(input_dir, "/atac_track_", chrom, "_", start_val, "_", end_val, ".pdf"), plot = p, width = 5, height = 5)

# MAGIC
magic_object <- readRDS(paste0(input_dir, "/seurat_magic.rds"))

magic_object[["P0_umap"]] <- CreateDimReducObject(embeddings = umap_coords@cell.embeddings[colnames(magic_object),],
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
gene <- "ENSG00000182674"
ggplot(df) +
  theme_void() +
  geom_point(aes(x = magic_object@reductions$P0_umap@cell.embeddings[,1],
                 y = magic_object@reductions$P0_umap@cell.embeddings[,2],
                 color = magic_object@assays$MAGIC_RNA$data[gene,]), shape = 16,
             size = 0.5, alpha = 0.2) +
  coord_fixed(ratio = 1) +
  scale_color_gradientn(colors = color_gradient,
                        values = color_values) +
  NoLegend()
ggsave(paste0(input_dir, "/P0_", gene, "_revision.png"), width = 5, height = 5, dpi = 2400)

min(magic_object@assays$MAGIC_RNA$data[gene,])
max(magic_object@assays$MAGIC_RNA$data[gene,])
gene <- "ENSG00000182674"
min_val <- min(magic_object@assays$MAGIC_RNA$data[gene,])
round(min_val,2)
max_val <- max(magic_object@assays$MAGIC_RNA$data[gene,])
round(max_val,2)
round(0.08*(max_val-min_val) + min_val,2)

# Time
CHOIR_time <- read.csv(paste0(input_dir, "/intermediate_files/clusters/compiled_time_CHOIR_parameter_list_ATAC_batch.csv"))
CHOIR_time$method <- "CHOIR"
time_df <- CHOIR_time
time_df <- time_df %>%
  mutate(default = ifelse(parameter_index == 1, TRUE, FALSE))
time_df %>%
  ggplot(aes(x = method, y = time_secs/3600, color = default)) +
  geom_beeswarm() +
  ylim(0,3)
ggsave(paste0(input_dir, "/Wang_time_revision.pdf"))

# Number of clusters
num_limits = c(0,82)
num_breaks = c(seq(0,80,5), 82)
num_labels = c(seq(0,80,5), ">80")
cap_max = 80
cap_loc = 82

cluster_parameter_list <- 
  data.frame(method = NULL,
             index = NULL,
             n_clusters = NULL)
for (i in c(1:4, 41:47, 55:61, 69:75, 83:171)) {
  method_i <- sub("_[[:digit:]]", "", colnames(seurat_object@meta.data)[i])
  method_i <- sub("[[:digit:]]", "", method_i)
  index_i <- sub(paste0(method_i, "_"), "", colnames(seurat_object@meta.data)[i])
  n_clusters_i <- dplyr::n_distinct(seurat_object@meta.data[,i])
  cluster_parameter_list <- rbind(cluster_parameter_list, data.frame(method = method_i,
                                                                     index = index_i,
                                                                     n_clusters = n_clusters_i))
}

cluster_parameter_list %>%
  mutate(n_clusters_capped = ifelse(n_clusters > cap_max, cap_loc, n_clusters)) %>%
  ggplot(aes(x = method, y = n_clusters_capped, alpha = default)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_beeswarm(shape = 21, fill = "white", cex = 0.25) + 
  geom_point(data = . %>% filter(default == TRUE), shape = 23, alpha = 1, fill = "gold") +
  xlab("Method") +
  ylab("Number of clusters") +
  scale_x_discrete(limits = c("CHOIR", "ArchR", "Signac", "CHOIR_ATAC", "ArchR_ATAC","Signac_ATAC", 
                              "CHOIR_RNA", "CHOIR_RNA_countsplit", "Seurat_RNA"),
                   breaks = c("CHOIR", "ArchR", "Signac", "CHOIR_ATAC", "ArchR_ATAC","Signac_ATAC", 
                              "CHOIR_RNA", "CHOIR_RNA_countsplit", "Seurat_RNA")) +
  NoLegend() +
  scale_alpha_manual(values = c(1,0)) +
  scale_y_continuous(breaks = num_breaks, labels = num_labels, limits = num_limits)
ggsave(paste0(input_dir, "/n_clusters_revision.pdf"), width = 5, height = 4, units = "in")
