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

# Set color palette
color_palette <- ClusteringBetaTest::CHOIRpalette(100)

# Import data
input_dir <- "/gladstone/corces/lab/users/cpetersen/cluster_benchmarking/Wang_2022_multiome"

raw_matrix <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/RNA_matrix_raw.rds"))
raw_matrix <- Azimuth:::ConvertEnsembleToSymbol(mat = raw_matrix, species = "human")

# Create Seurat Object
options(Seurat.object.assay.version = 'v5')
seurat_object <- CreateSeuratObject(counts = raw_matrix)

# Add metadata
cell_metadata <- read.csv(paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata.csv"))
seurat_object <- subset(seurat_object, cells = cell_metadata$CellID)
rownames(cell_metadata) <- cell_metadata$CellID
cell_metadata <- cell_metadata[colnames(seurat_object),]
seurat_object@meta.data <- cbind(seurat_object@meta.data, cell_metadata)

# Normalize & find variable features
seurat_object <- NormalizeData(seurat_object)

# Import CHOIR UMAP coordinates
parameter_list <- "CHOIR_parameter_list_multiome_batch"
choir_data <- readRDS(paste0(input_dir, "/intermediate_files/clusters/", parameter_list, "/CHOIR_output_1.rds"))
umap_coords <- RunUMAP(choir_data$reduction$P0_reduction)
# Add UMAP coordinates to Seurat object
umap_coords@cell.embeddings <- umap_coords@cell.embeddings[colnames(seurat_object),]
seurat_object[["P0_umap"]] <- CreateDimReducObject(embeddings = umap_coords@cell.embeddings,
                                                   assay = "RNA")

# Import clusters
# CHOIR
parameter_list <- "CHOIR_parameter_list_multiome_batch"
CHOIR_clusters <- read.csv(paste0(input_dir, 
                                  "/intermediate_files/clusters/compiled_clusters_", parameter_list, ".csv"))
rownames(CHOIR_clusters) <- CHOIR_clusters$CellID
colnames(CHOIR_clusters) <- sub("Parameters", "CHOIR", colnames(CHOIR_clusters))
CHOIR_clusters <- CHOIR_clusters[seurat_object@meta.data$CellID,]

# Add clustering results to Seurat object
seurat_object@meta.data$CHOIR_1 <- CHOIR_clusters[,-1]

# UMAPs for clusters
current_label <- "CHOIR_1"
plt <- DimPlot(seurat_object, group.by = current_label, reduction = "P0_umap",
               shuffle = TRUE, label = TRUE, raster = FALSE) +
  theme_void() + NoLegend() + ggtitle(current_label) +
  scale_color_manual(values = color_palette) +
  coord_fixed(ratio = 1) +
  xlim(-18,18) + ylim(-18,18)
plt[[1]]$layers[[1]]$aes_params$alpha = 0.2
plt[[1]]$layers[[1]]$aes_params$size = 0.2
plt[[1]]$layers[[1]]$aes_params$shape = 16
plt

# Plot with accuracy scores
plotCHOIR(seurat_object, accuracy_scores = TRUE, reduction = "P0_umap")

# Feature importance values for rods vs. cones
vals <- dplyr::filter(choir_data$records$feature_importance_records, 
              Cluster1 == "P3_L3_1")[,3:(ncol(choir_data$records$feature_importance_records)-2)] %>% c() %>% unlist()
# Top 10
sort(vals, decreasing=TRUE)[1:10]

# Log fold change gene expression for rods vs. cones
mks <- FindMarkers(seurat_object, ident.1 = 1, ident.2 = 4, group.by = "CHOIR_1")
mks$name <- rownames(mks)

# Comparison of feature importances vs. log fold change of gene expression
gene_key <- choir_data$var_features$P0_var_features[[2]] %>% data.frame()
gene_key$gene_seq_id <- paste0(gene_key$seqnames, " ", gene_key$start, " _")
vals_rna <- vals[names(vals) %in% gene_key$gene_seq_id]
vals_rna_df <- data.frame(gene_seq_id = names(vals_rna),
                          feature_importance = vals_rna)
vals_rna_df_merged <- merge(vals_rna_df, gene_key, by = "gene_seq_id", all.x = TRUE)
vals_rna_df_merged_mks <- merge(vals_rna_df_merged, mks, by = "name", all = TRUE)
# Top features
vals_rna_df_merged_mks %>% arrange(-feature_importance) %>% head()

# Plot
vals_rna_df_merged_mks %>%
  ggplot(aes(x = feature_importance, y = abs(avg_log2FC))) +
  geom_point(shape = 21, fill = "white") +
  theme_classic()

# Filter feature importances to ATAC features only
atac_vals <- vals[grepl("00 _", names(vals))]
atac_vals <- atac_vals[!grepl("100 _", names(atac_vals))]
atac_vals <- atac_vals[!grepl("800 _", names(atac_vals))]
sort(atac_vals, decreasing=TRUE)[1:10]

# ArchR Project
archr_proj <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/ATAC_obj_ArchR.rds"))
setwd(paste0(input_dir, "/intermediate_files/arrow_files"))
archr_proj_sub <- subsetArchRProject(archr_proj, cells = colnames(seurat_object), force = TRUE)
archr_proj_sub@cellColData$CHOIR_1 <- paste0("C", seurat_object@meta.data[rownames(archr_proj_sub@cellColData), "CHOIR_1"])

# Plot browser track for most importance ATAC features
chrom <- "chr1"
start_val <- 228457500
gr <- makeGRangesFromDataFrame(data.frame(chr=chrom, start=start_val - 20000, end=start_val + 20500))
gr_highlight <- makeGRangesFromDataFrame(data.frame(chr=chrom, start=start_val, end=start_val + 500))
p <- plotBrowserTrack(archr_proj_sub, region = gr, groupBy = "CHOIR_1",
                      highlight = gr_highlight)
grid::grid.newpage()
grid::grid.draw(p)

# MAGIC-imputed feature plots
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
gene <- "ENSG00000144278"
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

min(magic_object@assays$MAGIC_RNA$data[gene,])
max(magic_object@assays$MAGIC_RNA$data[gene,])

# Clustering tree
clustree(choir_data$clusters$full_tree, prefix = "L")

# Computational time
CHOIR_time <- read.csv(paste0(input_dir, "/intermediate_files/clusters/compiled_time_CHOIR_parameter_list_multiome_batch.csv"))
CHOIR_time$method <- "CHOIR"
CHOIR_time %>%
  ggplot(aes(x = method, y = time_secs/3600, color = default)) +
  geom_beeswarm() +
  ylim(0,1.5)
