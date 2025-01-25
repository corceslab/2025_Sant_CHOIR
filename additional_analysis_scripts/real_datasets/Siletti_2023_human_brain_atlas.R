# ---------------------------------------------------------------------------
# This script contains the analysis and plots for the Siletti et al. 2023
# human brain atlas dataset
# ---------------------------------------------------------------------------
library(CHOIR)
library(Seurat)
library(dplyr)
library(ggplot2)

input_dir <- "/gladstone/corces/lab/users/cpetersen/cluster_benchmarking/Siletti_2023"

# General
object <- readRDS(paste0(input_dir, "/object_combineTrees.rds.rds"))
color_palette <- CHOIRpalette(1302)
object <- runCHOIRumap(object)

object@meta.data$CellID <- rownames(object@meta.data)
object@meta.data$Parent <- all_cluster_ids[rownames(object@meta.data),]$Parent
selected_cols <- object@meta.data %>% dplyr::select(CellID, Batch, Parent, CHOIR_clusters_0.05, supercluster_term, cluster_id, subcluster_id)
write.csv(selected_cols, paste0(input_dir, "/cluster_labels_for_table.csv"))

# Plot
current_label <- "CHOIR_clusters_0.05"
plt <- DimPlot(object, group.by = current_label, reduction = "CHOIR_P0_reduction_UMAP",
               shuffle = TRUE, label = FALSE, raster = FALSE) +
  theme_void() + NoLegend() + ggtitle("") +
  scale_color_manual(values = color_palette) +
  coord_fixed(ratio = 1)
plt[[1]]$layers[[1]]$aes_params$alpha = 0.2
plt[[1]]$layers[[1]]$aes_params$size = 0.05
plt[[1]]$layers[[1]]$aes_params$shape = 16
plt
ggsave(paste0(input_dir, "/", current_label, ".png"),
       height = 5, width = 5, dpi = 2400)

current_label <- "CHOIR_parent"
plt <- DimPlot(seurat_object, group.by = current_label, reduction = "CHOIR_P0_reduction_UMAP",
               shuffle = TRUE, label = TRUE, raster = FALSE) +
  theme_void() + NoLegend() + ggtitle("") +
  scale_color_manual(values = CHOIRpalette(100)) +
  coord_fixed(ratio = 1)
plt[[1]]$layers[[1]]$aes_params$alpha = 0.2
plt[[1]]$layers[[1]]$aes_params$size = 0.05
plt[[1]]$layers[[1]]$aes_params$shape = 16
plt
ggsave(paste0(input_dir, "/", current_label, ".png"),
       height = 5, width = 5, dpi = 2400)

# UMAPs for each parent cluster
for (i in 1:33) {
  print(i)
  current_label <- i
  seurat_object$CHOIR_cluster <- seurat_object$CHOIR_parent == i
  plt <- DimPlot(seurat_object, group.by = "CHOIR_cluster", reduction = "CHOIR_P0_reduction_UMAP",
                 order = TRUE, shuffle = FALSE, label = FALSE, raster = FALSE) +
    theme_void() + NoLegend() + ggtitle("") +
    scale_color_manual(values = c("#DDDDDD", CHOIRpalette(100)[i])) +
    coord_fixed(ratio = 1)
  plt[[1]]$layers[[1]]$aes_params$alpha = 0.2
  plt[[1]]$layers[[1]]$aes_params$size = 0.05
  plt[[1]]$layers[[1]]$aes_params$shape = 16
  plt
  ggsave(paste0(input_dir, "/CHOIR_parent_", i, "_", parent_cluster_key$Parent[i], ".png"),
         height = 5, width = 5, dpi = 2400)
}

# Plot marker genes for parent cluster P24
object <- readRDS(paste0(input_dir, "/seurat_magic_P24.rds"))

parent <- "P24"
choir_data <- readRDS(paste0(input_dir, "/intermediate_files/clusters/Subtrees_final_set/CHOIR_output_", parent, ".rds"))

umap_coords <- RunUMAP(choir_data$reduction$P0_reduction)

color_gradient <- c("#DDDDDD", "#DDDDDD","#FFFC5A",
                    "#FFD32B","#FF9965",
                    "#FD619D","#C732D5","#6E34FC",
                    "#1632FB","#021EA9", "#000436")
color_values <- c(0, 0.08, 0.081, 
                  0.195, 0.31,
                  0.425, 0.54, 0.655,
                  0.77, 0.885, 1)

df <- NULL
gene <- "ENSG00000108821"
ggplot(df) +
  theme_void() +
  geom_point(aes(x = umap_coords@cell.embeddings[,1],
                 y = umap_coords@cell.embeddings[,2],
                 color = object@assays$MAGIC_RNA$data[gene, rownames(umap_coords@cell.embeddings)]), shape = 16,
             size = 0.5, alpha = 0.2) +
  coord_fixed(ratio = 1) +
  scale_color_gradientn(colors = color_gradient,
                        values = color_values) +
  NoLegend()
ggsave(paste0(input_dir, "/P24_", gene, "_magic.png"), width = 5, height = 5, dpi = 2400)

# Min/max
min(object@assays$MAGIC_RNA$data[gene,])
max(object@assays$MAGIC_RNA$data[gene,])

0.08*(max(object@assays$MAGIC_RNA$data[gene,])-
        min(object@assays$MAGIC_RNA$data[gene,])) + 
  min(object@assays$MAGIC_RNA$data[gene,])

# Time & memory usage
compiled_times <- read.csv(paste0(input_dir, "/compiled_times.csv"))

ggplot(compiled_times, aes(x = step, y = time)) +
  theme_classic() +
         geom_point() +
  scale_y_continuous(limits = c(0,130), breaks = c(0,24,48,72,96,120))
ggsave(paste0(input_dir, "/time.pdf"), width = 5, height = 5)

usage_all <- data.frame(total_cpu_time=NULL,
                        maxvmem=NULL,
                        mem=NULL,
                        vmem=NULL,
                        io=NULL)
subtree_set <- c("P1s0", "P1s1", "P1s2", "P2s0", "P2s1", "P2s2", "P2s3", "P3", "P4", "P5", "P6", "P7", "P8", 
                 "P13s0", "P13s1", "P13s2", "P13s3", "P14", "P19", "P20", "P21", "P22", "P23", "P24", "P25")
for (parameter_index in subtree_set) {
  usage_file_i <- list.files(path = paste0(input_dir, "/intermediate_files/clusters/Subtrees_final_set"), pattern = paste0("^usage_parameter_set_", parameter_index))
  usage_i <- read.table(paste0(input_dir, "/intermediate_files/clusters/Subtrees_final_set/", usage_file_i)) %>%
    transmute(total_cpu_time = sub(",", "", sub("cpu=", "", V3)),
              maxvmem = as.numeric(sub("T", "", sub("M", "", sub("G", "", sub("maxvmem=", "", V10))))),
              maxvmem_unit = ifelse(grepl("M", V10), "M", 
                                    ifelse(grepl("G", V10), "G", "T")),
              mem = as.numeric(sub("mem=", "", V4)),
              mem_unit = V5,
              vmem = as.numeric(sub(",", "", sub("M", "", sub("G", "", sub("vmem=", "", V9))))),
              vmem_unit = ifelse(grepl("M", V9), "M", "G"),
              io = as.numeric(sub("io=", "", V7)),
              io_unit = sub(",", "", V8)) %>%
    mutate(mem = ifelse(mem_unit == "MB", mem/1000, mem),
           io = ifelse(io_unit == "MB", io/1000, io),
           maxvmem = ifelse(maxvmem_unit == "M", maxvmem/1000, 
                            ifelse(maxvmem_unit == "T", maxvmem*1000, maxvmem)),
           vmem = ifelse(vmem_unit == "M", vmem/1000, vmem)) %>%
    dplyr::select(-c(mem_unit, io_unit, maxvmem_unit, vmem_unit))
  usage_all <- rbind(usage_all, usage_i)
}
write.csv(usage_all, paste0(input_dir, "/time_usage.csv"))

# Overall confusion matrix
conf_matrix <- table(object@meta.data$subcluster_id,object@meta.data$CHOIR_clusters_0.05)
dim(conf_matrix)
write.csv(conf_matrix, paste0(input_dir, "/confusion_matrix_all.csv"))

# P24 confusion matrix
P24_object <- subset(object, subset = CHOIR_parent == "P24")
con_matrx <- table(P24_object$subcluster_id[rownames(umap_coords@cell.embeddings)], P24_object$CHOIR_clusters_0.05[rownames(umap_coords@cell.embeddings)])
row_sums <- rowSums(con_matrx)
filter <- row_sums > 50
con_matrx <- con_matrx[filter,]
con_matrx <- con_matrx/rowSums(con_matrx)

color_gradient <- c("#FFFFFF", "#000436")
color_values <- c(0, 1)

# First order the columns by the ColSum
con_matrx <- con_matrx[,rev(order(colSums(con_matrx)))]

# Now order the rows by each column value in turn, but only regard values over like 0.2
con_matrx_bi <- con_matrx > 0.5

order_rows <- data.frame(con_matrx_bi) %>% 
  dplyr::arrange(X9, X41, X97, X96, X316, 
                 X586, X523, X468, X601, X103, 
                 X246, X763, X962, X1147, X888, 
                 X1218, X1252, X1296, X1210) %>%
  rownames()
con_matrx <- con_matrx[order_rows,]
superheat::superheat(con_matrx,
                     heat.pal = color_gradient,
                     heat.pal.values = color_values,
                     grid.hline = FALSE,
                     grid.vline = FALSE,
                     smooth.heat = FALSE)
