# ---------------------------------------------------------------------------
# Run CHOIR combineTrees on Siletti et al. 2023 human brain atlas dataset
# ---------------------------------------------------------------------------
library(CHOIR)
library(Seurat)
library(dplyr)

start_time <- Sys.time()

input_dir <- "/gladstone/corces/lab/users/cpetersen/cluster_benchmarking/Siletti_2023"

subtree_set <- c("P1s0", "P1s1", "P1s2", "P2s0", "P2s1", "P2s2", "P2s3", "P3", "P4", "P5", "P6", "P7", "P8", 
                 "P13s0", "P13s1", "P13s2", "P13s3", "P14", "P19", "P20", "P21", "P22", "P23", "P24", "P25")
other_parents_set <- c(9, 10, 11, 12, 15, 16, 17, 18)

# Import object
seurat_object <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/object_parent_clusters.rds"))

# Import & compile subtree records
subtree_records_list <- vector(mode = "list", length = 33)
for (s in 1:length(subtree_set)) {
  subtree_records_list[[s]] <- readRDS(paste0(input_dir, "/intermediate_files/clusters/", group, "/CHOIR_output_", subtree_set[s], ".rds"))
}

# For parent clusters that could not be subclustered, create record structure
for (s in 1:length(other_parents_set)) {
  current_cluster <- other_parents_set[s]
  cell_ids <- rownames(seurat_object@misc$CHOIR$clusters$P0_tree)[seurat_object@misc$CHOIR$clusters$P0_tree[, ncol(seurat_object@misc$CHOIR$clusters$P0_tree)] == current_cluster]
  subtree_x <- subset(seurat_object, subset = CellID %in% cell_ids)
  # Add cluster labels for a single cluster
  subtree_x@misc$CHOIR_subtree$clusters$CHOIR_clusters_0.05 <- data.frame(CellID = cell_ids,
                                                                          CHOIR_clusters_0.05 = 1,
                                                                          Record_cluster_label = "P0_L0_1")
  # Add parameter records matching remaining subtrees
  subtree_x@misc$CHOIR_subtree$parameters <- subtree_records_list[[1]]parameters
  subtree_records_list[[25+s]] <- subtree_x@misc$CHOIR_subtree
}

# Combine trees
seurat_object <- combineTrees(seurat_object,
                              subtree_list = subtree_records_list)

end_time <- Sys.time()
diff_time <- difftime(end_time, start_time, units = "secs")

saveRDS(seurat_object@misc$CHOIR, paste0(input_dir, "/intermediate_files/preprocessed/CHOIR_object_combineTrees.rds"))
saveRDS(seurat_object, paste0(input_dir, "/intermediate_files/preprocessed/object_combineTrees.rds"))
saveRDS(diff_time, paste0(input_dir, "/combineTrees_difftime.rds"))
