# ---------------------------------------------------------------------------
# This script generates simulated datasets for clustering method benchmarking 
# using Splatter
# ---------------------------------------------------------------------------
library(splatter)
library(Seurat)
library(scater)
library(dplyr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

# Define args
groups = as.numeric(strsplit(args[1], "_")[[1]][1])
size = as.numeric(strsplit(args[1], "_")[[1]][2])
i = as.numeric(strsplit(args[1], "_")[[1]][3])

# ---------------------------------------------------------------------------
# K = 1
# ---------------------------------------------------------------------------

if (groups == 1) {
  message("[K = 1] Generating dataset ", i, " of size ", c(500, 1000, 1500, 2000, 2500)[size])
  
  # Seed time
  seed_time <- Sys.time()
  
  # Seed
  set.seed(as.integer(i*size + as.numeric(seed_time)))
  
  # Sample parameters randomly from range
  batchCells_param <- c(500, 1000, 1500, 2000, 2500)[size]
  batchCells_param_batch <- rep(batchCells_param/2, 2)
  batch.facLoc_param <- sample(seq(0.05, 0.1, 0.01), 1)
  batch.facScale_param <- sample(seq(0.05, 0.1, 0.01), 1)
  out.prob_param <- sample(seq(0.01, 0.02, 0.001), 1)
  lib.loc_param <- sample(seq(10.5, 11, 0.1), 1)
  lib.scale_param <- sample(seq(0.1, 0.3, 0.05), 1)
  seed_param <- as.integer(i*size + as.numeric(seed_time))
  
  # Start time
  start <- Sys.time()
  
  # Simulate counts
  sim_single <- splatSimulate(nGenes = 10000,
                              batchCells = batchCells_param,
                              out.prob = out.prob_param,
                              lib.loc = lib.loc_param,
                              lib.scale = lib.scale_param,
                              method = "single",
                              seed = seed_param,
                              verbose = FALSE)
  
  # Save raw counts as RDS & as CSV (for PanoView)
  saveRDS(sim_single@assays@data$counts, paste0("/gladstone/mucke/sequencing/Cathrine/cluster_benchmarking/datasets/Splat_1/sim_", batchCells_param, "_", i, "/intermediate_files/preprocessed/RNA_matrix_raw.rds"))
  write.csv(sim_single@assays@data$counts, paste0("/gladstone/mucke/sequencing/Cathrine/cluster_benchmarking/datasets/Splat_1/sim_", batchCells_param, "_", i, "/intermediate_files/preprocessed/RNA_matrix_raw.csv"))
  
  # Write metadata files & save as CSV
  cell_metadata <- data.frame(CellID = sim_single$Cell, Ground_truth = 1)
  write.csv(cell_metadata, paste0("/gladstone/mucke/sequencing/Cathrine/cluster_benchmarking/datasets/Splat_1/sim_", batchCells_param, "_", i, "/intermediate_files/preprocessed/cell_metadata.csv"), row.names = FALSE)
  metadata <- data.frame(variable = c("name", "type", "strategy", "file_type", "tech", "cellranger",
                                      "species", "cells_nuclei", "loaded", "sample_ids", "cell_metadata"),
                         value = c(paste0("Splat_1_", batchCells_param, "_", i),
                                   "RNA_simulated", "ground_truth", "mtx", "RNA", "no",
                                   "human", "cells", "unknown", "sim", "yes"))
  write.csv(metadata, paste0("/gladstone/mucke/sequencing/Cathrine/cluster_benchmarking/datasets/Splat_1/sim_", batchCells_param, "_", i, "/metadata.csv"), row.names = FALSE)
  
  # Log normalization
  sim_single_seurat <- CreateSeuratObject(sim_single@assays@data$counts)
  sim_single_seurat_log <- sim_single_seurat %>%
    NormalizeData()
  # Save log-normalized counts as RDS & as CSV (for PanoView)
  saveRDS(sim_single_seurat_log@assays$RNA@data, paste0("/gladstone/mucke/sequencing/Cathrine/cluster_benchmarking/datasets/Splat_1/sim_", batchCells_param, "_", i, "/intermediate_files/preprocessed/RNA_matrix_norm_log.rds"))
  write.csv(sim_single_seurat_log@assays$RNA@data, paste0("/gladstone/mucke/sequencing/Cathrine/cluster_benchmarking/datasets/Splat_1/sim_", batchCells_param, "_", i, "/intermediate_files/preprocessed/RNA_matrix_norm_log.csv"))
  # Preview UMAP
  sim_single_seurat_log@meta.data$Ground_truth <- 1
  sim_single_seurat_log %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30) %>%
    DimPlot(group.by = "Ground_truth", label = TRUE) + NoLegend()
  ggsave(paste0("/gladstone/mucke/sequencing/Cathrine/cluster_benchmarking/datasets/Splat_1/sim_", batchCells_param, "_", i, "/quality_checks/umap_norm_log.jpg"),
         width = 5, height = 5)
  
  # End time
  end <- Sys.time()
  time_elapsed <- difftime(end, start, units = "secs")
  
  # Save parameter values as CSV
  parameter_df <- data.frame(Parameter = c("nGenes", "batchCells",
                                           "out.prob", "lib.loc", "lib.scale",
                                           "method", "seed", "time_elapsed_sec"),
                             Value = c(10000, batchCells_param, out.prob_param,
                                       lib.loc_param, lib.scale_param,
                                       "single", seed_param, time_elapsed))
  write.csv(parameter_df, paste0("/gladstone/mucke/sequencing/Cathrine/cluster_benchmarking/datasets/Splat_1/sim_", batchCells_param, "_", i, "/quality_checks/parameters.csv"), row.names = FALSE)
}

# ---------------------------------------------------------------------------
# K = 5, 10, or 20
# ---------------------------------------------------------------------------
if (groups > 1) {
  message("[K = ", groups, "] Generating dataset ", i, " of size ", c(5000, 10000, 15000, 20000, 25000)[size])
  
  # Seed time
  seed_time <- Sys.time()
  
  # Seed
  set.seed(as.integer(i*size + as.numeric(seed_time)))
  
  # Sample parameters randomly from range
  batchCells_param <- c(5000, 10000, 15000, 20000, 25000)[size]
  batchCells_param_batch <- rep(batchCells_param/2, 2)
  batch.facLoc_param <- sample(seq(0.025, 0.075, 0.001), 1)
  batch.facScale_param <- sample(seq(0.025, 0.075, 0.001), 1)
  out.prob_param <- sample(seq(0.01, 0.02, 0.001), 1)
  lib.loc_param <- sample(seq(10.5, 11, 0.1), 1)
  lib.scale_param <- sample(seq(0.1, 0.3, 0.05), 1)
  de.prob_param <- sample(seq(0.1, 0.25, 0.001), groups, replace = TRUE)
  de.downProb_param <- sample(seq(0.1, 0.25, 0.001), groups, replace = TRUE)
  de.facLoc_param <- sample(seq(0.05, 0.2, 0.001), groups, replace = TRUE) + 0.0075*((groups/5)-1)
  de.facScale_param <- sample(seq(0.05, 0.2, 0.001), groups, replace = TRUE) + 0.0075*((groups/5)-1)
  
  group.prob_param_pre <- rbeta(groups, 0.75, 0.75)
  # All groups consist of at least 100 cells & consist of at least 2.5% of the data
  for (j in 1:length(group.prob_param_pre)) {
    if ((group.prob_param_pre[j]/sum(group.prob_param_pre)) < max(100/batchCells_param, 0.025)) {
      group.prob_param_pre[j] <- max(((sum(group.prob_param_pre) - group.prob_param_pre[j])*100)/(batchCells_param - 100),
                                     ((sum(group.prob_param_pre) - group.prob_param_pre[j])*0.025)/0.975)
    }
  }
  group.prob_param <- group.prob_param_pre/sum(group.prob_param_pre)
  
  # Boost DE effect of rare groups
  for (j in 1:length(de.prob_param)) {
    if (group.prob_param[j] < 0.05) {
      de.facLoc_param[j] <- de.facLoc_param[j] + 0.05*((groups/5)/size) + 0.05*((groups/5)-1)/size
      de.prob_param[j] <- de.prob_param[j] + 0.05*((groups/5)/size) + 0.05*((groups/5)-1)/size
    }
  }
  
  seed_param <- as.integer(i*size + as.numeric(seed_time))
  
  # Start time
  start <- Sys.time()
  
  # Simulate counts
  sim_groups <- splatSimulate(nGenes = 10000,
                              batchCells = batchCells_param,
                              out.prob = out.prob_param,
                              lib.loc = lib.loc_param,
                              lib.scale = lib.scale_param,
                              group.prob = group.prob_param,
                              de.prob = de.prob_param,
                              de.downProb = de.downProb_param,
                              de.facLoc = de.facLoc_param,
                              de.facScale = de.facScale_param,
                              method = "groups",
                              seed = seed_param,
                              verbose = TRUE)
  print(paste0("Finished simulation ", i))
  
  # Save raw counts as RDS & as CSV (for PanoView)
  saveRDS(sim_groups@assays@data$counts, paste0("/gladstone/mucke/sequencing/Cathrine/cluster_benchmarking/datasets/Splat_", groups, "/sim_", batchCells_param, "_", i, "/intermediate_files/preprocessed/RNA_matrix_raw.rds"))
  write.csv(sim_groups@assays@data$counts, paste0("/gladstone/mucke/sequencing/Cathrine/cluster_benchmarking/datasets/Splat_", groups, "/sim_", batchCells_param, "_", i, "/intermediate_files/preprocessed/RNA_matrix_raw.csv"))
  
  # Write metadata file & save as CSV
  cell_metadata <- data.frame(CellID = sim_groups$Cell, Ground_truth = sim_groups$Group)
  write.csv(cell_metadata, paste0("/gladstone/mucke/sequencing/Cathrine/cluster_benchmarking/datasets/Splat_", groups, "/sim_", batchCells_param, "_", i, "/intermediate_files/preprocessed/cell_metadata.csv"), row.names = FALSE)
  metadata <- data.frame(variable = c("name", "type", "strategy", "file_type", "tech", "cellranger",
                                      "species", "cells_nuclei", "loaded", "sample_ids", "cell_metadata"),
                         value = c(paste0("Splat_", groups, "_", batchCells_param, "_", i),
                                   "RNA_simulated", "ground_truth", "mtx", "RNA", "no",
                                   "human", "cells", "unknown", "sim", "yes"))
  write.csv(metadata, paste0("/gladstone/mucke/sequencing/Cathrine/cluster_benchmarking/datasets/Splat_", groups, "/sim_", batchCells_param, "_", i, "/metadata.csv"), row.names = FALSE)
  
  # Log normalization
  sim_groups_seurat <- CreateSeuratObject(sim_groups@assays@data$counts)
  sim_groups_seurat_log <- sim_groups_seurat %>%
    NormalizeData()
  # Save log-normalized counts as RDS & CSV (required by PanoView)
  saveRDS(sim_groups_seurat_log@assays$RNA@data, paste0("/gladstone/mucke/sequencing/Cathrine/cluster_benchmarking/datasets/Splat_", groups, "/sim_", batchCells_param, "_", i, "/intermediate_files/preprocessed/RNA_matrix_norm_log.rds"))
  write.csv(sim_groups_seurat_log@assays$RNA@data, paste0("/gladstone/mucke/sequencing/Cathrine/cluster_benchmarking/datasets/Splat_", groups, "/sim_", batchCells_param, "_", i, "/intermediate_files/preprocessed/RNA_matrix_norm_log.csv"))
  # Preview UMAP
  sim_groups_seurat_log@meta.data$Ground_truth <- sim_groups$Group
  sim_groups_seurat_log %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30) %>%
    DimPlot(group.by = "Ground_truth", label = TRUE) + NoLegend()
  ggsave(paste0("/gladstone/mucke/sequencing/Cathrine/cluster_benchmarking/datasets/Splat_", groups, "/sim_", batchCells_param, "_", i, "/quality_checks/umap_norm_log.jpg"),
         width = 5, height = 5)
  
  # End time
  end <- Sys.time()
  time_elapsed <- difftime(end, start, units = "secs")
  
  # Save parameter values
  parameter_df <- data.frame(Parameter = c("nGenes", "batchCells",
                                           "out.prob", "lib.loc", "lib.scale",
                                           "group.prob", "de.prob",
                                           "de.downProb", "de.facLoc", "de.facScale",
                                           "method", "seed", "time_elapsed_sec"),
                             Value = c(10000, batchCells_param, 
                                       out.prob_param,
                                       lib.loc_param, lib.scale_param,
                                       paste(group.prob_param, collapse = " "),
                                       paste(de.prob_param, collapse = " "),
                                       paste(de.downProb_param, collapse = " "),
                                       paste(de.facLoc_param, collapse = " "),
                                       paste(de.facScale_param, collapse = " "),
                                       "groups", seed_param, time_elapsed))
  write.csv(parameter_df, paste0("/gladstone/mucke/sequencing/Cathrine/cluster_benchmarking/datasets/Splat_", groups, "/sim_", batchCells_param, "_", i, "/quality_checks/parameters.csv"), row.names = FALSE)
  
}
