# ---------------------------------------------------------------------------
# This script generates simulated datasets with varied parameter modulations 
# using Splatter or scDesign3
# ---------------------------------------------------------------------------
library(splatter)
library(Seurat)
library(scater)
library(dplyr)
library(stringr)
library(scDesign3)

# Define args
groups = 5
size = 5
actual_size = 25000
type_size = "25k"
type <- "batch"

# With batch effects
if (type == "batch") {
  balanced = TRUE
  for (i in 1:5) {
    # Import parameters
    parameter_data <- read.csv(paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", 
                                      groups, "/sim_", actual_size, "_", i, "/quality_checks/parameters.csv"))
    
    message("[K = ", groups, "] Generating dataset ", i, " of size ", c(5000, 10000, 15000, 20000, 25000)[size])
    
    seed_time <- as.numeric(parameter_data$Value[12])
    
    # Seed
    set.seed(seed_time)
    
    # Sample parameters randomly from range
    batchCells_param <- c(5000, 10000, 15000, 20000, 25000)[size]
    if (balanced == TRUE) {
      batchCells_param_batch <- rep(batchCells_param/2, 2)
    } else {
      batchCells_param_batch <- c(round(batchCells_param/3), batchCells_param - round(batchCells_param/3))
    }
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
    
    seed_param <- as.integer(seed_time)
    
    # ---------------------------------------------------------------------------
    # 2 batches
    # ---------------------------------------------------------------------------
    
    # Start time
    start <- Sys.time()
    
    # Simulate counts
    sim_groups <- splatSimulate(nGenes = 10000,
                                batchCells = batchCells_param_batch,
                                out.prob = out.prob_param,
                                lib.loc = lib.loc_param,
                                lib.scale = lib.scale_param,
                                batch.facLoc = batch.facLoc_param,
                                batch.facScale = batch.facScale_param,
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
    saveRDS(sim_groups@assays@data$counts, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_batch/sim_", batchCells_param, "_", i, "/intermediate_files/preprocessed/RNA_matrix_raw.rds"))
    write.csv(sim_groups@assays@data$counts, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_batch/sim_", batchCells_param, "_", i, "/intermediate_files/preprocessed/RNA_matrix_raw.csv"))
    
    # Write metadata file & save as CSV
    cell_metadata <- data.frame(CellID = sim_groups$Cell, Ground_truth = sim_groups$Group, Batch = sim_groups$Batch)
    write.csv(cell_metadata, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_batch/sim_", batchCells_param, "_", i, "/intermediate_files/preprocessed/cell_metadata.csv"), row.names = FALSE)
    metadata <- data.frame(variable = c("name", "type", "strategy", "file_type", "tech", "cellranger",
                                        "species", "cells_nuclei", "loaded", "sample_ids", "cell_metadata"),
                           value = c(paste0("Splatter_", groups, "_", batchCells_param, "_", i),
                                     "RNA_simulated", "ground_truth", "mtx", "RNA", "no",
                                     "human", "cells", "unknown", "sim", "yes"))
    write.csv(metadata, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_batch/sim_", batchCells_param, "_", i, "/metadata.csv"), row.names = FALSE)
    
    # Log normalization
    sim_groups_seurat <- CreateSeuratObject(sim_groups@assays@data$counts)
    sim_groups_seurat_log <- sim_groups_seurat %>%
      NormalizeData()
    # Save log-normalized counts as RDS & CSV (required by PanoView)
    saveRDS(sim_groups_seurat_log@assays$RNA$data, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_batch/sim_", batchCells_param, "_", i, "/intermediate_files/preprocessed/RNA_matrix_norm_log.rds"))
    write.csv(sim_groups_seurat_log@assays$RNA$data, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_batch/sim_", batchCells_param, "_", i, "/intermediate_files/preprocessed/RNA_matrix_norm_log.csv"))
    # Preview UMAP
    sim_groups_seurat_log@meta.data$Ground_truth <- sim_groups$Group
    sim_groups_seurat_log %>%
      FindVariableFeatures() %>%
      ScaleData() %>%
      RunPCA() %>%
      RunUMAP(dims = 1:30) %>%
      DimPlot(group.by = "Ground_truth", label = TRUE) + NoLegend()
    ggsave(paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_batch/sim_", batchCells_param, "_", i, "/quality_checks/umap_norm_log.jpg"),
           width = 5, height = 5)
    
    # End time
    end <- Sys.time()
    time_elapsed <- difftime(end, start, units = "secs")
    
    # Save parameter values
    parameter_df <- data.frame(Parameter = c("nGenes", "batchCells",
                                             "out.prob", "lib.loc", "lib.scale",
                                             "batch.loc", "batch.scale",
                                             "group.prob", "de.prob",
                                             "de.downProb", "de.facLoc", "de.facScale",
                                             "method", "seed", "time_elapsed_sec"),
                               Value = c(10000, batchCells_param, 
                                         out.prob_param,
                                         lib.loc_param, lib.scale_param,
                                         batch.facLoc_param, batch.facScale_param,
                                         paste(group.prob_param, collapse = " "),
                                         paste(de.prob_param, collapse = " "),
                                         paste(de.downProb_param, collapse = " "),
                                         paste(de.facLoc_param, collapse = " "),
                                         paste(de.facScale_param, collapse = " "),
                                         "groups", seed_param, time_elapsed))
    write.csv(parameter_df, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_batch/sim_", batchCells_param, "_", i, "/quality_checks/parameters.csv"), row.names = FALSE)
    
  }
}

# With DE modulation
if (type == "DE") {
  for (i in c(1:5)) {
    for (d in c(0.0125, 0.025, 0.05, 0.1, 0.2, 0.4)) { 
      
      # Import parameters
      parameter_data <- read.csv(paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", 
                                        groups, "/sim_", actual_size, "_", i, "/quality_checks/parameters.csv"))
      
      message("[K = ", groups, "] Generating dataset ", i, " of size ", c(5000, 10000, 15000, 20000, 25000)[size])
      
      seed_time <- as.numeric(parameter_data$Value[12])
      
      # Seed
      set.seed(seed_time)
      
      # Sample parameters randomly from range
      batchCells_param <- c(5000, 10000, 15000, 20000, 25000)[size]
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
      
      seed_param <- as.integer(seed_time)
      
      # Set de.prob and de.downprob of group 1 & 2
      de.prob_param[1] <- d
      de.downProb_param[1] <- d
      de.prob_param[2] <- d
      de.downProb_param[2] <- d
      
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
      print(paste0("Finished simulation ", d))
      
      # Save raw counts as RDS & as CSV (for PanoView)
      saveRDS(sim_groups@assays@data$counts, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_DE/sim_", batchCells_param, "_", i, "_", d, "/intermediate_files/preprocessed/RNA_matrix_raw.rds"))
      write.csv(sim_groups@assays@data$counts, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_DE/sim_", batchCells_param, "_", i, "_", d, "/intermediate_files/preprocessed/RNA_matrix_raw.csv"))
      
      # Write metadata file & save as CSV
      cell_metadata <- data.frame(CellID = sim_groups$Cell, Ground_truth = sim_groups$Group, Batch = sim_groups$Batch)
      write.csv(cell_metadata, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_DE/sim_", batchCells_param, "_", i, "_", d, "/intermediate_files/preprocessed/cell_metadata.csv"), row.names = FALSE)
      metadata <- data.frame(variable = c("name", "type", "strategy", "file_type", "tech", "cellranger",
                                          "species", "cells_nuclei", "loaded", "sample_ids", "cell_metadata"),
                             value = c(paste0("Splat_", groups, "_", batchCells_param, "_", i),
                                       type_size, "ground_truth", "mtx", "RNA", "no",
                                       "human", "cells", "unknown", "sim", "yes"))
      write.csv(metadata, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_DE/sim_", batchCells_param, "_", i, "_", d, "/metadata.csv"), row.names = FALSE)
      
      # Log normalization
      sim_groups_seurat <- CreateSeuratObject(sim_groups@assays@data$counts)
      sim_groups_seurat_log <- sim_groups_seurat %>%
        NormalizeData()
      # Save log-normalized counts as RDS & CSV (required by PanoView)
      saveRDS(sim_groups_seurat_log@assays$RNA$data, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_DE/sim_", batchCells_param, "_", i, "_", d, "/intermediate_files/preprocessed/RNA_matrix_norm_log.rds"))
      write.csv(sim_groups_seurat_log@assays$RNA$data, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_DE/sim_", batchCells_param, "_", i, "_", d, "/intermediate_files/preprocessed/RNA_matrix_norm_log.csv"))
      # Preview UMAP
      sim_groups_seurat_log@meta.data$Ground_truth <- sim_groups$Group
      sim_groups_seurat_log %>%
        FindVariableFeatures() %>%
        ScaleData() %>%
        RunPCA() %>%
        RunUMAP(dims = 1:30) %>%
        DimPlot(group.by = "Ground_truth", label = TRUE) + NoLegend()
      ggsave(paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_DE/sim_", batchCells_param, "_", i, "_", d, "/quality_checks/umap_norm_log.jpg"),
             width = 5, height = 5)
      
      # End time
      end <- Sys.time()
      time_elapsed <- difftime(end, start, units = "secs")
      
      # Save parameter values
      parameter_df <- data.frame(Parameter = c("nGenes", "batchCells",
                                               "out.prob", "lib.loc", "lib.scale",
                                               "batch.loc", "batch.scale",
                                               "group.prob", "de.prob",
                                               "de.downProb", "de.facLoc", "de.facScale",
                                               "method", "seed", "time_elapsed_sec"),
                                 Value = c(10000, batchCells_param, 
                                           out.prob_param,
                                           lib.loc_param, lib.scale_param,
                                           batch.facLoc_param, batch.facScale_param,
                                           paste(group.prob_param, collapse = " "),
                                           paste(de.prob_param, collapse = " "),
                                           paste(de.downProb_param, collapse = " "),
                                           paste(de.facLoc_param, collapse = " "),
                                           paste(de.facScale_param, collapse = " "),
                                           "groups", seed_param, time_elapsed))
      write.csv(parameter_df, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_DE/sim_", batchCells_param, "_", i, "_", d, "/quality_checks/parameters.csv"), row.names = FALSE)
    }
  }
}

# With library size modulation
if (type == "library_size") {
  for (i in c(1,2,3,4,5)) {
    for (d in c(11, 10.8, 10.6, 10.4, 10.2, 10, 9.8, 9.6, 9.4, 9.2)) { 
      print(d)
      # Import parameters
      parameter_data <- read.csv(paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", 
                                        groups, "/sim_", actual_size, "_", i, "/quality_checks/parameters.csv"))
      message("[K = ", groups, "] Generating dataset ", i, " of size ", c(5000, 10000, 15000, 20000, 25000)[size])
      
      seed_time <- as.numeric(parameter_data$Value[12])
      
      # Seed
      set.seed(seed_time)
      
      # Sample parameters randomly from range
      batchCells_param <- c(5000, 10000, 15000, 20000, 25000)[size]
      batch.facLoc_param <- sample(seq(0.025, 0.075, 0.001), 1)
      batch.facScale_param <- sample(seq(0.025, 0.075, 0.001), 1)
      out.prob_param <- as.numeric(parameter_data$Value[3])
      lib.loc_param <- d
      lib.scale_param <- as.numeric(parameter_data$Value[5])
      de.prob_param <- as.numeric(strsplit(parameter_data$Value[7], " ")[[1]])
      de.downProb_param <- as.numeric(strsplit(parameter_data$Value[8], " ")[[1]])
      de.facLoc_param <- as.numeric(strsplit(parameter_data$Value[9], " ")[[1]])
      de.facScale_param <- as.numeric(strsplit(parameter_data$Value[10], " ")[[1]])
      
      group.prob_param <- as.numeric(strsplit(parameter_data$Value[6], " ")[[1]])
      
      seed_param <- as.integer(seed_time)
      
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
      print(paste0("Finished simulation ", d))
      
      # Save raw counts as RDS & as CSV (for PanoView)
      saveRDS(sim_groups@assays@data$counts, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_library_size/sim_", batchCells_param, "_", i, "_", d, "/intermediate_files/preprocessed/RNA_matrix_raw.rds"))
      write.csv(sim_groups@assays@data$counts, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_library_size/sim_", batchCells_param, "_", i, "_", d, "/intermediate_files/preprocessed/RNA_matrix_raw.csv"))
      
      # Write metadata file & save as CSV
      cell_metadata <- data.frame(CellID = sim_groups$Cell, Ground_truth = sim_groups$Group, Batch = sim_groups$Batch)
      write.csv(cell_metadata, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_library_size/sim_", batchCells_param, "_", i, "_", d, "/intermediate_files/preprocessed/cell_metadata.csv"), row.names = FALSE)
      metadata <- data.frame(variable = c("name", "type", "strategy", "file_type", "tech", "cellranger",
                                          "species", "cells_nuclei", "loaded", "sample_ids", "cell_metadata"),
                             value = c(paste0("Splat_", groups, "_", batchCells_param, "_", i),
                                       type_size, "ground_truth", "mtx", "RNA", "no",
                                       "human", "cells", "unknown", "sim", "yes"))
      write.csv(metadata, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_library_size/sim_", batchCells_param, "_", i, "_", d, "/metadata.csv"), row.names = FALSE)
      
      # Log normalization
      sim_groups_seurat <- CreateSeuratObject(sim_groups@assays@data$counts)
      sim_groups_seurat_log <- sim_groups_seurat %>%
        NormalizeData()
      # Save log-normalized counts as RDS & CSV (required by PanoView)
      saveRDS(sim_groups_seurat_log@assays$RNA$data, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_library_size/sim_", batchCells_param, "_", i, "_", d, "/intermediate_files/preprocessed/RNA_matrix_norm_log.rds"))
      write.csv(sim_groups_seurat_log@assays$RNA$data, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_library_size/sim_", batchCells_param, "_", i, "_", d, "/intermediate_files/preprocessed/RNA_matrix_norm_log.csv"))
      # Preview UMAP
      sim_groups_seurat_log@meta.data$Ground_truth <- sim_groups$Group
      sim_groups_seurat_log %>%
        FindVariableFeatures() %>%
        ScaleData() %>%
        RunPCA() %>%
        RunUMAP(dims = 1:30) %>%
        DimPlot(group.by = "Ground_truth", label = TRUE) + NoLegend()
      ggsave(paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_library_size/sim_", batchCells_param, "_", i, "_", d, "/quality_checks/umap_norm_log.jpg"),
             width = 5, height = 5)
      
      # End time
      end <- Sys.time()
      time_elapsed <- difftime(end, start, units = "secs")
      
      # Save parameter values
      parameter_df <- data.frame(Parameter = c("nGenes", "batchCells",
                                               "out.prob", "lib.loc", "lib.scale",
                                               "batch.loc", "batch.scale",
                                               "group.prob", "de.prob",
                                               "de.downProb", "de.facLoc", "de.facScale",
                                               "method", "seed", "time_elapsed_sec"),
                                 Value = c(10000, batchCells_param, 
                                           out.prob_param,
                                           lib.loc_param, lib.scale_param,
                                           batch.facLoc_param, batch.facScale_param,
                                           paste(group.prob_param, collapse = " "),
                                           paste(de.prob_param, collapse = " "),
                                           paste(de.downProb_param, collapse = " "),
                                           paste(de.facLoc_param, collapse = " "),
                                           paste(de.facScale_param, collapse = " "),
                                           "groups", seed_param, time_elapsed))
      write.csv(parameter_df, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_library_size/sim_", batchCells_param, "_", i, "_", d, "/quality_checks/parameters.csv"), row.names = FALSE)
    }
  }
}

# Downsampled
if (type == "downsample") {
  # Set up values
  n_cells <- actual_size
  type <- type_size
  n_groups <- groups
  iteration <- 1
  n_downsample <- c(500, 400, 300, 200, 100)
  key <- n_downsample
  seed <- 44226
  set.seed(seed)
  
  input_dir <- paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", 
                      n_groups, "/sim_", n_cells, "_", iteration)
  output_dir <- paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", 
                       n_groups, "_", n_cells, "_downsample2")
  
  # Import data
  cell_metadata <- read.csv(paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata.csv"))
  count_matrix <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/RNA_matrix_raw.rds"))
  
  # Find smallest group truth group
  smallest_group <- names(table(cell_metadata$Ground_truth))[which(table(cell_metadata$Ground_truth) == 
                                                                     min(table(cell_metadata$Ground_truth)))]
  
  # Get cell IDs
  smallest_group_cell_ids <- dplyr::filter(cell_metadata, Ground_truth == smallest_group)$CellID
  other_group_cell_ids <- dplyr::filter(cell_metadata, Ground_truth != smallest_group)$CellID
  
  # Downsample iteratively
  cell_ids_list <- vector(mode = "list", length = length(n_downsample))
  cell_ids_list[[1]] <- sample(smallest_group_cell_ids, n_downsample[1])
  for (i in 2:length(n_downsample)) {
    cell_ids_list[[i]] <- sample(cell_ids_list[[i-1]], n_downsample[i])
  }
  
  # Write new files
  for (i in 1:length(cell_ids_list)) {
    message(Sys.time(), " : Downsampling ", i, "..")
    cell_ids_current <- cell_ids_list[[i]]
    cell_metadata_current <- cell_metadata %>% 
      dplyr::filter(CellID %in% c(other_group_cell_ids, cell_ids_current))
    count_matrix_current <- count_matrix[, cell_metadata_current$CellID]
    norm_log_matrix_current <- NormalizeData(count_matrix_current)
    # Write metadata files & save as CSV
    write.csv(cell_metadata_current, paste0(output_dir, "/sim_", n_cells, "_", iteration, "_", key[i], 
                                            "/intermediate_files/preprocessed/cell_metadata.csv"), row.names = FALSE)
    metadata <- data.frame(variable = c("name", "type", "strategy", "file_type", "tech", "cellranger",
                                        "species", "cells_nuclei", "loaded", "sample_ids", "cell_metadata", "seed"),
                           value = c(paste0("Splatter_", n_groups, "_", n_cells, "_", iteration, "_", key[i]),
                                     type, "ground_truth", "mtx", "RNA", "no",
                                     "human", "cells", "unknown", "sim", "yes", seed))
    write.csv(metadata, paste0(output_dir, "/sim_", n_cells, "_", iteration, "_", key[i], "/metadata.csv"), row.names = FALSE)
    # Save raw counts as RDS & as CSV (for PanoView)
    saveRDS(count_matrix_current, paste0(output_dir, "/sim_", n_cells, "_", iteration, "_", key[i], 
                                         "/intermediate_files/preprocessed/RNA_matrix_raw.rds"))
    write.csv(count_matrix_current, paste0(output_dir, "/sim_", n_cells, "_", iteration, "_", key[i], 
                                           "/intermediate_files/preprocessed/RNA_matrix_raw.csv"))
    # Save normalized counts as RDS & as CSV (for PanoView)
    saveRDS(norm_log_matrix_current, paste0(output_dir, "/sim_", n_cells, "_", iteration, "_", key[i], 
                                            "/intermediate_files/preprocessed/RNA_matrix_norm_log.rds"))
    write.csv(norm_log_matrix_current, paste0(output_dir, "/sim_", n_cells, "_", iteration, "_", key[i], 
                                              "/intermediate_files/preprocessed/RNA_matrix_norm_log.csv"))
  }
}

# Trajectory
if (type == "trajectory") {
  
  message("[K = ", groups, "] Generating dataset ", i, " of size ", c(5000, 10000, 15000, 20000, 25000)[size])
  
  # Seed time
  seed_time <- Sys.time()
  
  # Seed
  set.seed(as.integer(i*size + as.numeric(seed_time)))
  
  # Sample parameters randomly from range
  batchCells_param <- c(5000, 10000, 15000, 20000, 25000)[size]
  out.prob_param <- sample(seq(0.01, 0.02, 0.001), 1)
  lib.loc_param <- sample(seq(10.5, 11, 0.1), 1)
  lib.scale_param <- sample(seq(0.1, 0.3, 0.05), 1)
  de.prob_param <- sample(seq(0.1, 0.25, 0.001), groups, replace = TRUE)
  de.downProb_param <- sample(seq(0.1, 0.25, 0.001), groups, replace = TRUE)
  de.facLoc_param <- sample(seq(0.05, 0.2, 0.001), groups, replace = TRUE) + 0.0075*((groups/5)-1)
  de.facScale_param <- sample(seq(0.05, 0.2, 0.001), groups, replace = TRUE) + 0.0075*((groups/5)-1)
  
  #origin <- c(0,0, sample(c(0,1,2), groups-2, replace = TRUE))
  origin <- c(0,1,2,3,4)
  steps <- round(runif (groups, min = 50, max = 100))
  
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
  
  seed_param <- as.integer(seed_time)
  
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
                              method = "paths",
                              path.from = origin,
                              path.skew = 0,
                              path.nSteps = steps,
                              seed = seed_param,
                              verbose = TRUE)
  print(paste0("Finished simulation ", i))
  
  # Log normalization
  sim_groups_seurat <- CreateSeuratObject(sim_groups@assays@data$counts)
  sim_groups_seurat_log <- sim_groups_seurat %>%
    NormalizeData()
  # Preview UMAP
  sim_groups_seurat_log@meta.data$Ground_truth <- sim_groups$Group
  sim_groups_seurat_log %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30) %>%
    DimPlot(group.by = "Ground_truth", label = TRUE) + NoLegend()
  ggsave(paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_trajectory/sim_", batchCells_param, "_", i, "/quality_checks/umap_norm_log.jpg"),
         width = 5, height = 5)
  
  # Save raw counts as RDS & as CSV (for PanoView)
  saveRDS(sim_groups@assays@data$counts, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_trajectory/sim_", batchCells_param, "_", i, "/intermediate_files/preprocessed/RNA_matrix_raw.rds"))
  write.csv(sim_groups@assays@data$counts, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_trajectory/sim_", batchCells_param, "_", i, "/intermediate_files/preprocessed/RNA_matrix_raw.csv"))
  
  # Write metadata file & save as CSV
  cell_metadata <- data.frame(CellID = sim_groups$Cell, Ground_truth = sim_groups$Group, Batch = sim_groups$Batch)
  write.csv(cell_metadata, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_trajectory/sim_", batchCells_param, "_", i, "/intermediate_files/preprocessed/cell_metadata.csv"), row.names = FALSE)
  metadata <- data.frame(variable = c("name", "type", "strategy", "file_type", "tech", "cellranger",
                                      "species", "cells_nuclei", "loaded", "sample_ids", "cell_metadata"),
                         value = c(paste0("Splat_", groups, "_", batchCells_param, "_", i),
                                   type_size, "ground_truth", "mtx", "RNA", "no",
                                   "human", "cells", "unknown", "sim", "yes"))
  write.csv(metadata, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_trajectory/sim_", batchCells_param, "_", i, "/metadata.csv"), row.names = FALSE)
  
  # Save log-normalized counts as RDS & CSV (required by PanoView)
  saveRDS(sim_groups_seurat_log@assays$RNA$data, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_trajectory/sim_", batchCells_param, "_", i, "/intermediate_files/preprocessed/RNA_matrix_norm_log.rds"))
  write.csv(sim_groups_seurat_log@assays$RNA$data, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_trajectory/sim_", batchCells_param, "_", i, "/intermediate_files/preprocessed/RNA_matrix_norm_log.csv"))
  
  # End time
  end <- Sys.time()
  time_elapsed <- difftime(end, start, units = "secs")
  
  # Save parameter values
  parameter_df <- data.frame(Parameter = c("nGenes", "batchCells",
                                           "out.prob", "lib.loc", "lib.scale",
                                           "group.prob", "de.prob",
                                           "de.downProb", "de.facLoc", "de.facScale",
                                           "method", "seed", "time_elapsed_sec", "steps", "origin"),
                             Value = c(10000, batchCells_param, 
                                       out.prob_param,
                                       lib.loc_param, lib.scale_param,
                                       paste(group.prob_param, collapse = " "),
                                       paste(de.prob_param, collapse = " "),
                                       paste(de.downProb_param, collapse = " "),
                                       paste(de.facLoc_param, collapse = " "),
                                       paste(de.facScale_param, collapse = " "),
                                       "groups", seed_param, time_elapsed,
                                       paste(steps, collapse = " "),
                                       paste(origin, collapse = " ")))
  write.csv(parameter_df, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "_trajectory/sim_", batchCells_param, "_", i, "/quality_checks/parameters.csv"), row.names = FALSE)
  
}

# scDesign3
if (type == "scDesign3") {
  key <- "sce_full_Zhengmix4uneq"
  sim_key <- "sim_2500_1"
  ncells <- actual_size
  
  sce <- get(key)(metadata = FALSE)
  colData(sce)$cell_type = as.factor(colData(sce)$phenoid)
  colData(sce)$library = colSums(counts(sce))
  
  # Simulate with scDesign3
  seed <- as.integer(Sys.time())
  set.seed(seed)
  start_time <- Sys.time()
  example_simu <- scdesign3(
    sce = sce,
    assay_use = "counts",
    celltype = "cell_type",
    pseudotime = NULL,
    spatial = NULL,
    other_covariates = "library",
    mu_formula = "cell_type + offset(log(library))",
    sigma_formula = "1",
    family_use = "nb",
    n_cores = 2,
    usebam = FALSE,
    corr_formula = "1",
    copula = "gaussian",
    DT = TRUE,
    pseudo_obs = FALSE,
    return_model = FALSE,
    nonzerovar = FALSE,
    parallelization = "pbmcmapply", ncell = ncells) # ) #
  end_time <- Sys.time()
  end_time-start_time
  # Preview UMAP
  count_matrix_sim <- example_simu$new_count
  seurat_sim <- CreateSeuratObject(count_matrix_sim)
  seurat_sim <- seurat_sim %>%
    NormalizeData()
  seurat_sim@meta.data$Ground_truth <- example_simu$new_covariate$cell_type
  seurat_sim %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30) %>%
    DimPlot(group.by = "Ground_truth", label = TRUE) + NoLegend()
  ggsave(paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Zheng_2017/", 
                sim_key, "/quality_checks/umap_norm_log.jpg"),
         width = 5, height = 5)
  
  # Save simulated data
  saveRDS(seurat_sim@assays$RNA$counts, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Zheng_2017/", 
                                               sim_key, "/intermediate_files/preprocessed/RNA_matrix_raw.rds"))
  write.csv(seurat_sim@assays$RNA$counts, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Zheng_2017/", 
                                                 sim_key, "/intermediate_files/preprocessed/RNA_matrix_raw.csv"))
  saveRDS(seurat_sim@assays$RNA$data, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Zheng_2017/", 
                                             sim_key, "/intermediate_files/preprocessed/RNA_matrix_norm_log.rds"))
  write.csv(seurat_sim@assays$RNA$data, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Zheng_2017/", 
                                               sim_key, "/intermediate_files/preprocessed/RNA_matrix_norm_log.csv"))
  
  cell_metadata_sim <- data.frame(CellID = colnames(count_matrix_sim), Ground_truth = example_simu$new_covariate$cell_type)
  write.csv(cell_metadata_sim, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Zheng_2017/", 
                                      sim_key, "/intermediate_files/preprocessed/cell_metadata.csv"), row.names = FALSE)
  metadata_sim <- data.frame(variable = c("name", "type", "strategy", "file_type", "tech", "cellranger",
                                          "species", "cells_nuclei", "loaded", "sample_ids", "cell_metadata"),
                             value = c(sim_key,
                                       "5k", "ground_truth", "mtx", "RNA", "no",
                                       "human", "cells", "unknown", "sim", "yes"))
  write.csv(metadata_sim, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Zheng_2017/", 
                                 sim_key, "/metadata.csv"), row.names = FALSE)
  
  # Save parameter values
  parameter_df <- data.frame(Parameter = c("assay_use",
                                           "celltype",
                                           "pseudotime",
                                           "spatial",
                                           "other_covariates",
                                           "mu_formula",
                                           "sigma_formula",
                                           "family_use",
                                           "n_cores",
                                           "usebam",
                                           "corr_formula",
                                           "copula",
                                           "DT",
                                           "pseudo_obs",
                                           "return_model",
                                           "nonzerovar",
                                           "parallelization",
                                           "seed", "ncell"), # ), #
                             Value = c("counts",
                                       "cell_type",
                                       "NULL",
                                       "NULL",
                                       "library",
                                       "cell_type + offset(log(library))",
                                       "1",
                                       "nb",
                                       2,
                                       FALSE,
                                       "1",
                                       "gaussian",
                                       TRUE,
                                       FALSE,
                                       FALSE,
                                       FALSE,
                                       "pbmcmapply",
                                       seed, ncells)) # )) #
  write.csv(parameter_df, paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Zheng_2017/", 
                                 sim_key, "/quality_checks/parameters.csv"), row.names = FALSE)
}