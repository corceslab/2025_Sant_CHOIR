# ---------------------------------------------------------------------------
# Run Seurat clustering
# ---------------------------------------------------------------------------
library(Seurat)
library(BPCells)
library(SeuratObject)
library(harmony)

run_Seurat <- function(input_dir,
                       cluster_parameter_index,
                       cluster_parameters) {
  
  # Set up ---------------------------
  
  # Data type
  input_data1 <- cluster_parameters$input_data1
  input_data2 <- cluster_parameters$input_data2
  
  # Parameters
  algorithm <- as.numeric(cluster_parameters$parameter1_value)
  resolution <- as.numeric(cluster_parameters$parameter2_value)
  random.seed <- as.numeric(cluster_parameters$parameter3_value)
  batch_correction_method <- cluster_parameters$parameter4_value
  
  if (grepl("Kinker", input_dir) & cluster_parameters$parameter6_value == "no") {
    # Create cluster assignments dataframe
    metadata_file <- cluster_parameters$parameter5_value
    cell_metadata <- read.csv(paste0(input_dir, "/intermediate_files/preprocessed/", metadata_file, ".csv"))
    
    clusters <- data.frame(CellID = cell_metadata$CellID,
                           Clusters = NA)
    colnames(clusters) <- c("CellID", paste0("Parameters_",
                                             cluster_parameter_index))
    
    # Time
    time <- NA
  } else {
    
    options(Seurat.object.assay.version = 'v5')
    
    # Input
    if (grepl("bpcells", input_data1)) {
      options(Seurat.object.assay.version = 'v5')
      raw_counts1 <- open_matrix_dir(dir = paste0(input_dir, "/intermediate_files/preprocessed/", input_data1))
      # Kinker et al. downsampling
      if (grepl("Kinker", input_dir)) {
        metadata_file <- cluster_parameters$parameter5_value
        cell_metadata <- read.csv(paste0(input_dir, "/intermediate_files/preprocessed/", metadata_file, ".csv"))
        raw_counts1 <- raw_counts1[,cell_metadata$CellID]
      }
      if (input_data1 == "bpcells_CCA_MNN") {
        norm_counts1 <- raw_counts1
      } else {
        norm_counts1 <- LogNormalize(raw_counts1)
      }
      # If multi-omic
      if (!is.na(input_data2)) {
        raw_counts2 <- open_matrix_dir(dir = paste0(input_dir, "/intermediate_files/preprocessed/", input_data2))
        if (input_data2 == "bpcells_CCA_MNN2") {
          norm_counts2 <- raw_counts2
        } else {
          norm_counts2 <- LogNormalize(raw_counts2)
        }
        norm_counts <- list(norm_counts1, norm_counts2)
        multi_modal <- TRUE
      } else {
        norm_counts <- list(norm_counts1)
        multi_modal <- FALSE
      }
    } else {
      norm_counts1 <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/", input_data1, ".rds"))
      
      # If multi-modal
      if (!is.na(input_data2)) {
        norm_counts2 <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/", input_data2, ".rds"))
        norm_counts <- list(norm_counts1, norm_counts2)
        multi_modal <- TRUE
      } else {
        norm_counts <- list(norm_counts1)
        multi_modal <- FALSE
      }
    }
    
    if (batch_correction_method == "Harmony") {
      cell_metadata <- read.csv(paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata.csv"))
      rownames(cell_metadata) <- cell_metadata$CellID
      shared_cells <- intersect(rownames(cell_metadata), colnames(norm_counts[[1]]))
      cell_metadata <- cell_metadata[shared_cells, ]
      
    }
    
    # Seurat ---------------------------
    
    # Prep
    
    # Start time
    start_time <- Sys.time()
    
    # Dimensionality reduction
    reductions <- list()
    for (i in 1:length(norm_counts)) {
      norm_counts_i <- norm_counts[[i]]
      if (batch_correction_method == "Harmony") {
        norm_counts_i <- norm_counts_i[, shared_cells]
      }
      
      # Scale data
      if ((i == 1 & grepl("SCTransform", input_data1)) | (i == 2 & grepl("SCTransform", input_data2))) {
        var_features <- rownames(norm_counts_i)
        scaled_features <- norm_counts_i[var_features,]
      } else {
        # Find variable features
        var_features <- Seurat::FindVariableFeatures(norm_counts_i, verbose = FALSE) 
        if ("vst.variance.standardized" %in% colnames(var_features)) {
          var_features <- var_features %>%
            dplyr::arrange(-vst.variance.standardized) %>%
            head(2000) %>%
            rownames()
        } else {
          var_features <- var_features %>%
            dplyr::arrange(-variance.standardized) %>%
            head(2000) %>%
            rownames()
        }
        scaled_features <- Seurat::ScaleData(norm_counts_i[var_features,], verbose = FALSE)
      }
      
      # Run PCA
      reduction_coords <- Seurat::RunPCA(object = scaled_features,
                                         assay = "RNA",
                                         features = var_features,
                                         seed.use = random.seed)@cell.embeddings
      # Run Harmony
      if (batch_correction_method == "Harmony") {
        reduction_coords <- do.call(harmony::HarmonyMatrix, c(list("data_mat" = reduction_coords,
                                                                   "meta_data" = cell_metadata,
                                                                   "vars_use" = "Batch",
                                                                   "do_pca" = FALSE)))
      }
      reductions[[i]] <- reduction_coords
    }
    
    # Find neighbors
    if (multi_modal == FALSE) {
      nearest_neighbors <- Seurat::FindNeighbors(object = reductions[[1]])
    } else if (multi_modal == TRUE) {
      # Prep for finding neighbors
      tmp <- matrix(stats::rnorm(nrow(reductions[[1]]) * 3, 10), ncol = nrow(reductions[[1]]), nrow = 3)
      colnames(tmp) <- rownames(reductions[[1]])
      rownames(tmp) <- paste0("t",seq_len(nrow(tmp)))
      tmp_seurat <- Seurat::CreateSeuratObject(tmp, min.cells = 0, min.features = 0, assay = 'tmp')
      dim_list = vector("list", length = length(reductions))
      for (i in 1:length(reductions)) {
        tmp_seurat[[paste0("DR_", i)]] <- Seurat::CreateDimReducObject(embeddings = reductions[[1]],
                                                                       key = paste0("DR_", i, "_"), assay = 'tmp')
        dim_list[[i]] <- 1:ncol(reductions[[1]])
      }
      nearest_neighbors <- Seurat::FindMultiModalNeighbors(object = tmp_seurat,
                                                           reduction.list = list(paste0("DR_", seq(1, length(reductions)))),
                                                           dim.list = dim_list,
                                                           knn.graph.name = "nn",
                                                           snn.graph.name = "snn")@graphs
    }
    
    # Cluster
    if (ncol(nearest_neighbors[["snn"]]) > 30000 & algorithm == 4) {
      method <- "igraph"
    } else {
      method <- "matrix"
    }
    
    cluster_output <- FindClusters(nearest_neighbors[["snn"]],
                                   algorithm = algorithm,
                                   method = method,
                                   resolution = resolution,
                                   random.seed = random.seed)
    
    # End time
    end_time <- Sys.time()
    
    # Output ---------------------------
    
    # Create cluster assignments dataframe
    clusters <- data.frame(CellID = rownames(cluster_output),
                           Clusters = cluster_output[,1])
    colnames(clusters) <- c("CellID", paste0("Parameters_",
                                             cluster_parameter_index))
    
    # Time
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
  }
  # Return cluster assignments dataframe
  return(list("clusters" = clusters,
              "time" = time))
}

