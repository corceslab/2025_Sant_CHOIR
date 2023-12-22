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
  
  options(Seurat.object.assay.version = 'v5')

  # Input
  if (grepl("bpcells", input_data1)) {
    options(Seurat.object.assay.version = 'v5')
    raw_counts <- open_matrix_dir(dir = paste0(input_dir, "/intermediate_files/preprocessed/", input_data1))
    norm_counts <- LogNormalize(raw_counts)
  } else {
    norm_counts <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/", input_data1, ".rds"))
  }
  
  if (batch_correction_method == "Harmony") {
    cell_metadata <- read.csv(paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata.csv"))
    rownames(cell_metadata) <- cell_metadata$CellID
    cell_metadata <- cell_metadata[colnames(norm_counts[[1]]), ]
  }

  # Seurat ---------------------------

  # Prep

  # Start time
  start_time <- Sys.time()

  # Dimensionality reduction
  # Find variable features
  var_features <- Seurat::FindVariableFeatures(norm_counts, verbose = FALSE) 
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
  scaled_features <- Seurat::ScaleData(norm_counts[var_features,], verbose = FALSE)
  
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
  
  # Find neighbors
  nearest_neighbors <- Seurat::FindNeighbors(object = reduction_coords)
  
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

  # Return cluster assignments dataframe
  return(list("clusters" = clusters,
              "time" = time))
}

