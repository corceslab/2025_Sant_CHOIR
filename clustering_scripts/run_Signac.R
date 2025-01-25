# ---------------------------------------------------------------------------
# Run Signac clustering
# ---------------------------------------------------------------------------
library(Seurat)
library(Signac)
library(dplyr)

run_Signac <- function(input_dir,
                       cluster_parameter_index,
                       cluster_parameters) {

  # Set up ---------------------------

  # Data type
  input_data1 <- cluster_parameters$input_data1
  input_data2 <- cluster_parameters$input_data2

  # Parameters
  reduction_method <- strsplit(cluster_parameters$parameter1_value, "_")[[1]]
  batch_correction_method <- cluster_parameters$parameter2_value
  algorithm <- as.numeric(cluster_parameters$parameter3_value)
  resolution <- as.numeric(cluster_parameters$parameter4_value)
  random.seed <- as.numeric(cluster_parameters$parameter5_value)  # Add throughout

  # Input
  counts1 <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/", input_data1, ".rds"))
  # Subset to cell_metadata
  cell_metadata <- read.csv(paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata.csv"))
  counts1 <- counts1[, cell_metadata$CellID]
  # If multi-modal
  if (!is.na(input_data2)) {
    counts2 <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/", input_data2, ".rds"))
    # Subset to cell_metadata
    counts2 <- counts2[, cell_metadata$CellID]
    multi_modal <- TRUE
  } else {
    counts <- list(counts1)
    multi_modal <- FALSE
  }

  # Signac ---------------------------

  # Prep

  # If multi-modal, subset to shared cell IDs
  if (multi_modal == TRUE) {
    shared_cell_ids <- intersect(colnames(counts1), colnames(counts2))
    counts1 <- counts1[, shared_cell_ids]
    counts2 <- counts2[, shared_cell_ids]
    counts <- list(counts1, counts2)
  }

  # Start time
  start_time <- Sys.time()

  # Dimensionality reduction
  reductions <- list()
  # For each modality, run the specified reduction method
  for (i in 1:length(counts)) {
    counts_i <- counts[[i]]
    reduction_method_i <- reduction_method[i]

    if (reduction_method_i == "PCA") {
      # Scale data
      if ((i == 1 & grepl("SCTransform", input_data1)) | (i == 2 & grepl("SCTransform", input_data2))) {
        var_features <- rownames(counts_i)
        scaled_features <- counts_i[var_features,]
      } else {
        # Find variable features
        var_features <- Seurat::FindVariableFeatures(counts_i, verbose = FALSE) %>%
          dplyr::arrange(-vst.variance.standardized) %>%
          head(2000) %>%
          rownames()
        scaled_features <- Seurat::ScaleData(counts_i[var_features,], verbose = FALSE)
      }
      # Run PCA
      reduction_coords <- Seurat::RunPCA(object = scaled_features,
                                         assay = "RNA",
                                         features = var_features,
                                         seed.use = random.seed)@cell.embeddings
      reductions[[i]] <- reduction_coords
    } else if (reduction_method_i == "LSI") {
      # Scale data
      scaled_features <- RunTFIDF(counts_i)
      # Find top features
      var_features <- FindTopFeatures(scaled_features, min.cutoff = 'q0') %>%
        arrange(percentile) %>% head(25000)
      # Run LSI
      reduction_coords <- RunSVD(scaled_features, features = var_features)@cell.embeddings
      # Batch correct?
      if (batch_correction_method == "Harmony") {
        reduction_coords <- harmony::HarmonyMatrix(data_mat = reduction_coords,
                                                   meta_data = cell_metadata,
                                                   vars_use = "Batch",
                                                   do_pca = FALSE)
      }
      reductions[[i]] <- reduction_coords
    }
  }

  # Find neighbors
  if (multi_modal == FALSE) {
    nearest_neighbors <- Seurat::FindNeighbors(object = reductions[[1]])
  } else if (multi_modal == TRUE) {
    # Prep for finding neighbors
    tmp <- matrix(stats::rnorm(nrow(reductions[[1]]) * 3, 10), ncol = nrow(reductions[[1]]), nrow = 3)
    colnames(tmp) <- rownames(reductions[[1]])
    rownames(tmp) <- paste0("t",seq_len(nrow(tmp)))
    tmp_seurat <- Seurat::CreateSeuratObject(tmp, min.cells = 0, min.features = 0, assay = 'RNA')
    second_assay <- CreateAssay5Object(counts = tmp)
    tmp_seurat[["ATAC"]] <- second_assay
    dim_list = vector("list", length = length(reductions))
    tmp_seurat[["PCA"]] <- Seurat::CreateDimReducObject(embeddings = reductions[[2]],
                                                        key = "PCA_", assay = "RNA")
    tmp_seurat[["LSI"]] <- Seurat::CreateDimReducObject(embeddings = reductions[[1]],
                                                        key = "LSI_", assay = "ATAC")
    dim_list <- list(1:50,1:50)
    
    nearest_neighbors <- Seurat::FindMultiModalNeighbors(object = tmp_seurat,
                                                         reduction.list = list("PCA", "LSI"),
                                                         dims.list = dim_list,
                                                         knn.graph.name = "nn",
                                                         snn.graph.name = "snn",
                                                         verbose = TRUE)@graphs
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

  # Return cluster assignments dataframe
  return(list("clusters" = clusters,
              "time" = time))
}

