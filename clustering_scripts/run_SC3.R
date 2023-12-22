# ---------------------------------------------------------------------------
# Run SC3
# ---------------------------------------------------------------------------
library(SingleCellExperiment)
library(SC3)

run_SC3 <- function(input_dir,
                    cluster_parameter_index,
                    cluster_parameters) {

  # Set up ---------------------------

  # Data type
  input_data1 <- cluster_parameters$input_data1
  input_data2 <- cluster_parameters$input_data2

  # Parameters
  gene_filter <- as.logical(cluster_parameters$parameter1_value)
  d_region_min <- as.numeric(cluster_parameters$parameter2_value)
  d_region_max <- as.numeric(cluster_parameters$parameter3_value)
  rand_seed <- as.numeric(cluster_parameters$parameter4_value)
  n_cores <- as.numeric(cluster_parameters$parameter5_value)

  # Input
  raw_counts <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/", input_data1, ".rds"))
  norm_counts <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/", input_data2, ".rds"))

  # SC3 ---------------------------

  # Prep

  # Subset genes in raw_counts matrix to those present in norm_counts matrix
  raw_counts <- raw_counts[rownames(norm_counts), ]

  # Create SingleCellExperiment object
  # 'counts' slot is used for the gene dropout rate
  # 'logcounts' slot is used in the main clustering algorithm
  object <- SingleCellExperiment(assays = list(counts = as.matrix(raw_counts),
                                               logcounts = as.matrix(norm_counts)))

  rowData(object)$feature_symbol <- rownames(object)

  # Start time
  start_time <- Sys.time()

  object <- sc3_prepare(object, n_cores = n_cores, rand_seed = rand_seed)

  # Estimate # of clusters k
  try(object_w_k <- sc3_estimate_k(object))
  if (exists("object_w_k")) {
    k_num <- metadata(object_w_k)$sc3$k_estimation
    
    # Cluster
    object <- sc3(object_w_k,
                  ks = k_num,
                  gene_filter = gene_filter,
                  d_region_min = d_region_min,
                  d_region_max = d_region_max,
                  n_cores = n_cores,
                  rand_seed = rand_seed)
    
    # End time
    end_time <- Sys.time()
    
    # Output ---------------------------
    
    # Create cluster assignments dataframe
    clusters <- data.frame(CellID = colnames(object),
                           Clusters = c(colData(object)[, paste0("sc3_", k_num, "_clusters")]))
    
    colnames(clusters) <- c("CellID", paste0("Parameters_", cluster_parameter_index))
    
    # Time
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  } else {
    # Output ---------------------------
    
    # Create cluster assignments dataframe
    clusters <- data.frame(CellID = colnames(norm_counts),
                           Clusters = rep(NA, ncol(norm_counts)))
    
    colnames(clusters) <- c("CellID", paste0("Parameters_", cluster_parameter_index))
    
    # Time
    time <- NA
  }

  # Return cluster assignments dataframe
  return(list("clusters" = clusters,
              "time" = time))
}

