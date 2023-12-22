# ---------------------------------------------------------------------------
# Run CIDR
# ---------------------------------------------------------------------------
library(cidr)
library(dplyr)

run_CIDR <- function(input_dir,
                     cluster_parameter_index,
                     cluster_parameters) {

  # Set up ---------------------------

  # Data type
  input_data1 <- cluster_parameters$input_data1
  input_data2 <- cluster_parameters$input_data2

  # Parameters
  cMethod <- cluster_parameters$parameter1_value
  set_seed <- as.numeric(cluster_parameters$parameter2_value)
  threads <- as.numeric(cluster_parameters$parameter3_value)

  # Input
  raw_counts <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/", input_data1, ".rds"))

  # CIDR ---------------------------

  # Prep
  set.seed(set_seed)

  raw_counts <- as.matrix(raw_counts)

  # Start time
  start_time <- Sys.time()

  object <- scDataConstructor(raw_counts, tagType = "raw")

  object <- object %>%
    determineDropoutCandidates() %>%
    wThreshold()

  # Batch correction
  if (!is.na(input_data2)) {
    batch_corrected_counts <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/", input_data2, ".rds"))
    object@nData <- as.matrix(batch_corrected_counts)
  }

  # Cluster
  object <- object %>%
    scDissim(threads = threads) %>%
    scPCA(plotPC = FALSE) %>%
    nPC()
  try(object_clustered <- scCluster(object, nPC = object@nPC, cMethod = cMethod))
  if (exists('object_clustered')) {
    # End time
    end_time <- Sys.time()

    # Output ---------------------------

    # Create cluster assignments dataframe
    clusters <- data.frame(CellID = colnames(object_clustered@tags),
                           Clusters = object_clustered@clusters)

    colnames(clusters) <- c("CellID", paste0("Parameters_", cluster_parameter_index))

    # Time
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
  } else {
    # Output ---------------------------

    # Create cluster assignments dataframe
    clusters <- data.frame(CellID = colnames(object@tags),
                           Clusters = rep(NA, ncol(object@tags)))

    colnames(clusters) <- c("CellID", paste0("Parameters_", cluster_parameter_index))

    # Time
    time <- NA
  }

  # Return cluster assignments dataframe
  return(list("clusters" = clusters,
              "time" = time))
}

