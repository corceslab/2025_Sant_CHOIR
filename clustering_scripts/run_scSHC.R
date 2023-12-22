# ---------------------------------------------------------------------------
# Run sc-SHC
# ---------------------------------------------------------------------------
library(scSHC)

run_scSHC <- function(input_dir,
                      cluster_parameter_index,
                      cluster_parameters) {

  # Set up ---------------------------

  # Data type
  input_data <- cluster_parameters$input_data1

  # Parameters
  alpha_value <- as.numeric(cluster_parameters$parameter1_value)
  set_seed <- as.numeric(cluster_parameters$parameter2_value)
  cores_value <- as.numeric(cluster_parameters$parameter3_value)
  batch_correct <- cluster_parameters$parameter4_value

  # Input
  raw_counts <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/", input_data, ".rds"))

  # sc-SHC ---------------------------

  set.seed(set_seed)

  # Prep
  raw_counts <- as.matrix(raw_counts)
  
  if (batch_correct == "no") {
    batch_value <- NULL
  } else {
    # Get batch labels from cell metadata
    cell_metadata <- read.csv(paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata.csv"))
    cell_metadata <- cell_metadata[match(cell_metadata$CellID, colnames(raw_counts)),]
    batch_value <- cell_metadata$Batch
  }

  # Start time
  start_time <- Sys.time()

  # Cluster
  clusters <- scSHC(raw_counts,
                    alpha = alpha_value,
                    cores = cores_value,
                    batch = batch_value)

  # End time
  end_time <- Sys.time()

  # Output ---------------------------

  # Create cluster assignments dataframe
  clusters <- data.frame(CellID = names(clusters[[1]]),
                         Clusters = clusters[[1]])
  colnames(clusters) <- c("CellID", paste0("Parameters_", cluster_parameter_index))

  # Time
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Return cluster assignments dataframe
  return(list("clusters" = clusters,
              "time" = time))
}