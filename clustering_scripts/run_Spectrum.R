# ---------------------------------------------------------------------------
# Run Spectrum
# ---------------------------------------------------------------------------
library(Spectrum)

run_Spectrum <- function(input_dir,
                         cluster_parameter_index,
                         cluster_parameters) {

  # Set up ---------------------------

  # Data type
  input_data1 <- cluster_parameters$input_data1
  input_data2 <- cluster_parameters$input_data2

  # Parameters
  method <- as.numeric(cluster_parameters$parameter1_value)
  maxk <- as.numeric(cluster_parameters$parameter2_value)
  set_seed <- as.numeric(cluster_parameters$parameter3_value)

  # Input
  norm_counts <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/", input_data1, ".rds"))
  norm_counts <- as.matrix(norm_counts)

  # Spectrum ---------------------------

  set.seed(set_seed)

  # Start time
  start_time <- Sys.time()

  # Cluster
  try(object <- Spectrum(norm_counts,
                     showres = FALSE,
                     method = method,
                     maxk = maxk))

  if (exists('object')) {
    # End time
    end_time <- Sys.time()

    # Output ---------------------------

    # Create cluster assignments dataframe
    clusters <- data.frame(CellID = colnames(norm_counts),
                           Clusters = object$assignments)

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

