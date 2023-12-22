# ---------------------------------------------------------------------------
# Run scCAN
# ---------------------------------------------------------------------------
library(scCAN)

run_scCAN <- function(input_dir,
                      cluster_parameter_index,
                      cluster_parameters) {

  # Set up ---------------------------

  # Data type
  input_data <- cluster_parameters$input_data1

  # Parameters
  k_max <- as.numeric(cluster_parameters$parameter1_value)
  samp.size <- as.numeric(cluster_parameters$parameter2_value)
  r.seed <- as.numeric(cluster_parameters$parameter3_value)
  ncores <- as.numeric(cluster_parameters$parameter4_value)

  # Input
  norm_counts <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/", input_data, ".rds"))

  # scCAN ---------------------------

  # Prep
  norm_counts_t <- t(as.matrix(norm_counts))

  # Start time
  start_time <- Sys.time()

  try(object <- scCAN(norm_counts_t,
                  samp.size = samp.size,
                  r.seed = r.seed,
                  ncores = ncores,
                  k = 2:k_max))

  if (exists('object')) {
    # End time
    end_time <- Sys.time()

    # Output ---------------------------

    # Create cluster assignments dataframe
    clusters <- data.frame(CellID = colnames(norm_counts),
                           Clusters = object$cluster)

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

