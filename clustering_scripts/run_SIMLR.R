# ---------------------------------------------------------------------------
# Run SIMLR
# ---------------------------------------------------------------------------
library(SIMLR)
library(igraph)

run_SIMLR <- function(input_dir,
                      cluster_parameter_index,
                      cluster_parameters) {

  # Set up ---------------------------

  # Data type
  input_data <- cluster_parameters$input_data1

  # Parameters
  NUMC_max <- as.numeric(cluster_parameters$parameter1_value)
  set_seed <- as.numeric(cluster_parameters$parameter2_value)

  # Input
  norm_counts <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/", input_data, ".rds"))

  # SIMLR ---------------------------

  set.seed(set_seed)

  # Prep
  norm_counts <- as.matrix(norm_counts)

  # Start time
  start_time <- Sys.time()

  # Estimate number of clusters
  NUMC <- 2:NUMC_max
  estimate_n_clusters <- SIMLR_Estimate_Number_of_Clusters(norm_counts,
                                                           NUMC = NUMC,
                                                           cores.ratio = 0)
  K1_n_clusters <- NUMC[which.min(estimate_n_clusters$K1)]
  K2_n_clusters <- NUMC[which.min(estimate_n_clusters$K2)]
  n_clusters <- round(mean(K1_n_clusters, K2_n_clusters))

  # Cluster
  object <- SIMLR(X = norm_counts,
                  c = n_clusters,
                  cores.ratio = 0)

  # End time
  end_time <- Sys.time()

  # Output ---------------------------

  # Create cluster assignments dataframe
  clusters <- data.frame(CellID = colnames(norm_counts),
                         Clusters = object$y$cluster)

  colnames(clusters) <- c("CellID", paste0("Parameters_", cluster_parameter_index))

  # Time
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Return cluster assignments dataframe
  return(list("clusters" = clusters,
              "time" = time))
}

