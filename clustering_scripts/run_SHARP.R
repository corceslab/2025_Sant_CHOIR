# ---------------------------------------------------------------------------
# Run SHARP
# ---------------------------------------------------------------------------
library(SHARP)

run_SHARP <- function(input_dir,
                      cluster_parameter_index,
                      cluster_parameters) {

  # Set up ---------------------------

  # Data type
  input_data <- cluster_parameters$input_data1

  # Parameters
  flashmark <- as.logical(cluster_parameters$parameter1_value)
  ensize.K <- as.numeric(cluster_parameters$parameter2_value)
  rN.seed <- as.numeric(cluster_parameters$parameter3_value)
  n.cores <- as.numeric(cluster_parameters$parameter4_value)

  # Input
  counts <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/", input_data, ".rds"))

  # SHARP ---------------------------

  # Prep
  counts <- as.matrix(counts)

  if (grepl("log", input_data)) {
    logflag <- TRUE
  } else {
    logflag <- FALSE
  }

  # Start time
  start_time <- Sys.time()

  # Cluster
  try(object <- SHARP(counts,
                  prep = FALSE,
                  logflag = logflag,
                  flashmark = flashmark,
                  ensize.K = ensize.K,
                  n.cores = n.cores,
                  rN.seed = rN.seed))
  
  if (exists("object")) {
    # End time
    end_time <- Sys.time()
    
    # Output ---------------------------
    
    # Create cluster assignments dataframe
    clusters <- data.frame(CellID = colnames(counts),
                           Clusters = object$pred_clusters)
    
    colnames(clusters) <- c("CellID", paste0("Parameters_", cluster_parameter_index))
    
    # Time
    time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
  } else {
    # Output ---------------------------
    
    # Create cluster assignments dataframe
    clusters <- data.frame(CellID = colnames(counts),
                           Clusters = rep(NA, ncol(counts)))
    
    colnames(clusters) <- c("CellID", paste0("Parameters_", cluster_parameter_index))
    
    # Time
    time <- NA
  }

  # Return cluster assignments dataframe
  return(list("clusters" = clusters,
              "time" = time))
}

