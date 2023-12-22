# ---------------------------------------------------------------------------
# Run RaceID3
# ---------------------------------------------------------------------------
library(RaceID)

run_RaceID3 <- function(input_dir,
                        cluster_parameter_index,
                        cluster_parameters) {

  # Set up ---------------------------

  # Data type
  input_data <- cluster_parameters$input_data1

  # Parameters
  FUNcluster <- cluster_parameters$parameter1_value
  clustnr <- as.numeric(cluster_parameters$parameter2_value)
  batch_correct <- cluster_parameters$parameter3_value
  rseed <- as.numeric(cluster_parameters$parameter4_value)
  no_cores <- as.numeric(cluster_parameters$parameter5_value)

  # Input
  raw_counts <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/", input_data, ".rds"))

  # RaceID3 ---------------------------

  # Prep
  object <- SCseq(as.matrix(raw_counts))

  # Batch correction
  if (batch_correct == "yes") {
    # Cell metadata
    metadata <- read.csv(paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata.csv"))
    # Split matrix into batches
    batch_IDs <- unique(metadata$Batch)
    batch_cell_IDs <- list()
    for (i in 1:length(batch_IDs)) {
      cell_IDs_i <- metadata$CellID[which(metadata$Batch == batch_IDs[i])]
      batch_cell_IDs[[i]] <- cell_IDs_i
    }
    # Start time
    start_time <- Sys.time()
    # Filter
    object <- filterdata(object, mintotal=1, LBatch = batch_cell_IDs)
  } else {
    # Start time
    start_time <- Sys.time()
    # Filter
    object <- filterdata(object, mintotal=1)
  }

  # Cluster
  object <- compdist(object,
                     no_cores = no_cores)
  object <- clustexp(object,
                     clustnr = clustnr,
                     rseed = rseed,
                     FUNcluster = FUNcluster)

  # End time
  end_time <- Sys.time()

  # Output ---------------------------

  # Create cluster assignments dataframe
  clusters <- data.frame(CellID = names(object@cluster$kpart),
                         Clusters = object@cluster$kpart)

  colnames(clusters) <- c("CellID", paste0("Parameters_", cluster_parameter_index))

  # Time
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Return cluster assignments dataframe
  return(list("clusters" = clusters,
              "time" = time))
}

