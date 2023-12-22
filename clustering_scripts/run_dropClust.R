# ---------------------------------------------------------------------------
# Run dropClust
# ---------------------------------------------------------------------------
library(SingleCellExperiment)
library(dropClust)

run_dropClust <- function(input_dir,
                          cluster_parameter_index,
                          cluster_parameters) {

  # Set up ---------------------------

  # Data type
  input_data <- cluster_parameters$input_data1

  # Parameters
  method <- cluster_parameters$parameter1_value
  deepSplit <- cluster_parameters$parameter2_value
  batch_correct <- cluster_parameters$parameter3_value
  set_seed <- as.numeric(cluster_parameters$parameter4_value)

  # Input
  raw_counts <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/", input_data, ".rds"))

  # dropClust ---------------------------

  set.seed(set_seed)

  # Depending on batch correction
  if (batch_correct == "no") {
    # Create SingleCellExperiment object
    object <- SingleCellExperiment(assays = list(counts = raw_counts))
    # Start time
    start_time <- Sys.time()
    # Prep
    object <- object %>%
      FilterGenes() %>%
      CountNormalize() %>%
      RankGenes() %>%
      Sampling() %>%
      RankPCAGenes()
  } else if (batch_correct == "yes") {
    # Cell metadata
    metadata <- read.csv(paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata.csv"))
    # Split matrix into batches
    batch_IDs <- unique(metadata$Batch)
    objects <- list()
    for (i in 1:length(batch_IDs)) {
      cell_IDs_i <- metadata$CellID[which(metadata$Batch == batch_IDs[i])]
      objects[[i]] <- SingleCellExperiment(assays = list(counts = raw_counts[, cell_IDs_i]))
    }
    # Start time
    start_time <- Sys.time()
    # Merge
    object <- Merge(objects, use.de.genes = FALSE)
    # Correct
    object <- Correction(object)
  }

  # Cluster
  if (method == "default") {
    object <- Cluster(object, method = method, conf = 0)
  } else if (method == "hclust") {
    object <- Cluster(object, method = method, conf = 0, deepSplit = as.numeric(deepSplit))
  }

  # End time
  end_time <- Sys.time()

  # Output ---------------------------

  # Create cluster assignments dataframe
  clusters <- data.frame(CellID = colnames(object),
                         Clusters = object$ClusterIDs)

  colnames(clusters) <- c("CellID", paste0("Parameters_", cluster_parameter_index))

  # Time
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Return cluster assignments dataframe
  return(list("clusters" = clusters,
              "time" = time))
}

