# ---------------------------------------------------------------------------
# This script runs the subclustering analysis on Simulated Datasets 46-50
# ---------------------------------------------------------------------------
library(dplyr)
library(reticulate)

# Set up ---------------------------

args <- commandArgs(trailingOnly = TRUE)

args <- strsplit(args[1], "\\.")[[1]]
# Define args
input_dataset <- args[1]
method <- args[2]
parameter_key <- args[3]
parameter_key_name <- args[4]
cluster_parameter_index <- args[5]
cluster_number <- args[6]

# Set variables
input_dir <- paste0("/gladstone/mucke/sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_5/",
                    input_dataset)
output_dir <- paste0("/gladstone/mucke/sequencing/Cathrine/cluster_benchmarking/datasets/Subsetted_Splatter/",
                     input_dataset)
scripts_dir <- "/wynton/group/muckelab/user/cpetersen/cluster_benchmarking/snakemake_scripts/cluster"
temp_dir <- "/wynton/group/muckelab/user/cpetersen/cluster_benchmarking/temp"
if (method == "CHOIR") {
  parameter_file <- paste0("/gladstone/mucke/sequencing/Cathrine/cluster_benchmarking/parameter_lists/CHOIR/",
                           parameter_key, ".csv")
} else {
  parameter_file <- paste0("/gladstone/mucke/sequencing/Cathrine/cluster_benchmarking/parameter_lists/other_methods/",
                           parameter_key, ".csv")
}

# Import method list & get parameters for this run
cluster_parameter_list <- read.csv(parameter_file)
cluster_parameter_list$Index <- seq(1:nrow(cluster_parameter_list))
cluster_parameters <- cluster_parameter_list %>% dplyr::filter(Index == as.numeric(cluster_parameter_index))
cluster_method <- cluster_parameters$method

# Import
original_clusters <- read.csv(paste0(input_dir, "/intermediate_files/clusters/compiled_clusters_", parameter_key,".csv"))
cluster_metrics <- read.csv(paste0(input_dir, "/intermediate_files/clusters/compiled_clusters_metrics.csv"))
current_index_metrics <- cluster_metrics %>% 
  dplyr::filter(list == parameter_key_name,
                index == cluster_parameter_index)

if (current_index_metrics$correct_n_clusters == TRUE & 
    current_index_metrics$ARI_0.9 == TRUE &
    current_index_metrics$status == "complete") {
  proceed <- TRUE
  print("Proceed")
  
  # Generate subsetted matrix
  cell_IDs <- original_clusters[original_clusters[, paste0("Parameters_", cluster_parameter_index)] == cluster_number,]$CellID
  if (method %in% c("RaceID3", "scSHC_new", "SC3", "CIDR", "dropClust", "SHARP")) {
    count_matrix <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/RNA_matrix_raw.rds"))
    count_matrix <- count_matrix[,cell_IDs]
    saveRDS(count_matrix, paste0(output_dir, "/intermediate_files/preprocessed/RNA_matrix_raw_", 
                                method, "_", cluster_parameter_index, "_", cluster_number, ".rds"))
    input_data <- paste0("RNA_matrix_raw_", method, "_", cluster_parameter_index, "_", cluster_number)
    cluster_parameters$input_data1 <- input_data
  } 
  if (!method %in% c("RaceID3", "scSHC_new", "CIDR", "dropClust", "SHARP")) {
    norm_matrix <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/RNA_matrix_norm_log.rds"))
    norm_matrix <- norm_matrix[,cell_IDs]
    if (method %in% c("Cytocipher", "PanoView", "GiniClust3", "SCCAF")) {
      write.csv(as.matrix(norm_matrix), paste0(output_dir, "/intermediate_files/preprocessed/RNA_matrix_norm_log_", 
                                    method, "_", cluster_parameter_index, "_", cluster_number, ".csv"))
    } else {
      saveRDS(norm_matrix, paste0(output_dir, "/intermediate_files/preprocessed/RNA_matrix_norm_log_", 
                                  method, "_", cluster_parameter_index, "_", cluster_number, ".rds"))
    }
    input_data <- paste0("RNA_matrix_norm_log_", method, "_", cluster_parameter_index, "_", cluster_number)
  }
  cluster_parameters$input_data1 <- input_data
  if (method == "SC3") {
    cluster_parameters$input_data1 <- paste0("RNA_matrix_raw_", method, "_", cluster_parameter_index, "_", cluster_number)
    cluster_parameters$input_data1 <- paste0("RNA_matrix_norm_log_", method, "_", cluster_parameter_index, "_", cluster_number)
  }

} else {
  print("Skip this index")
  proceed <- FALSE
}

if (proceed == TRUE) {
  # Source function for running selected method
  if (cluster_method %in% c("PanoView", "GiniClust3", "SCCAF", "Cytocipher")) {
    source_python(paste0(scripts_dir, "/run_", cluster_method ,".py"))
  } else {
    source(paste0(scripts_dir, "/run_", cluster_method ,".R"))
  }
  
  print("Sourced python.")
  
  # Run selected method ---------------------------
  
  try(output <- switch(cluster_method,
                   "CHOIR" = run_CHOIR(output_dir, cluster_parameter_index, cluster_parameters, temp_dir, parameter_key),
                   "CIDR" = run_CIDR(output_dir, cluster_parameter_index, cluster_parameters),
                   "Cytocipher" = run_Cytocipher(as.character(output_dir), as.character(cluster_parameter_index),
                                                 as.character(input_data),
                                                 as.numeric(cluster_parameters$parameter1_value),
                                                 as.numeric(cluster_parameters$parameter2_value),
                                                 as.integer(cluster_parameters$parameter3_value),
                                                 as.character(cluster_parameters$parameter4_value),
                                                 as.integer(cluster_parameters$parameter5_value),
                                                 as.integer(cluster_parameters$parameter6_value),
                                                 as.character(cluster_parameters$parameter7_value)),
                   "dropClust" = run_dropClust(output_dir, cluster_parameter_index, cluster_parameters),
                   "GiniClust3" = run_GiniClust3(as.character(output_dir), as.character(cluster_parameter_index),
                                                 as.character(input_data),
                                                 cluster_parameters$parameter1_value,
                                                 as.numeric(cluster_parameters$parameter2_value),
                                                 as.numeric(cluster_parameters$parameter3_value),
                                                 as.numeric(cluster_parameters$parameter4_value),
                                                 cluster_parameters$parameter5_value,
                                                 as.integer(cluster_parameters$parameter6_value)),
                   "PanoView" = run_PanoView(output_dir, cluster_parameter_index,
                                             input_data,
                                             as.numeric(cluster_parameters$parameter1_value),
                                             as.numeric(cluster_parameters$parameter2_value),
                                             as.integer(cluster_parameters$parameter3_value)),
                   "RaceID3" = run_RaceID3(output_dir, cluster_parameter_index, cluster_parameters),
                   "SC3" = run_SC3(output_dir, cluster_parameter_index, cluster_parameters),
                   "SCCAF" = run_SCCAF(as.character(output_dir), as.character(cluster_parameter_index),
                                       as.character(input_data),
                                       as.numeric(cluster_parameters$parameter1_value),
                                       as.numeric(cluster_parameters$parameter2_value),
                                       as.integer(cluster_parameters$parameter3_value),
                                       as.character(cluster_parameters$parameter4_value)),
                   "scCAN" = run_scCAN(output_dir, cluster_parameter_index, cluster_parameters),
                   "scSHC" = run_scSHC(output_dir, cluster_parameter_index, cluster_parameters),
                   "Seurat" = run_Seurat(output_dir, cluster_parameter_index, cluster_parameters),
                   "SHARP" = run_SHARP(output_dir, cluster_parameter_index, cluster_parameters),
                   "SIMLR" = run_SIMLR(output_dir, cluster_parameter_index, cluster_parameters),
                   "Spectrum" = run_Spectrum(output_dir, cluster_parameter_index, cluster_parameters)))
  
  # Output ---------------------------
  
  if (cluster_method %in% c("PanoView", "GiniClust3", "SCCAF", "Cytocipher")) {
    clusters <- data.frame("CellID" = output[[1]],
                           "Cluster" = output[[2]])
    colnames(clusters) <- c("CellID", paste0("Parameters_", cluster_parameter_index))
    
    # Timekeeping records
    time <- output[[3]]
    time_output <- data.frame(parameter_index = cluster_parameter_index,
                              time_secs = as.numeric(time))
  } else {
    clusters <- output[["clusters"]]
    # Timekeeping records
    time_output <- data.frame(parameter_index = cluster_parameter_index,
                              time_secs = output[["time"]])
  }
  
  # Write CSVs
  write.csv(clusters,
            paste0(output_dir, "/intermediate_files/clusters/", parameter_key, "/clusters_parameter_set_", cluster_parameter_index, "_", cluster_number, ".csv"),
            row.names = FALSE)
  write.csv(time_output,
            paste0(output_dir, "/intermediate_files/clusters/", parameter_key, "/time_parameter_set_", cluster_parameter_index, "_", cluster_number, ".csv"),
            row.names = FALSE)
} else {
  write.csv("Not run", paste0(output_dir, "/intermediate_files/clusters/", parameter_key, "/clusters_parameter_set_", cluster_parameter_index, "_", cluster_number, "_not_run.csv"))
}
