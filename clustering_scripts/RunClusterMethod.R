# ---------------------------------------------------------------------------
# Run clustering methods
# ---------------------------------------------------------------------------
library(dplyr)
library(reticulate)

# Set up ---------------------------

args <- commandArgs(trailingOnly = TRUE)

# Define args
scripts_dir <- args[1]
input_dir <- args[2]
temp_dir <- args[3]
parameter_file <- args[4]
parameter_key <- args[5]
cluster_parameter_index <- args[6]

# Import method list & get parameters for this run
cluster_parameter_list <- read.csv(parameter_file)
cluster_parameter_list$Index <- seq(1:nrow(cluster_parameter_list))
cluster_parameters <- cluster_parameter_list %>% dplyr::filter(Index == as.numeric(cluster_parameter_index))
cluster_method <- cluster_parameters$method

# Source function for running selected method
if (cluster_method %in% c("PanoView", "GiniClust3", "SCCAF", "Cytocipher")) {
  source_python(paste0(scripts_dir, "/run_", cluster_method ,".py"))
} else {
  source(paste0(scripts_dir, "/run_", cluster_method ,".R"))
}

# Run selected method ---------------------------

output <- switch(cluster_method,
                 "ArchR" = run_ArchR(input_dir, cluster_parameter_index, cluster_parameters, temp_dir),
                 "CHOIR" = run_CHOIR(input_dir, cluster_parameter_index, cluster_parameters, temp_dir, parameter_key),
                 "CIDR" = run_CIDR(input_dir, cluster_parameter_index, cluster_parameters),
                 "Cytocipher" = run_Cytocipher(input_dir, cluster_parameter_index,
                                               cluster_parameters$input_data1,
                                               as.numeric(cluster_parameters$parameter1_value),
                                               as.numeric(cluster_parameters$parameter2_value),
                                               as.integer(cluster_parameters$parameter3_value),
                                               as.character(cluster_parameters$parameter4_value),
                                               as.integer(cluster_parameters$parameter5_value),
                                               as.integer(cluster_parameters$parameter6_value),
                                               as.character(cluster_parameters$parameter7_value),
                                               as.character(cluster_parameters$parameter8_value),
                                               as.character(cluster_parameters$parameter9_value)),
                 "dropClust" = run_dropClust(input_dir, cluster_parameter_index, cluster_parameters),
                 "GiniClust3" = run_GiniClust3(input_dir, cluster_parameter_index,
                                               cluster_parameters$input_data1,
                                               cluster_parameters$parameter1_value,
                                               as.numeric(cluster_parameters$parameter2_value),
                                               as.numeric(cluster_parameters$parameter3_value),
                                               as.numeric(cluster_parameters$parameter4_value),
                                               cluster_parameters$parameter5_value,
                                               as.integer(cluster_parameters$parameter6_value),
                                               as.character(cluster_parameters$parameter7_value),
                                               as.character(cluster_parameters$parameter8_value)),
                 "PanoView" = run_PanoView(input_dir, cluster_parameter_index,
                                           cluster_parameters$input_data1,
                                           as.numeric(cluster_parameters$parameter1_value),
                                           as.numeric(cluster_parameters$parameter2_value),
                                           as.integer(cluster_parameters$parameter3_value)),
                 "RaceID3" = run_RaceID3(input_dir, cluster_parameter_index, cluster_parameters),
                 "SAFEclustering" = run_SAFEclustering(input_dir, cluster_parameter_index, cluster_parameters),
                 "SC3" = run_SC3(input_dir, cluster_parameter_index, cluster_parameters),
                 "SCCAF" = run_SCCAF(input_dir, cluster_parameter_index,
                                     cluster_parameters$input_data1,
                                     as.numeric(cluster_parameters$parameter1_value),
                                     as.numeric(cluster_parameters$parameter2_value),
                                     as.integer(cluster_parameters$parameter3_value),
                                     as.character(cluster_parameters$parameter4_value),
                                     as.character(cluster_parameters$parameter5_value),
                                     as.character(cluster_parameters$parameter6_value)),
                 "scCAN" = run_scCAN(input_dir, cluster_parameter_index, cluster_parameters),
                 "scSHC" = run_scSHC(input_dir, cluster_parameter_index, cluster_parameters),
                 "Seurat" = run_Seurat(input_dir, cluster_parameter_index, cluster_parameters),
                 "SHARP" = run_SHARP(input_dir, cluster_parameter_index, cluster_parameters),
                 "Signac" = run_Signac(input_dir, cluster_parameter_index, cluster_parameters),
                 "SIMLR" = run_SIMLR(input_dir, cluster_parameter_index, cluster_parameters),
                 "Spectrum" = run_Spectrum(input_dir, cluster_parameter_index, cluster_parameters))

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
          paste0(input_dir, "/intermediate_files/clusters/", parameter_key, "/clusters_parameter_set_", cluster_parameter_index, "_", cluster_parameters$snakemake_category[1], ".csv"),
          row.names = FALSE)
write.csv(time_output,
          paste0(input_dir, "/intermediate_files/clusters/", parameter_key, "/time_parameter_set_", cluster_parameter_index, "_", cluster_parameters$snakemake_category[1], ".csv"),
          row.names = FALSE)
