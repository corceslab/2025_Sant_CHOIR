# ---------------------------------------------------------------------------
# Compile cluster results
# ---------------------------------------------------------------------------
library(dplyr)

# Set up ---------------------------

args <- commandArgs(trailingOnly = TRUE)

# Define args
input_dir <- args[1]
parameter_key <- args[2]

# Cluster files
files <- list.files(path = paste0(input_dir, "/intermediate_files/clusters/", parameter_key), pattern = "^clusters_")
if (length(files) > 0) {
  # Get cell IDs
  cluster_df <- data.frame(CellID = read.csv(paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata.csv"))$CellID)
  # Version of cell IDs with no special characters
  cluster_df$CellID_alt <- make.names(cluster_df$CellID)

  for (i in 1:length(files)) {
    file_i <- read.csv(paste0(input_dir, "/intermediate_files/clusters/", parameter_key, "/", files[i]))
    colnames(file_i) <- c("CellID", paste0("Parameters_", sub("_[AB]\\.csv", "", sub("clusters_parameter_set_", "", files[i]))))
    if (nrow(file_i) == length(intersect(file_i$CellID, cluster_df$CellID))) {
      cluster_df <- merge(cluster_df, file_i, by = "CellID", all = TRUE)
    } else if (nrow(file_i) == length(intersect(file_i$CellID, cluster_df$CellID_alt))) {
      cluster_df <- merge(cluster_df, file_i, by.x = "CellID_alt", by.y = "CellID", all = TRUE)
    } else {
      stop(paste0("Error with cell IDs for file ", files[i]))
    }
  }
  cluster_df <- cluster_df %>% dplyr::select(-CellID_alt)
  write.csv(cluster_df, paste0(input_dir, "/intermediate_files/clusters/compiled_clusters_", parameter_key, ".csv"), row.names = FALSE)
}

# Time files
time_files <- list.files(path = paste0(input_dir, "/intermediate_files/clusters/", parameter_key), pattern = "^time_")
if (length(files) > 0) {
  time_df <- data.frame(parameter_index = NULL,
                        time_secs = NULL,
                        total_cpu_time = NULL,
                        maxvmem = NULL)

  for (i in 1:length(files)) {
    time_i <- read.csv(paste0(input_dir, "/intermediate_files/clusters/", parameter_key, "/", time_files[i]))

    parameter_index <- time_i$parameter_index
    usage_file_i <- list.files(path = paste0(input_dir, "/intermediate_files/clusters/", parameter_key), pattern = paste0("^usage_parameter_set_", parameter_index, "_"))
    usage_i <- read.table(paste0(input_dir, "/intermediate_files/clusters/", parameter_key, "/", usage_file_i)) %>%
      transmute(total_cpu_time = sub(",", "", sub("cpu=", "", V3)),
                maxvmem = as.numeric(sub("M", "", sub("G", "", sub("maxvmem=", "", V10)))),
                maxvmem_unit = ifelse(grepl("M", V10), "M", "G"),
                mem = as.numeric(sub("mem=", "", V4)),
                mem_unit = V5,
                vmem = as.numeric(sub(",", "", sub("M", "", sub("G", "", sub("vmem=", "", V9))))),
                vmem_unit = ifelse(grepl("M", V9), "M", "G"),
                io = as.numeric(sub("io=", "", V7)),
                io_unit = sub(",", "", V8)) %>%
      mutate(mem = ifelse(mem_unit == "MB", mem/1000, mem),
             io = ifelse(io_unit == "MB", io/1000, io),
             maxvmem = ifelse(maxvmem_unit == "M", maxvmem/1000, maxvmem),
             vmem = ifelse(vmem_unit == "M", vmem/1000, vmem)) %>%
      dplyr::select(-c(mem_unit, io_unit, maxvmem_unit, vmem_unit))
    time_i <- cbind(time_i, usage_i)
    time_df <- rbind(time_df, time_i)
  }
  write.csv(time_df, paste0(input_dir, "/intermediate_files/clusters/compiled_time_", parameter_key, ".csv"), row.names = FALSE)
}
