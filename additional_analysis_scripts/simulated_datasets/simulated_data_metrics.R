# ---------------------------------------------------------------------------
# This script assesses benchmarking metrics for the clustering of each 
# simulated dataset
# ---------------------------------------------------------------------------
library(pdfCluster)
library(dplyr)
library(Seurat)
library(pbmcapply)

for (groups in c(1,5,10,20)) {
  for (size in 1:5) {
    if (groups == 1) {
      size_val <- c(500,1000,1500,2000,2500)[size]
    } else {
      size_val <- c(5000,10000,15000,20000,25000)[size]
    }
    for (i in 1:5) {
      print(paste0(groups, "_", size_val, "_", i))
      parameter_dir <- "/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/parameter_lists"
      input_dir <- paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", 
                          groups, "/sim_", size_val, "_", i)
      type <- "RNA_simulated"
      ground_truth <- TRUE
      
      cluster_parameter_list_compiled <- list()
      for (file in c(1:2)) {
        if (file == 1) {
          cluster_parameter_list <- read.csv(paste0(parameter_dir, "/CHOIR/CHOIR_parameter_list_RNA_simulated.csv"))
          cluster_df <- read.csv(paste0(input_dir, "/intermediate_files/clusters/compiled_clusters_CHOIR_parameter_list_RNA_simulated.csv"))
          time_df <- read.csv(paste0(input_dir, "/intermediate_files/clusters/compiled_time_CHOIR_parameter_list_RNA_simulated.csv"))
          cluster_parameter_list$list <- "CHOIR"
        } else if (file == 2) {
          cluster_parameter_list <- read.csv(paste0(parameter_dir, "/other_methods/cluster_parameter_list_RNA_simulated.csv"))
          cluster_df <- read.csv(paste0(input_dir, "/intermediate_files/clusters/compiled_clusters_cluster_parameter_list_RNA_simulated.csv"))
          time_df <- read.csv(paste0(input_dir, "/intermediate_files/clusters/compiled_time_cluster_parameter_list_RNA_simulated.csv"))
          cluster_parameter_list$list <- "other_methods"
        }
        cluster_parameter_list$index <- seq(1:nrow(cluster_parameter_list))
        
        metadata <- read.csv(paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata.csv"))
        cluster_parameter_list <- cluster_parameter_list %>%
          mutate(time = NA,
                 total_cpu_time = NA,
                 maxvmem = NA,
                 mem = NA,
                 vmem = NA,
                 io = NA,
                 status = NA,
                 n_clusters = NA,
                 correct_n_clusters = NA,
                 ARI = NA,
                 ARI_0.9 = NA,
                 time = NA,
                 min_DEGs = NA)
        
        if (ground_truth == TRUE) {
          n_ground_truth <- dplyr::n_distinct(metadata$Ground_truth)
          ground_truth_clusters <- unique(metadata$Ground_truth)
        }
        
        for (i in 1:(ncol(cluster_df)-1)) {
          index_i <- as.numeric(sub("Parameters_", "", colnames(cluster_df)[i + 1]))
          
          time_i <- time_df[i,2]
          cluster_parameter_list$time[index_i] <- time_i
          
          total_cpu_time_i <- time_df[i,3]
          cluster_parameter_list$total_cpu_time[index_i] <- total_cpu_time_i
          
          maxvmem_i <- time_df[i,4]
          cluster_parameter_list$maxvmem[index_i] <- maxvmem_i
          
          mem_i <- time_df[i,5]
          cluster_parameter_list$mem[index_i] <- mem_i
          
          vmem_i <- time_df[i,6]
          cluster_parameter_list$vmem[index_i] <- vmem_i
          
          io_i <- time_df[i,7]
          cluster_parameter_list$io[index_i] <- io_i
          
          status_i <- ifelse(is.na(time), "out of time",
                             ifelse(!is.na(cluster_df[1, i + 1]), "complete", "error"))
          cluster_parameter_list$status[index_i] <- status_i
          
          if (status_i == "complete") {
            clusters_i <- cluster_df[,c(1, i + 1)]
            clusters_i <- clusters_i[match(metadata$CellID, clusters_i$CellID), ]
            M_clusters <- unique(clusters_i[,2])
            M <- dplyr::n_distinct(clusters_i[,2])
            
            cluster_parameter_list$n_clusters[index_i] <- M
            
            cluster_parameter_list$time[index_i] <- dplyr::filter(time_df, parameter_index == index_i)$time_secs
            
            if (ground_truth == TRUE) {
              cluster_parameter_list$correct_n_clusters[index_i] <- M == n_ground_truth
              cluster_parameter_list$ARI[index_i] <- pdfCluster::adj.rand.index(as.numeric(as.factor(clusters_i[,2])),
                                                                                as.numeric(as.factor(metadata$Ground_truth)))
              cluster_parameter_list$ARI_0.9[index_i] <- cluster_parameter_list$ARI[index_i] >= 0.9
            }
          }
        }
        cluster_parameter_list_compiled[[file]] <- cluster_parameter_list
      }
      cluster_parameter_list <- do.call(rbind, cluster_parameter_list_compiled)
      write.csv(cluster_parameter_list,
                paste0(input_dir, "/intermediate_files/clusters/compiled_clusters_metrics.csv"),
                row.names = FALSE)
    }
  }
}