# ---------------------------------------------------------------------------
# This script assesses benchmarking metrics for the subclustering analysis of
# Simulated Datasets 46-50
# ---------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(Seurat)

input_dir <- "/gladstone/corces/lab/users/cpetersen/cluster_benchmarking/Subsetted_Splatter"

# Metrics
datasets <- c("sim_25000_1", "sim_25000_2", "sim_25000_3", "sim_25000_4", "sim_25000_5")
parameter_lists <- c("CHOIR_parameter_list_RNA_simulated",
                     "cluster_parameter_list_RNA_simulated")

# Create dataframe
cluster_stats <- data.frame(dataset = NULL,
                            parameter_list = NULL,
                            method = NULL,
                            index = NULL,
                            cluster = NULL,
                            run = NULL,
                            n_clusters = NULL,
                            ground_truth_correction = NULL)

# For each dataset
for (dataset_i in 1:5) {
  cell_metadata_i <- read.csv(paste0(input_dir, "/", datasets[dataset_i], 
                                     "/intermediate_files/preprocessed/cell_metadata.csv"))
  # For each parameter list
  for (parameter_list_j in 1:6) {
    print(parameter_list_j)
    # Import all cluster results and find the number of clusters
    # and whether the misclassified cells are limited to those cells originally misclassified
    current_directory <- paste0(input_dir, "/", datasets[dataset_i], 
                                "/intermediate_files/clusters/",
                                parameter_lists[parameter_list_j])
    current_files <- list.files(current_directory)
    current_files <- current_files[grepl(paste(paste0("clusters_parameter_set_", 
                                                      index_lists[[parameter_list_j]], "_"),
                                               collapse = "|"),
                                         current_files)]
    # For each file
    for (f in current_files) {
      if (grepl("not_run", f)) {
        cluster_stats_f <- data.frame(dataset = datasets[dataset_i],
                                      parameter_list = parameter_lists[parameter_list_j],
                                      index = sub("_", "", sub("_\\d", "", substr(sub("clusters_parameter_set_", "", f), 1, 3))),
                                      cluster = substr(sub("clusters_parameter_set_\\d*_", "", f), 1, 1),
                                      run = FALSE,
                                      n_clusters = NA,
                                      ground_truth_correction = NA)
      } else {
        clusters_f <- read.csv(paste0(current_directory, "/", f))
        if (any(is.na(clusters_f[,2]))) {
          n_clust <- NA
        } else {
          n_clust <- dplyr::n_distinct(clusters_f[,2])
        }
        percent_max_cluster <- NA
        if (!is.na(n_clust)) {
          if (n_clust != 1) {
            # If we remove the primary cluster, are the remaining cells those that were already misclassified?
            max_cluster <- names(table(clusters_f[,2])[which(table(clusters_f[,2]) == max(table(clusters_f[,2])))])[1]
            max_cluster_cells <- clusters_f[clusters_f[,2] == max_cluster,]$CellID
            minor_cluster_cells <- clusters_f[clusters_f[,2] != max_cluster,]$CellID
            # Major ground truth group for max cluster
            ground_truth_table <- table(cell_metadata_i[cell_metadata_i$CellID %in% max_cluster_cells, "Ground_truth"])
            max_ground_truth <- names(ground_truth_table[which(ground_truth_table == max(ground_truth_table))])
            # Do any of the minor cluster cells belong to the max cluster?
            percent_max_cluster <- sum(cell_metadata_i[cell_metadata_i$CellID %in% minor_cluster_cells, "Ground_truth"] == max_ground_truth)/
              length(minor_cluster_cells)
          }
        }
        cluster_stats_f <- data.frame(dataset = datasets[dataset_i],
                                      parameter_list = parameter_lists[parameter_list_j],
                                      index = sub("_", "", sub("_\\d", "", substr(sub("clusters_parameter_set_", "", f), 1, 3))),
                                      cluster = substr(sub("clusters_parameter_set_\\d*_", "", f), 1, 1),
                                      run = TRUE,
                                      n_clusters = n_clust,
                                      ground_truth_correction = percent_max_cluster)
      }
      cluster_stats <- rbind(cluster_stats, cluster_stats_f)
    }
  }
}

# Plot
num_limits = c(0,22)
num_breaks = c(seq(5,20), 22)
num_labels = c(seq(5,20), ">20")
cap_max = 20
cap_loc = 22
cluster_stats %>% group_by(method, index, dataset, default) %>%
  dplyr::filter(dataset == "sim_25000_5") %>%
  summarise(sum_clusters = sum(n_clusters)) %>%
  mutate(n_clusters_capped = ifelse(sum_clusters > cap_max, cap_loc, sum_clusters)) %>%
  ggplot(aes(x = method, y = n_clusters_capped, alpha = default)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept = 5, color = "red") +
  geom_beeswarm(shape = 21, fill = "white", cex = 0.25) + 
  geom_point(data = . %>% filter(default == TRUE), shape = 23, alpha = 1, fill = "gold") +
  xlab("Method") +
  ylab("Number of clusters") +
  scale_x_discrete(limits = c("CHOIR", "CIDR", "Cytocipher","GiniClust3",
                              "SCCAF", "scSHC", "Seurat"),
                   breaks = c("CHOIR", "CIDR", "Cytocipher","GiniClust3",
                              "SCCAF", "scSHC", "Seurat")) +
  NoLegend() +
  scale_alpha_manual(values = c(1,0)) +
  scale_y_continuous(breaks = num_breaks, labels = num_labels, limits = num_limits)