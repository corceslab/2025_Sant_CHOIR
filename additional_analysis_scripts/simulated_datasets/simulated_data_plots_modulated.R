# ---------------------------------------------------------------------------
# This script generated plots for the simulated datasets that have features
# that have been incrementally modulated
# ---------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggbeeswarm)

simulation_set <- "Splatter_5_DE"

iteration <- 1
alt_val_set <- c(0.4,0.2,0.1,0.05,0.025,0.0125)

input_dir <- paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/", simulation_set)

# Import original metrics
original_metrics_pre <- read.csv(paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_5/sim_25000_", iteration, "/intermediate_files/clusters/compiled_clusters_metrics.csv"))

original_metrics_pre %>% dplyr::filter(status != "complete") %>% 
  dplyr::mutate(mem_na = is.na(mem)) %>%
  dplyr::group_by(method, status, mem_na) %>% summarise(n = n())

# Filter by performance
original_metrics <- original_metrics_pre %>% 
  dplyr::filter(status == "complete") %>%
  dplyr::filter(correct_n_clusters == TRUE, ARI_0.9 == TRUE)

original_metrics$alt_value <- "orig"
original_metrics$method_index <- paste0(original_metrics$method, original_metrics$index)

compiled_metrics <- original_metrics

# Import metrics across values
for (alt_v in 1:length(alt_val_set)) {
  current_metrics <- read.csv(paste0(input_dir, "/sim_25000_", iteration, "_", alt_val_set[alt_v], "/intermediate_files/clusters/compiled_clusters_metrics.csv"))
  current_metrics <- current_metrics %>% 
    dplyr::filter(status == "complete") %>%
    dplyr::filter(correct_n_clusters == TRUE, ARI_0.9 == TRUE)
  if (nrow(current_metrics) > 0) {
    current_metrics$alt_value <- alt_val_set[alt_v]
    current_metrics$method_index <- paste0(current_metrics$method, current_metrics$index)
    if (alt_v == 1) {
      current_metrics <- current_metrics %>% dplyr::filter(method_index %in% original_metrics$method_index)
    } else {
      current_metrics <- current_metrics %>% dplyr::filter(method_index %in% dplyr::filter(compiled_metrics, alt_value == alt_val_set[alt_v - 1])$method_index)
    }
    compiled_metrics <- rbind(compiled_metrics, current_metrics)
  }
}

results <- compiled_metrics %>% 
  dplyr::filter(alt_value != "orig") %>%
  group_by(method, method_index, default) %>%
  summarise(min_de = min(alt_value)) %>%
  arrange(default)
additional_results <- original_metrics %>% dplyr::filter(!(method_index %in% results$method_index))
if (dataset_group == "Splatter_5_DE" & nrow(additional_results) > 0) {
  additional_results$min_de <- 0.4
} else if (dataset_group == "Splatter_5_library_size" & nrow(additional_results) > 0) {
  additional_results$min_de <- 11
} else if (dataset_group == "Splatter_5_downsample" & nrow(additional_results) > 0) {
  additional_results$min_de <- 500
} 

original_metrics_pre$method_index <- paste0(original_metrics_pre$method, original_metrics_pre$index)
original_results <- original_metrics_pre %>% 
  dplyr::filter(!(method_index %in% results$method_index)) %>% 
  dplyr::filter(!(method_index %in% additional_results$method_index))
original_results$min_de <- 1

if (dataset_group == "Splatter_5_DE" & nrow(additional_results) > 0) {
  results$min_de <- as.numeric(results$min_de)/2
} else if (dataset_group == "Splatter_5_library_size" & nrow(additional_results) > 0) {
  results$min_de <- as.numeric(results$min_de) - 0.2
} else if (dataset_group == "Splatter_5_downsample" & nrow(additional_results) > 0) {
  results$min_de <- as.numeric(results$min_de) - 1000
} 

if (nrow(additional_results) > 0) {
  results <- rbind(results, dplyr::select(additional_results, method, method_index, default, min_de))
}
if (nrow(original_results) > 0) {
  results <- rbind(results, dplyr::select(original_results, method, method_index, default, min_de))
}

if (dataset_group == "Splatter_5_DE" & nrow(additional_results) > 0) {
  results <- results %>% mutate(points = ifelse(min_de == 1, 0,
                                                ifelse(min_de == 0.4, 1,
                                                       ifelse(min_de == 0.2, 2,
                                                              ifelse(min_de == 0.1, 3,
                                                                     ifelse(min_de == 0.05, 4,
                                                                            ifelse(min_de == 0.025, 5,
                                                                                   ifelse(min_de == 0.0125, 6,
                                                                                          ifelse(min_de == 0.00625, 7, 8)))))))))
} else if (dataset_group == "Splatter_5_library_size" & nrow(additional_results) > 0) {
  results <- results %>% mutate(points = ifelse(min_de == 20, 0,
                                                ifelse(min_de == 11, 1,
                                                       ifelse(min_de == 10.8, 2,
                                                              ifelse(min_de == 10.6, 3,
                                                                     ifelse(min_de == 10.4, 4,
                                                                            ifelse(min_de == 10.2, 5,
                                                                                   ifelse(min_de == 10, 6,
                                                                                          ifelse(min_de == 9.8, 7, 
                                                                                                 ifelse(min_de == 9.6, 8, 
                                                                                                        ifelse(min_de == 9.4, 9,
                                                                                                               ifelse(min_de == 9.2, 10, 11))))))))))))
} else if (dataset_group == "Splatter_5_downsample" & nrow(additional_results) > 0) {
  results <- results %>% mutate(points = ifelse(min_de == 1000, 0,
                                                ifelse(min_de == 500, 1,
                                                       ifelse(min_de == 400, 2,
                                                              ifelse(min_de == 300, 3,
                                                                     ifelse(min_de == 200, 4,
                                                                            ifelse(min_de == 100, 5, 6)))))))
} 

# Repeat everything above for each iteration
results_1 <- results %>% dplyr::mutate(iteration = 1)

# Then merge
all_results <- rbind(results_1, results_2)
all_results <- rbind(all_results, results_3)
all_results <- rbind(all_results, results_4)
all_results <- rbind(all_results, results_5)

# Identify best performer for each method
sum_results <- all_results %>% group_by(method, method_index, default) %>%
  summarise(sum_points = sum(points))
sum_results <- sum_results %>% group_by(method) %>% slice_max(sum_points, n = 1, with_ties = FALSE)

# Plot
if (dataset_group == "Splatter_5_DE" & nrow(additional_results) > 0) {
  break_set <- factor(c(0.00625, rev(alt_val_set), 1))
} else if (dataset_group == "Splatter_5_library_size" & nrow(additional_results) > 0) {
  break_set <- factor(c(9.0, rev(alt_val_set), 20))
} else if (dataset_group == "Splatter_5_downsample" & nrow(additional_results) > 0) {
  break_set <- factor(c(1000, alt_val_set, 0))
} 

results_1 %>% dplyr::filter(default == TRUE | method_index %in% sum_results$method_index) %>%
  ggplot(aes(x = method, y = as.factor(as.numeric(min_de)), alpha = !default)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_beeswarm(shape = 21, fill = "white", cex = 0.25) + 
  geom_point(data = . %>% filter(default == TRUE), shape = 23, alpha = 1, fill = "gold") +
  xlab("Method") +
  ylab("Min DE") +
  scale_y_discrete(limits = break_set, breaks = break_set) +
  scale_x_discrete(limits = c("CHOIR", "CIDR", "Cytocipher","dropClust", "GiniClust3", "PanoView", "RaceID3", "SAFEclustering", "SC3",
                              "SCCAF", "scCAN", "scSHC_new", "Seurat", "SHARP", "SIMLR", "Spectrum"),
                   breaks = c("CHOIR", "CIDR", "Cytocipher","dropClust", "GiniClust3", "PanoView", "RaceID3", "SAFEclustering", "SC3",
                              "SCCAF", "scCAN", "scSHC_new", "Seurat", "SHARP", "SIMLR", "Spectrum")) +
  Seurat::NoLegend()
