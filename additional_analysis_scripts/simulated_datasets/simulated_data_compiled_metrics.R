# ---------------------------------------------------------------------------
# This script compiles metrics for the 100 simulated datasets
# ---------------------------------------------------------------------------
library(superheat)
library(dplyr)

# Import & compile metrics
metrics_i <- read.csv(paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", 1,"/sim_", 
                             500, "_", 1, "/intermediate_files/clusters/compiled_clusters_metrics.csv")) %>%
  mutate(dataset = paste0("sim_", 1, "_", 500, "_", 1))

metrics <- data.frame(matrix(ncol = length(metrics_i), nrow = 0))
colnames(metrics) <- colnames(metrics_i)

# Adjust
sizes <- seq(500,2500,500)
sizes <- seq(5000,25000,5000)
groups <- 20

for (size in 1:5) {
  print(size)
  for (i in 1:5) {
    metrics_i <- read.csv(paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups,"/sim_", 
                                 sizes[size], "_", i, "/intermediate_files/clusters/compiled_clusters_metrics.csv")) %>%
      mutate(dataset = paste0("sim_", groups, "_", sizes[size], "_", i))
    metrics <- rbind(metrics, metrics_i)
    rm(metrics_i)
  }
}

# Zero out memory and time if complete and not measured
metrics <- metrics %>% 
  mutate(maxvmem = ifelse(status == "complete" & is.na(maxvmem), 0, maxvmem),
         vmem = ifelse(status == "complete" & is.na(vmem), 0, vmem),
         mem = ifelse(status == "complete" & is.na(mem), 0, mem),
         io = ifelse(status == "complete" & is.na(io), 0, io),
         time = ifelse(status == "complete" & is.na(time), 0, time))
# Max out time if out of time
metrics <- metrics %>% 
  mutate(status = ifelse((is.na(time) & is.na(total_cpu_time)), "out of time", status),
         time = ifelse((is.na(time) & is.na(total_cpu_time)), 96*3600, time))

# Compile performance
metrics <- metrics %>% 
  mutate(both_correct = ifelse(((grepl("sim_1_", dataset) & 
                                   correct_n_clusters == TRUE) | 
                                  (correct_n_clusters == TRUE & ARI_0.9 == TRUE)), TRUE, FALSE),
         groups = as.numeric(sub("_","", substr(dataset, 5,6))))

best_indices <- metrics %>% group_by(method, list, index, default) %>% summarise(sum_correct = sum(both_correct, na.rm = TRUE)) %>%
  group_by(method, list) %>% arrange(method, -default) %>% slice_max(sum_correct, n = 1, with_ties = FALSE) %>% data.frame()
metrics$best <- NA
for (i in 1:nrow(metrics)) {
  metrics$best[i] <- ifelse(metrics$index[i] == best_indices[best_indices$method == metrics$method[i], "index"], TRUE, FALSE)
}

consolidated_metrics <- metrics %>% group_by(method, groups) %>% summarise(n = n()) %>% data.frame()
consolidated_metrics$best_sum_correct <- metrics %>% dplyr::filter(best == TRUE) %>% group_by(method, groups) %>% 
  summarise(best_sum_correct = sum(both_correct, na.rm = TRUE)) %>% data.frame() %>%
  dplyr::select(best_sum_correct) %>% unlist()
consolidated_metrics$default_sum_correct <- metrics %>% dplyr::filter(default == TRUE) %>% group_by(method, groups) %>% 
  summarise(default_sum_correct = sum(both_correct, na.rm = TRUE)) %>% data.frame() %>%
  dplyr::select(default_sum_correct) %>% unlist()

# Plot
superheat(consolidated_metrics[,-c(1,3)]/100,
          heat.pal = c("white",
                       "#FFFC5A","#FFD32B","#FF9965",
                       "#FD619D","#C732D5","#6E34FC",
                       "#1632FB","#021EA9", "#000436"),
          heat.pal.values = c(0, 0.0625,
                              0.125, 0.25, 0.375,
                              0.5, 0.625, 0.75, 
                              0.875, 0.9375, 1))