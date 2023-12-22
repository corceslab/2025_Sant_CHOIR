# ---------------------------------------------------------------------------
# This script generates plots for the simulated datasets
# ---------------------------------------------------------------------------
library(pdfCluster)
library(ggbeeswarm)
library(dplyr)
library(ggplot2)
library(Seurat)
library(pbmcapply)
library(CHOIR)

# Get max values
groups <- 1
size <- 500

max_maxvmem <- 0
max_mem <- 0
max_vmem <- 0
max_io <- 0
max_time <- 0
max_umap <- 0

for (i in c(1,2,3,4,5)) {
  print(i)
  input_dir <- paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "/sim_", size, "_", i)
  if (i %in% prob) {
    input_dir <- paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splat_", groups, "/sim_", size, "_", i)
  }
  cluster_parameter_list <- read.csv(paste0(input_dir, "/intermediate_files/clusters/compiled_clusters_metrics.csv"))

  # Max values
  max_maxvmem <- max(max_maxvmem, cluster_parameter_list$maxvmem, na.rm = TRUE)
  max_mem <- max(max_mem, cluster_parameter_list$mem, na.rm = TRUE)
  max_vmem <- max(max_vmem, cluster_parameter_list$vmem, na.rm = TRUE)
  max_io <- max(max_io, cluster_parameter_list$io, na.rm = TRUE)
  max_time <- max(max_time, cluster_parameter_list$time, na.rm = TRUE)
  
  # Max UMAP values
  # CHOIR dim reduction
  choir_data <- readRDS(paste0(input_dir, "/intermediate_files/clusters/CHOIR_parameter_list_RNA_simulated/CHOIR_output_1.rds"))
  umap_coords <- RunUMAP(choir_data$reduction$P0_reduction)
  saveRDS(umap_coords, paste0(input_dir, "/umap_coords.rds"))
  max_umap <- max(max_umap, abs(c(umap_coords@cell.embeddings[,1], umap_coords@cell.embeddings[,2])))
}

# Plots
for (iteration in c(1,2,3,4,5)) {
  print(iteration)
  # Import 
  input_dir <- paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", groups, "/sim_", size, "_", iteration)
  if (iteration %in% prob) {
    input_dir <- paste0("/Volumes/Mucke-Sequencing/Cathrine/cluster_benchmarking/datasets/Splat_", groups, "/sim_", size, "_", iteration)
  }
  cluster_parameter_list <- read.csv(paste0(input_dir, "/intermediate_files/clusters/compiled_clusters_metrics.csv"))
  
  # Zero out memory and time if complete and not measured
  cluster_parameter_list <- cluster_parameter_list %>% 
    mutate(maxvmem = ifelse(status == "complete" & is.na(maxvmem), 0, maxvmem),
           vmem = ifelse(status == "complete" & is.na(vmem), 0, vmem),
           mem = ifelse(status == "complete" & is.na(mem), 0, mem),
           io = ifelse(status == "complete" & is.na(io), 0, io),
           time = ifelse(status == "complete" & is.na(time), 0, time))
  
  # UMAP
  # cell metadata
  cell_metadata <- read.csv(paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata.csv"))
  rownames(cell_metadata) <- cell_metadata$CellID
  # CHOIR dim reduction
  umap_coords <- readRDS(paste0(input_dir, "/umap_coords.rds"))
  # temporary Seurat object
  tmp <- matrix(stats::rnorm(nrow(umap_coords) * 3, 10), ncol = nrow(umap_coords), nrow = 3)
  colnames(tmp) <- rownames(umap_coords)
  rownames(tmp) <- paste0("t",seq_len(nrow(tmp)))
  tmp_seurat <- Seurat::CreateSeuratObject(tmp, min.cells = 0, min.features = 0, assay = 'tmp', meta.data = cell_metadata)
  tmp_seurat[["umap"]] <- Seurat::CreateDimReducObject(embeddings = umap_coords@cell.embeddings, assay = 'tmp')

  if (max_umap < 3) {
    umap_lims = c(-3,3)
    umap_breaks = c(-3,-2,-1,0,1,2,3)
  } else if (max_umap < 4) {
    umap_lims = c(-4,4)
    umap_breaks = c(-4,-3,-2,-1,0,1,2,3,4)
  } else if (max_umap < 15) {
    umap_lims = c(-15,15)
    umap_breaks = c(-15,-10,-5,0,5,10,15)
  } else if (max_umap < 20) {
    umap_lims = c(-20,20)
    umap_breaks = c(-20,-15,-10,-5,0,5,10,15,20)
  } else if (max_umap < 25) {
    umap_lims = c(-25,25)
    umap_breaks = c(-25,-20,-15,-10,-5,0,5,10,15,20,25)
  }
  # UMAP labeled by ground truth groups
  tmp_seurat@meta.data$Ground_truth2 <- as.numeric(sub("Group", "", tmp_seurat@meta.data$Ground_truth))
  umap_plot <- DimPlot(tmp_seurat, group.by = "Ground_truth2", label = TRUE) +
    theme(plot.title = element_text(size = rel(0))) +
    scale_color_manual(values = CHOIRpalette(groups)[c(4,5,3,1,2)]) +
    scale_x_continuous(breaks = umap_breaks, limits = umap_lims) +
    scale_y_continuous(breaks = umap_breaks, limits = umap_lims) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    NoLegend()
  umap_plot[[1]]$layers[[1]]$aes_params$alpha = .2
  umap_plot
  ggsave(paste0("~/Desktop/CHOIR_simulated_plots/plot_", groups, "_", size, "_", iteration, "_umap.pdf"), width = 5, height = 5, units = "in")

  # UMAP labeled with default CHOIR clusters
  CHOIR_clusters <- read.csv(paste0(input_dir, "/intermediate_files/clusters/compiled_clusters_CHOIR_parameter_list_RNA_simulated.csv"))
  rownames(CHOIR_clusters) <- CHOIR_clusters$CellID
  tmp_seurat@meta.data$CHOIR_clusters <- CHOIR_clusters[colnames(tmp_seurat), "Parameters_1"]
  umap_plot <- DimPlot(tmp_seurat, group.by = "CHOIR_clusters", label = TRUE) +
    theme(plot.title = element_text(size = rel(0))) +
    scale_color_manual(values = CHOIRpalette(dplyr::n_distinct(tmp_seurat@meta.data$CHOIR_clusters))) +
    scale_x_continuous(breaks = umap_breaks, limits = umap_lims) +
    scale_y_continuous(breaks = umap_breaks, limits = umap_lims) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    NoLegend()
  umap_plot[[1]]$layers[[1]]$aes_params$alpha = .2
  umap_plot
  ggsave(paste0("~/Desktop/CHOIR_simulated_plots/plot_", groups, "_", size, "_", iteration, "_umap_CHOIR_1.pdf"), width = 5, height = 5, units = "in")
  
  # UMAP labeled with default Seurat clusters
  Seurat_clusters <- read.csv(paste0(input_dir, "/intermediate_files/clusters/compiled_clusters_cluster_parameter_list_Seurat_RNA_simulated_mat.csv"))
  rownames(Seurat_clusters) <- Seurat_clusters$CellID
  tmp_seurat@meta.data$Seurat_clusters <- Seurat_clusters[colnames(tmp_seurat), "Parameters_1"]
  umap_plot <- DimPlot(tmp_seurat, group.by = "Seurat_clusters", label = TRUE) +
    theme(plot.title = element_text(size = rel(0))) +
    scale_color_manual(values = CHOIRpalette(dplyr::n_distinct(tmp_seurat@meta.data$Seurat_clusters))) +
    scale_x_continuous(breaks = umap_breaks, limits = umap_lims) +
    scale_y_continuous(breaks = umap_breaks, limits = umap_lims) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    NoLegend()
  umap_plot[[1]]$layers[[1]]$aes_params$alpha = .2
  umap_plot
  ggsave(paste0("~/Desktop/CHOIR_simulated_plots/plot_", groups, "_", size, "_", iteration, "_umap_Seurat_1.pdf"), width = 5, height = 5, units = "in")
  
  # Time
  if (max_time/3600 < 0.1) {
    time_limits <- c(0,0.1)
    time_breaks <- c(0,0.02,0.04,0.06,0.08,0.1)
  } else if (max_time/3600 < 0.3) {
    time_limits <- c(0,0.3)
    time_breaks <- c(0,0.1,0.2,0.3)
  } else if (max_time/3600 < 0.7) {
    time_limits <- c(0,0.7)
    time_breaks <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7)
  } else if (max_time/3600 < 3) {
    time_limits <- c(0,3)
    time_breaks <- c(0,1,2,3)
  } else if (max_time/3600 < 15) {
    time_limits <- c(0,15)
    time_breaks <- c(0,5,10,15)
  } else if (max_time/3600 < 20) {
    time_limits <- c(0,20)
    time_breaks <- c(0,5,10,15,20)
  } else if (max_time/3600 < 72) {
    time_limits <- c(0,72)
    time_breaks <- c(0,24,48,72)
  } else {
    time_limits <- c(0,96)
    time_breaks <- c(0,24,48,72,96)
  }
  cluster_parameter_list %>% 
    dplyr::filter(status == "complete") %>%
    ggplot(aes(x = method, y = time/3600, alpha = default)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    geom_beeswarm(shape = 21, fill = "white", cex = 0.25) + 
    geom_point(data = . %>% filter(default == TRUE), shape = 23, alpha = 1, fill = "gold") +
    xlab("Method") +
    ylab("Time (h)") + 
    scale_y_continuous(limits = time_limits, breaks = time_breaks) +
    scale_x_discrete(limits = c("CHOIR", "CIDR", "Cytocipher","dropClust", "GiniClust3", "PanoView", "RaceID3", "SC3",
                                "SCCAF", "scCAN", "scSHC", "Seurat", "SHARP", "SIMLR", "Spectrum"),
                     breaks = c("CHOIR", "CIDR", "Cytocipher","dropClust", "GiniClust3", "PanoView", "RaceID3", "SC3",
                                "SCCAF", "scCAN", "scSHC", "Seurat", "SHARP", "SIMLR", "Spectrum")) +
    NoLegend() +
    scale_alpha_manual(values = c(1,0))
  ggsave(paste0("~/Desktop/CHOIR_simulated_plots/plot_", groups, "_", size, "_", iteration, "_time.pdf"), width = 5, height = 4, units = "in")
  
  # Number of clusters
  if (groups %in% c(1,5)) {
    num_limits = c(0,12)
    num_breaks = c(seq(0,10), 12)
    num_labels = c(seq(0,10), ">10")
    cap_max = 10
    cap_loc = 12
  } else if (groups == 10) {
    num_limits = c(0,22)
    num_breaks = c(seq(0,20), 22)
    num_labels = c(seq(0,20), ">20")
    cap_max = 20
    cap_loc = 22
  } else if (groups == 20) {
    num_limits = c(0,42)
    num_breaks = c(seq(0,40,5), 42)
    num_labels = c(seq(0,40,5), ">40")
    cap_max = 40
    cap_loc = 42
  }
  cluster_parameter_list %>%
    dplyr::filter(status == "complete") %>%
    mutate(n_clusters_capped = ifelse(n_clusters > cap_max, cap_loc, n_clusters)) %>%
    ggplot(aes(x = method, y = n_clusters_capped, alpha = default)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    geom_hline(yintercept = groups, color = "red") +
    geom_beeswarm(shape = 21, fill = "white", cex = 0.25) + 
    geom_point(data = . %>% filter(default == TRUE), shape = 23, alpha = 1, fill = "gold") +
    xlab("Method") +
    ylab("Number of clusters") +
    scale_x_discrete(limits = c("CHOIR", "CIDR", "Cytocipher","dropClust", "GiniClust3", "PanoView", "RaceID3", "SC3",
                                "SCCAF", "scCAN", "scSHC", "Seurat", "SHARP", "SIMLR", "Spectrum"),
                     breaks = c("CHOIR", "CIDR", "Cytocipher","dropClust", "GiniClust3", "PanoView", "RaceID3", "SC3",
                                "SCCAF", "scCAN", "scSHC", "Seurat", "SHARP", "SIMLR", "Spectrum")) +
    NoLegend() +
    scale_alpha_manual(values = c(1,0)) +
    scale_y_continuous(breaks = num_breaks, labels = num_labels, limits = num_limits)
  ggsave(paste0("~/Desktop/CHOIR_simulated_plots/plot_", groups, "_", size, "_", iteration, "_n_clusters.pdf"), width = 5, height = 4, units = "in")

  # ARI
  if (groups > 1) {
    if (min(cluster_parameter_list$ARI, na.rm = TRUE) < -0.1) {
      stop("ARI min")
    }
    cluster_parameter_list %>% 
      dplyr::filter(status == "complete") %>%
      ggplot(aes(x = method, y = ARI, alpha = default, color = correct_n_clusters)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      geom_beeswarm(shape = 21, fill = "white", cex = 0.25) + 
      geom_point(data = . %>% filter(default == TRUE), shape = 23, alpha = 1, fill = "gold") +
      xlab("Method") +
      ylab("Adjusted Rand Index") +
      scale_x_discrete(limits = c("CHOIR", "CIDR", "Cytocipher","dropClust", "GiniClust3", "PanoView", "RaceID3", "SC3",
                                  "SCCAF", "scCAN", "scSHC", "Seurat", "SHARP", "SIMLR", "Spectrum"),
                       breaks = c("CHOIR", "CIDR", "Cytocipher","dropClust", "GiniClust3", "PanoView", "RaceID3", "SC3",
                                  "SCCAF", "scCAN", "scSHC", "Seurat", "SHARP", "SIMLR", "Spectrum")) +
      NoLegend() +
      scale_alpha_manual(values = c(1,0)) +
      scale_color_manual(values = c("dark grey", "black")) +
      scale_y_continuous(breaks = seq(0,1, 0.2), limits = c(-0.1,1))
    ggsave(paste0("~/Desktop/CHOIR_simulated_plots/plot_", groups, "_", size, "_", iteration, "_ARI.pdf"), width = 5, height = 4, units = "in")
  }
}

