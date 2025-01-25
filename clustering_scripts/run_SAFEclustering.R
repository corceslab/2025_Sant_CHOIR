# ---------------------------------------------------------------------------
# Run SAFE-clustering
# ---------------------------------------------------------------------------
library(Seurat)
library(stringr)
library(SAFEclustering)

# Modified SAFE-clustering function to fix small bug (see https://github.com/yycunc/SAFEclustering/issues/6)
# Original: "cluster_number <- c(cluster_number, max(c(sc3OUTPUT)))"
# Fixed:    "cluster_number <- c(cluster_number, max(as.numeric(sc3OUTPUT)))"
individual_clustering <- function(inputTags, mt_filter = TRUE, mt.pattern = "^MT-", mt.cutoff = 0.1,
                                  SC3 = TRUE, gene_filter = FALSE, svm_num_cells = 5000, CIDR = TRUE, nPC.cidr = NULL,
                                  Seurat = TRUE, nGene_filter = TRUE, low.genes = 200, high.genes = 8000, nPC.seurat = NULL, resolution = 0.7,
                                  tSNE = TRUE, saver = FALSE, dimensions = 3, perplexity = 30, tsne_min_cells = 200, tsne_min_perplexity = 10, var_genes = NULL,
                                  SEED = 1){

  cluster_number <- NULL
  cluster_results <- NULL
  inputTags = as.matrix(inputTags)

  # Filter out cells that have mitochondrial genes percentage over 5%
  if (mt_filter == TRUE){
    mito.genes <- grep(pattern = mt.pattern, x = rownames(x = inputTags), value = TRUE)
    percent.mito <- Matrix::colSums(inputTags[mito.genes, ])/Matrix::colSums(inputTags)
    inputTags <- inputTags[,which(percent.mito <= mt.cutoff)]
  }

  ##### SC3
  if(SC3 == TRUE){
    message("Performing SC3 clustering...")

    sc3OUTPUT <- SAFEclustering:::sc3_SAFE(inputTags = inputTags, gene_filter = gene_filter,
                                           svm_num_cells = svm_num_cells, SEED = SEED)
    cluster_results <- rbind(cluster_results, matrix(c(sc3OUTPUT), nrow = 1, byrow = TRUE))
    cluster_number <- c(cluster_number, max(as.numeric(sc3OUTPUT))) ##### FIXED ERROR #####
  }


  ##### CIDR
  if(CIDR == TRUE){
    message("Performing CIDR clustering...")

    cidrOUTPUT <- SAFEclustering:::cidr_SAFE(inputTags = inputTags, nPC.cidr = nPC.cidr, SEED = SEED)

    if(is.null(nPC.cidr)) {
      nPC.cidr <- cidrOUTPUT@nPC
    }

    cluster_results <- rbind(cluster_results, matrix(c(cidrOUTPUT@clusters), nrow = 1, byrow = TRUE))
    cluster_number <- c(cluster_number,  cidrOUTPUT@nCluster)
  }


  ##### Seurat
  if (Seurat == TRUE){
    message("Performing Seurat clustering...")

    if(is.null(nPC.seurat)) {
      nPC.seurat <- nPC.cidr
    }

    seurat_output <- SAFEclustering:::seurat_SAFE(inputTags = inputTags, nGene_filter = nGene_filter, low.genes = low.genes, high.genes = high.genes,
                                                  nPC.seurat = nPC.seurat, resolution = resolution, SEED = SEED)
    cluster_results <- rbind(cluster_results, matrix(c(seurat_output), nrow = 1, byrow = TRUE))
    cluster_number <- c(cluster_number, max(!is.na(seurat_output)))
  }


  ##### tSNE+kmeans
  if(tSNE == TRUE){
    message("Performing tSNE + k-means clustering...")

    ### Dimensionality reduction by Rtsne
    if(length(inputTags[1,]) < tsne_min_cells) {
      perplexity = tsne_min_perplexity
    }

    tsne_kmeansOUTPUT <- SAFEclustering:::tSNE_kmeans_SAFE(inputTags = inputTags, saver = saver, dimensions = dimensions,
                                                           perplexity = perplexity, k.min = 2, k.max = max(cluster_number), var_genes = var_genes, SEED = SEED)
    cluster_results = rbind(cluster_results, matrix(c(tsne_kmeansOUTPUT$cluster), nrow = 1, byrow = TRUE))
  }

  return(cluster_results)
}

run_SAFEclustering <- function(input_dir,
                               cluster_parameter_index,
                               cluster_parameters) {
  # Set up ---------------------------
  if (grepl("Subsetted", input_dir)) {
    groups <- 5
    size_it <- substr(input_dir, 90,97)
  } else if (grepl("trajectory", input_dir)) {
    groups <- 5
    size_it <- substr(input_dir, 93,100)
  } else if (grepl("Zhengmix4uneq", input_dir)) {
    groups <- 4
    size_it <- substr(input_dir, 111,150)
  } else if (grepl("Zhengmix4eq", input_dir)) {
    groups <- 4
    size_it <- substr(input_dir, 109,150)
  } else {
    id = sub("/gladstone/mucke/sequencing/Cathrine/cluster_benchmarking/datasets/Splatter_", "", input_dir)
    groups <- sub("/", "", substr(id, 1, 2))
    if (grepl("downsample", id)) {
      size_it <- substr(id, 25, 35)
    } else if (grepl("de_2", id)) {
      size_it <- substr(id, 12, nchar(id))
    } else if (grepl("lib", id)) {
      size_it <- substr(id, 11, nchar(id))
    } else {
      size_it <- sub(paste0(groups, "/sim"), "", id)
    }
  }

  setwd(paste0("/wynton/group/muckelab/user/cpetersen/cluster_benchmarking/snakemake_output_sim_", groups, size_it, "_safe"))

  # Data type
  input_data1 <- cluster_parameters$input_data1

  # Parameters
  runSC3 <- as.logical(cluster_parameters$parameter1_value)
  gene_filter <- as.logical(cluster_parameters$parameter2_value)
  svm_num_cells <- as.numeric(cluster_parameters$parameter3_value)
  runCIDR <- as.logical(cluster_parameters$parameter4_value)
  runSeurat <- as.logical(cluster_parameters$parameter5_value)
  resolution <- as.numeric(cluster_parameters$parameter6_value)
  runtSNE <- as.logical(cluster_parameters$parameter7_value)
  dimensions <- as.numeric(cluster_parameters$parameter8_value)
  perplexity <- as.numeric(cluster_parameters$parameter9_value)
  runMCLA <- as.logical(cluster_parameters$parameter10_value)
  runCSPA <- as.logical(cluster_parameters$parameter11_value)
  runHGPA <- as.logical(cluster_parameters$parameter12_value)
  seed <- as.integer(cluster_parameters$parameter13_value)

  if (runCIDR == FALSE) {
    nPC.seurat <- 30
  } else {
    nPC.seurat <- NULL
  }

  # Input
  object <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/", input_data1, ".rds"))

  # SAFE-clustering ---------------------------

  # Start time
  start_time <- Sys.time()

  try(output <- individual_clustering(inputTags = object,
                                  mt_filter = FALSE,
                                  SC3 = runSC3,
                                  gene_filter = gene_filter,
                                  svm_num_cells = svm_num_cells,
                                  CIDR = runCIDR,
                                  Seurat = runSeurat,
                                  nPC.seurat = nPC.seurat,
                                  nGene_filter = FALSE,
                                  resolution = resolution,
                                  tSNE = runtSNE,
                                  dimensions = dimensions,
                                  perplexity = perplexity,
                                  SEED = seed))

  if (exists("output")) {
    output <- apply(output, 2, as.numeric)

    try(final_output <- SAFE(cluster_results = output,
                         program.dir = paste0("/wynton/group/muckelab/user/cpetersen/cluster_benchmarking/snakemake_output_sim_", groups, size_it, "_safe"),
                         MCLA = runMCLA,
                         CSPA = runCSPA,
                         HGPA = runHGPA,
                         SEED = seed))
    if (exists("final_output")) {
      # End time
      end_time <- Sys.time()

      # Output ---------------------------

      # Create cluster assignments dataframe
      clusters <- data.frame(CellID = colnames(object),
                             Clusters = final_output$optimal_clustering)

      colnames(clusters) <- c("CellID", paste0("Parameters_", cluster_parameter_index))

      # Time
      time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    } else {
      # Output ---------------------------

      # Create cluster assignments dataframe
      clusters <- data.frame(CellID = colnames(object),
                             Clusters = rep(NA, ncol(object)))

      colnames(clusters) <- c("CellID", paste0("Parameters_", cluster_parameter_index))

      # Time
      time <- NA
    }
  } else {
    # Output ---------------------------

    # Create cluster assignments dataframe
    clusters <- data.frame(CellID = colnames(object),
                           Clusters = rep(NA, ncol(object)))

    colnames(clusters) <- c("CellID", paste0("Parameters_", cluster_parameter_index))

    # Time
    time <- NA
  }

  # Return cluster assignments dataframe
  return(list("clusters" = clusters,
              "time" = time))
}
