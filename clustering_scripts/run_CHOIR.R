# ---------------------------------------------------------------------------
# Run CHOIR
# ---------------------------------------------------------------------------
library(Seurat)
library(ArchR)
library(stringr)
library(CHOIR)
library(BPCells)
library(SeuratObject)
library(stringr)

run_CHOIR <- function(input_dir,
                      cluster_parameter_index,
                      cluster_parameters,
                      temp_dir,
                      parameter_key) {

  # Set up ---------------------------

  # Data type
  input_data1 <- cluster_parameters$input_data1
  input_data2 <- cluster_parameters$input_data2

  # Parameters

  alpha <- as.numeric(cluster_parameters$parameter1_value)
  min_accuracy <- as.numeric(cluster_parameters$parameter2_value)
  batch_correction_method <- cluster_parameters$parameter3_value
  sample_max <- ifelse(cluster_parameters$parameter4_value == "Inf", Inf, 
                       as.numeric(cluster_parameters$parameter4_value))
  downsampling_rate <- ifelse(cluster_parameters$parameter5_value == "auto", "auto",
                              as.numeric(cluster_parameters$parameter5_value))
  random_seed <- as.numeric(cluster_parameters$parameter6_value)
  n_cores <- as.numeric(cluster_parameters$parameter7_value)

  # Input
  if (grepl("bpcells", input_data1)) {
    modality1 <- open_matrix_dir(dir = paste0(input_dir, "/intermediate_files/preprocessed/", input_data1))
    # If multi-omic
    if (!is.na(input_data2)) {
      modality2 <- open_matrix_dir(dir = paste0(input_dir, "/intermediate_files/preprocessed/", input_data2))
      multi_modal <- TRUE
    } else {
      multi_modal <- FALSE
    }
  } else {
    modality1 <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/", input_data1, ".rds"))
    # If multi-omic
    if (!is.na(input_data2)) {
      modality2 <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/", input_data2, ".rds"))
      multi_modal <- TRUE
    } else {
      multi_modal <- FALSE
    }
  }
  
  message("Multi-modal? ", multi_modal)

  # Process parameters
  if (grepl("ATAC", input_data1)) {
    atac <- TRUE
  } else {
    atac <- FALSE
  }

  # Process multi-modal parameters
  if (multi_modal == TRUE) {
    if (grepl("ATAC", input_data2)) {
      atac <- c(atac, TRUE)
    } else {
      atac <- c(atac, FALSE)
    }
    batch_correction_method <- strsplit(batch_correction_method, "_")[[1]]
  }
  
  # CHOIR ---------------------------

  # Prep

  # cell_metadata
  cell_metadata <- read.csv(paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata.csv"))
  message("Read cell metadata.")

  # Depends on object type & tech
  if (any(grepl("ArchR", c(input_data1, input_data2)))) {
    setwd(paste0(input_dir, "/intermediate_files/arrow_files"))
    # Load ArchR object
    # If multi-modal
    if (multi_modal == TRUE) {
      # Subset to shared cell IDs
      # Find shared cell IDs
      cell_ids_1 <- rownames(modality1@cellColData)
      cell_ids_2 <- colnames(modality2)
      shared_cell_ids <- intersect(cell_ids_1, cell_ids_2)
      # Subset to cell_metadata
      shared_cell_ids <- intersect(shared_cell_ids, cell_metadata$CellID)
      # Subset
      modality1 <- subsetArchRProject(ArchRProj = modality1,
                                      cells = shared_cell_ids,
                                      dropCells = FALSE,
                                      outputDirectory = paste0(temp_dir, "/choir_", cluster_parameter_index))
      modality2 <- modality2[, shared_cell_ids]
      # Add RNA to ArchR object
      # Import RNA features
      metadata <- read.csv(paste0(input_dir, "/metadata.csv"))
      sample_ids <- dplyr::filter(metadata, variable == "sample_ids")$value
      sample_ids <- unlist(str_split(sample_ids, " "))
      RNA_features <- data.frame(V1 = NULL,
                                 V2 = NULL,
                                 V3 = NULL,
                                 V4 = NULL,
                                 V5 = NULL,
                                 V6 = NULL)
      for (i in 1:length(sample_ids)) {
        sample_id <- sample_ids[i]
        features_i <- read.table(file = paste0(input_dir, "/aligned_reads/",
                                               sample_id, "/features.tsv"),
                                 sep = "\t", header = FALSE)
        features_i <- features_i %>% dplyr::filter(V3 == "Gene Expression", V1 %in% rownames(modality2))
        RNA_features <- rbind(RNA_features, features_i)
      }
      RNA_features <- unique(RNA_features)
      RNA_features <- RNA_features[match(RNA_features$V1, rownames(modality2)),]
      if (identical(RNA_features$V1, rownames(modality2))) {
        RNA_features$V4[RNA_features$V4 == ""] <- " "
        rowRanges <- GRanges(RNA_features$V4,
                             IRanges(start = RNA_features$V5, end = RNA_features$V6),
                             feature_id = RNA_features$V1)
        seRNA <- SummarizedExperiment(assays = SimpleList(counts = modality2),
                                      rowRanges = rowRanges)
      } else {
        stop("Cannot create GRanges object because feature names do not match.")
      }
      object <- addGeneExpressionMatrix(input = modality1, seRNA = seRNA, force = TRUE)
      ArchR_matrix <- c("TileMatrix", "GeneExpressionMatrix")
      ArchR_depthcol <- c("nFrags", "Gex_nUMI")
    } else {
      object <- modality1
      # Subset to cell_metadata
      object <- subsetArchRProject(ArchRProj = object,
                                   cells = cell_metadata$CellID,
                                   dropCells = FALSE,
                                   outputDirectory = paste0(temp_dir, "/choir_", cluster_parameter_index))
      ArchR_matrix <- "TileMatrix"
      ArchR_depthcol <- "nFrags"
    }

    if (batch_correction_method  == "Harmony") {
      # Cell metadata
      rownames(cell_metadata) <- cell_metadata$CellID
      cell_metadata <- cell_metadata[rownames(object@cellColData), ]
      object@cellColData$Batch <- as.character(cell_metadata$Batch)
      batch_labels <- "Batch"
    } else {
      batch_labels <- NULL
    }
    use_assay <- NULL
    use_slot <- NULL
  } else {
    # Create Seurat object
    if (input_data1 == "bpcells") {
      # Create Seurat Object
      object <- CreateSeuratObject(counts = modality1, assay = "RNA")
      # Normalize
      object <- NormalizeData(object)
      use_assay <- "RNA"
      use_slot <- "data"
    } else {
      object <- CreateSeuratObject(modality1, assay = "RNA")
      use_assay <- "RNA"
      use_slot <- "counts"
    }
    message("Created Seurat object.")
    if (batch_correction_method  == "Harmony") {
      # Cell metadata
      rownames(cell_metadata) <- cell_metadata$CellID
      cell_metadata <- cell_metadata[colnames(object),]
      object@meta.data$Batch <- as.character(cell_metadata$Batch)
      batch_labels <- "Batch"
    } else {
      batch_labels <- NULL
    }
    ArchR_matrix <- NULL
    ArchR_depthcol <- NULL
    message("Set batch correction, etc.")
  }
  message("Created object.")

  # Start time
  start_time <- Sys.time()
  
  # Cluster
  object <- CHOIR(object = object,
                  alpha = alpha,
                  p_adjust = p_adjust,
                  min_accuracy = min_accuracy,
                  sample_max = sample_max,
                  downsampling_rate = downsampling_rate,
                  batch_correction_method = batch_correction_method,
                  batch_labels = batch_labels,
                  use_assay = use_assay,
                  use_slot = use_slot,
                  ArchR_matrix = ArchR_matrix,
                  ArchR_depthcol = ArchR_depthcol,
                  atac = atac,
                  n_cores = n_cores,
                  random_seed = random_seed,
                  verbose = TRUE)
  
  # End time
  end_time <- Sys.time()

  # Output ---------------------------

  if (grepl("ArchR", input_data1)) {
    # Save CHOIR output
    saveRDS(object@projectMetadata$CHOIR,
            paste0(input_dir, "/intermediate_files/clusters/", parameter_key, "/CHOIR_output_", cluster_parameter_index, ".rds"))
  } else {
    # Save CHOIR output
    saveRDS(object@misc$CHOIR,
            paste0(input_dir, "/intermediate_files/clusters/", parameter_key, "/CHOIR_output_", cluster_parameter_index, ".rds"))
  }

  # Create cluster assignments dataframe
  if (is.null(ArchR_matrix)) {
    clusters <- data.frame(CellID = rownames(object@meta.data),
                           Clusters = object@meta.data[, paste0("CHOIR_clusters_", alpha)])
  } else {
    clusters <- data.frame(CellID = rownames(object@cellColData),
                           Clusters = object@cellColData[, paste0("CHOIR_clusters_", alpha)])
  }

  colnames(clusters) <- c("CellID", paste0("Parameters_", cluster_parameter_index))

  # Time
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Return cluster assignments dataframe
  return(list("clusters" = clusters,
              "time" = time))
}

