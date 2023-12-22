# ---------------------------------------------------------------------------
# Quality control for sn/scRNA-seq data
# ---------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(Matrix)
library(Seurat)
library(DropletUtils)
library(SingleCellExperiment)
library(scater)
library(flexmix)
library(BiocParallel)
library(scuttle)
library(DoubletFinder)

# Set up ---------------------------

args <- commandArgs(trailingOnly = TRUE)

# Define args
input_dir <- args[1]
sample_id <- args[2]
strategy <- args[3]
tech <- args[4]
cells_nuclei <- args[5]
cellranger <- args[6]
species <- args[7]
num_loaded <- args[8]
cell_metadata_avail <- args[9]

# Import read count matrix
mat_files <- list.files(path = paste0(input_dir, "/aligned_reads/", sample_id), pattern = "^matrix")
if ("matrix.mtx" %in% mat_files) {
  matrix_data <- readMM(paste0(input_dir, "/aligned_reads/", sample_id, "/matrix.mtx"))
} else {
  matrix_data <- read.table(paste0(input_dir, "/aligned_reads/", sample_id, "/matrix.tsv"))
}
barcodes <- read.table(file = paste0(input_dir, "/aligned_reads/",
                                     sample_id, "/barcodes.tsv"),
                       sep = "\t", header = FALSE)
features <- read.table(file = paste0(input_dir, "/aligned_reads/",
                                     sample_id, "/features.tsv"),
                       sep = "\t", header = FALSE)

if (ncol(features) == 1) {
  features$V2 <- features$V1
} else if (ncol(features) == 3) {
  features$V2 <- features$V3
}

if (nrow(matrix_data) == nrow(barcodes)) {
  matrix_data <- t(matrix_data)
}

# Remove peaks if multi-omic
if (tech == "RNA_ATAC") {
  subset_RNA <- features$V3 == "Gene Expression"
  features <- features[subset_RNA,]
  matrix_data <- matrix_data[subset_RNA,]
}

rownames(matrix_data) <- features$V1
colnames(matrix_data) <- barcodes$V1

# Empty droplets & ambient RNA ---------------------------

# Number of UMIs per droplet
n_UMIs <- colSums(matrix_data)
# Remove droplets with 0 UMIs
matrix_data <- matrix_data[,n_UMIs > 0]
n_UMIs <- n_UMIs[n_UMIs > 0]
# Check how many cells have <100 UMIs
less_than_100 <- sum(n_UMIs < 100)

# If data has been through cellranger:
if (cellranger == "yes") {
  # Import filtered barcodes
  filtered_barcodes <- read.table(file = paste0(input_dir, "/aligned_reads/",
                                                sample_id, "/filtered_barcodes.tsv"),
                                  sep = "\t", header = FALSE)$V1
  # Import cellranger clusters
  cellranger_clusters <- read.csv(paste0(input_dir, "/aligned_reads/",
                                         sample_id, "/cellranger_clusters.csv"))
  # Construct filtered count matrix
  filtered_counts <- matrix_data[,filtered_barcodes]
  # Add sample IDs to barcodes
  colnames(filtered_counts) <- paste(sample_id, colnames(filtered_counts), sep = "#")
} else if (less_than_100 >= 10) { # If there are at least 10 droplets with <100 UMIs
  print("Detect empty droplets with emptyDrops.")
  # emptyDrops method
  emptyDrops_check <- emptyDrops(matrix_data, lower = 100, test.ambient = FALSE)
  is_cell <- emptyDrops_check$FDR <= 0.01
  is_cell[is.na(is_cell)] <- FALSE
  filtered_counts <- matrix_data[,is_cell]
  filtered_counts <- as.matrix(filtered_counts)
  # Add sample IDs to barcodes
  colnames(filtered_counts) <- paste(sample_id, colnames(filtered_counts), sep = "#")
} else {
  # Assume that empty drops have already been filtered out
  print("Less than 10 droplets with <100 UMIs.")
  filtered_counts <- matrix_data
  # Add sample IDs to barcodes
  colnames(filtered_counts) <- paste(sample_id, colnames(filtered_counts), sep = "#")
}

# Low quality cells & doublets ---------------------------

# Create metadata dataframe
cell_metadata <- data.frame(CellID = colnames(filtered_counts))
if (cell_metadata_avail == "yes") {
  input_metadata <- read.csv(paste0(input_dir, "/aligned_reads/", sample_id, "/cell_metadata.csv"))
  input_metadata$CellID <- paste(sample_id, input_metadata$CellID, sep = "#")
  cell_metadata <- merge(cell_metadata, input_metadata, by = "CellID", all.x = TRUE)
}
if (!("Batch" %in% colnames(cell_metadata))) {
  cell_metadata$Batch <- sample_id
}

# Create SingleCellExperiment object
sce <- SingleCellExperiment(list(counts = filtered_counts))

# Get quality control metrics
# Mitochondrial genes
if (species == "human") {
  mt_genes <- dplyr::filter(features, grepl("^MT-",  V2))$V1
} else if (species == "mouse") {
  mt_genes <- dplyr::filter(features, grepl("^mt-",  V2))$V1
}
mt_genes_index <- rownames(sce) %in% mt_genes
features_mt <- list(mito = rownames(sce)[mt_genes_index])

# Identify # of recovered cells
num_recovered <- ncol(filtered_counts)
# Estimated multiplet rate
if (num_loaded == "unknown") {
  multiplet_rate = (num_recovered/1250)/100 # based on 10x estimates
} else {
  num_loaded <- as.numeric(num_loaded)
  multiplet_rate = (num_loaded/2062.5)/100 # based on 10x estimates
}

# Run QC
if (strategy != "ground_truth") {
  qc_metrics <- perCellQCMetrics(sce@assays@data[[1]], subsets = features_mt, BPPARAM = MulticoreParam(workers = 12))
  sce$total <- qc_metrics$sum
  sce$detected <- qc_metrics$detected
  sce$subsets_mito_percent <- qc_metrics$subsets_mito_percent
  # Plot QC features
  QC_p1 <- sce@colData %>% data.frame() %>%
    ggplot(aes(x = total/1e6)) +
    geom_histogram(bins = 100) +
    xlab("Library sizes (millions)") +
    ylab("Number of cells")
  QC_p2 <- sce@colData %>% data.frame() %>%
    ggplot(aes(x = detected)) +
    geom_histogram(bins = 100) +
    xlab("Number of expressed genes") +
    ylab("Number of cells")
  QC_p3 <- sce@colData %>% data.frame() %>%
    ggplot(aes(x = subsets_mito_percent)) +
    geom_histogram(bins = 100) +
    xlab("Mitochondrial percent (%)") +
    ylab("Number of cells")
  QC_p_all <- QC_p1 / QC_p2 / QC_p3
  ggsave(paste0(input_dir, "/quality_checks/", sample_id, "_QC_features.jpg"))
  
  ## MAD (with 0.5% minimum difference)
  # Get QC outliers
  qc_stats <- data.frame(libsize_drop = isOutlier(sce$total, nmads=3, type="lower", log=TRUE),
                         feature_drop = isOutlier(sce$detected, nmads=3, type="lower", log=TRUE),
                         mito_drop = isOutlier(sce$subsets_mito_percent, type="higher", min.diff=0.5))
  qc_stats$quality <- Reduce("|", qc_stats[,colnames(qc_stats)!="discard"])
  qc_merge <- qc_stats %>% transmute(CellID = colnames(sce),
                                     quality = ifelse(quality == TRUE, "low quality", "intact"))
  # Add to metadata
  cell_metadata <- merge(cell_metadata, qc_merge, by = "CellID", all.x = TRUE)
  # Plot filter thresholds
  MAD_p1 <- plotColData(sce, y="total",
                        colour_by=I(qc_stats$libsize_drop)) + ylab("Library size") +
    theme(legend.position = "none")
  MAD_p2 <- plotColData(sce, y="detected",
                        colour_by=I(qc_stats$feature_drop)) + ylab("Number of expressed genes") +
    theme(legend.position = "none")
  MAD_p3 <- plotColData(sce, y="subsets_mito_percent",
                        colour_by=I(qc_stats$mito_drop)) + ylab("Mitochondrial percent (%)")
  MAD_p_all <- MAD_p1 + MAD_p2 + MAD_p3
  ggsave(plot = MAD_p_all, filename = paste0(input_dir, "/quality_checks/", sample_id, "_MAD_QC_filters.jpg"))
  
  
  low_quality <- cell_metadata$CellID[cell_metadata[,"quality"] == "low quality"]
  qc_matrix <- sce@assays@data[[1]][ , !(colnames(sce@assays@data[[1]]) %in% low_quality)]
  
  if (multiplet_rate != 0) {
    # Create Seurat object & run standard workflow
    seurat_obj <- CreateSeuratObject(qc_matrix, min.cells = 1, min.features = 0) %>%
      NormalizeData() %>%
      FindVariableFeatures() %>%
      ScaleData() %>%
      RunPCA(npcs = 30) %>%
      RunUMAP(reduction = "pca", dims = 1:30)
    ## pK Identification (no ground-truth) ------------------------------------------
    sweep.res.list <- paramSweep_v3(seurat_obj, PCs = 1:30, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    pK_val <- as.numeric(as.character(filter(bcmvn,BCmetric == max(BCmetric))$pK[1]))
    ## Homotypic Doublet # Estimate ----------------------------------------
    nExp_poi <- round(multiplet_rate*nrow(seurat_obj@meta.data))
    
    print(paste0("nrows: ", nrow(seurat_obj@meta.data)))
    print(paste0("multiplet_rate: ", multiplet_rate))
    print(paste0("nExp_poi: ", nExp_poi))
    print(paste0("pK_val: ", pK_val))
    
    ## Run DoubletFinder  -----------------------------------------------------------
    seurat_obj <- doubletFinder_v3(seurat_obj, PCs = 1:30, pN = 0.25, pK = pK_val,
                                   nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    
    ncol_metadata <- ncol(seurat_obj@meta.data)
    doublet_data <- data.frame(CellID = rownames(seurat_obj@meta.data),
                               doublet_prob = seurat_obj@meta.data[, ncol_metadata-1],
                               doublet_class = seurat_obj@meta.data[, ncol_metadata])
    colnames(doublet_data) <- c("CellID",
                                "doublet_prob",
                                "doublet_class")
    cell_metadata <- merge(cell_metadata, doublet_data, by = "CellID", all.x = TRUE)
    
    ggplot(doublet_data, aes(x = doublet_data[,2], fill = doublet_data[,3])) +
      theme_classic() +
      geom_histogram(bins = 50) +
      xlab("Doublet probability score") +
      ylab("Count") +
      labs(fill = "Doublet class:")
    ggsave(paste0(input_dir, "/quality_checks/", sample_id, "_doublet_score.jpg"))
    
    DimPlot(seurat_obj, group.by = colnames(seurat_obj@meta.data)[ncol_metadata])
    ggsave(paste0(input_dir, "/quality_checks/", sample_id, "_doublet_umap.jpg"))
  } else {
    doublet_data <- data.frame(CellID = cell_metadata$CellID,
                               doublet_class = rep("Singlet", nrow(cell_metadata)))
    colnames(doublet_data) <- c("CellID",
                                "doublet_class")
    cell_metadata <- merge(cell_metadata, doublet_data, by = "CellID", all.x = TRUE)
  }
}

# Output ---------------------------

# Save read counts
saveRDS(sce, paste0(input_dir, "/intermediate_files/qc/", sample_id, "/qc_matrix.rds"))

# Save cell metadata
saveRDS(cell_metadata, paste0(input_dir, "/intermediate_files/qc/", sample_id, "/qc_cell_metadata.rds"))
