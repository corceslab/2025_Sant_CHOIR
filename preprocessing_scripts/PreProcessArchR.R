# ---------------------------------------------------------------------------
# Create ArchR object and run QC steps
# ---------------------------------------------------------------------------
library(ArchR)
library(parallel)
library(Seurat)

# Set up ---------------------------

args <- commandArgs(trailingOnly = TRUE)

# Define args
input_dir <- args[1]
species <- args[2]
cell_metadata_avail <- args[3]
strategy <- args[4]
sample_ids <- args[5:length(args)]

if (species == "human") {
  addArchRGenome("hg38")
} else if (species == "mouse") {
  addArchRGenome("mm10")
}

addArchRThreads(1)

setwd(paste0(input_dir, "/intermediate_files/arrow_files"))
ArrowFiles <- readRDS(paste0(input_dir, "/intermediate_files/arrow_files/arrow_files.rds"))

# Create ArchR project ---------------------------
proj <- ArchRProject(ArrowFiles, copyArrows = FALSE)

# Add metadata if it exists
if (cell_metadata_avail == "yes") {
  for (i in 1:length(sample_ids)) {
    sample_id <- sample_ids[i]
    input_metadata <- read.csv(paste0(input_dir, "/aligned_reads/", sample_id, "/cell_metadata.csv"))
    proj@cellColData$CellID <- rownames(proj@cellColData)
    proj@cellColData <- merge(proj@cellColData, input_metadata, by = "CellID", all.x = TRUE)
  }
} else {
  # Add CellID to metadata
  proj@cellColData$CellID <- rownames(proj@cellColData)
}
if (!("Batch" %in% colnames(proj@cellColData))) {
  proj@cellColData$Batch <- proj@cellColData$Sample
}

if (strategy == "ground_truth") {
  #Filter Cells
  proj <- proj[!is.na(proj$Ground_truth)]
} else {
  # QC plot
  plotFragmentSizes(ArchRProj = proj)
  ggsave(paste0(input_dir, "/quality_checks/ArchR_QC_fragments.pdf"))
  plotTSSEnrichment(ArchRProj = proj)
  ggsave(paste0(input_dir, "/quality_checks/ArchR_QC_TSS.pdf"))
  # Filter Cells
  proj <- proj[proj$TSSEnrichment > 2 & proj$nFrags > 2500]
  # Doublet Filtration
  proj <- addDoubletScores(proj)
  proj <- filterDoublets(proj)
}

# Save object
saveRDS(proj, paste0(input_dir, "/intermediate_files/preprocessed/ATAC_obj_ArchR.rds"))
# Save metadata
write.csv(proj@cellColData, paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata_ArchR.csv"), row.names = FALSE)

write.table(paste0(input_dir, "/intermediate_files/preprocessed/ATAC_obj_ArchR.rds"),
            file = paste0(input_dir, "/quality_checks/ATAC_preprocessed_file_list_ArchR.txt"))
