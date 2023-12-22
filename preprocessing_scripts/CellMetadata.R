# ---------------------------------------------------------------------------
# Finalize cell metadata
# ---------------------------------------------------------------------------

library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

# Define args
input_dir <- args[1]
tech <- args[2]

if (tech == "RNA") {
  cell_metadata <- read.csv(paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata_RNA.csv"))
  write.csv(cell_metadata, paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata.csv"), row.names = FALSE)
} else if (tech == "ATAC") {
  cell_metadata_ArchR <- read.csv(paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata_ArchR.csv"))
  write.csv(cell_metadata_ArchR, paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata.csv"), row.names = FALSE)
} else if (tech == "RNA_ATAC") {
  cell_metadata_RNA <- read.csv(paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata_RNA.csv"))
  cell_metadata_ArchR <- read.csv(paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata_ArchR.csv"))
  shared_cells <- intersect(cell_metadata_ArchR$CellID, cell_metadata_RNA$CellID)
  cell_metadata_RNA <- dplyr::filter(cell_metadata_RNA, CellID %in% shared_cells)
  cell_metadata_ArchR <-  dplyr::filter(cell_metadata_ArchR, CellID %in% shared_cells)
  cell_metadata <- merge(cell_metadata_ArchR, cell_metadata_RNA, by = c("CellID", "Batch"), all = TRUE)
  write.csv(cell_metadata, paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata.csv"), row.names = FALSE)
}
