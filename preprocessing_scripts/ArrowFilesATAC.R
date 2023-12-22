# ---------------------------------------------------------------------------
# Create Arrow files
# ---------------------------------------------------------------------------
library(ArchR)
library(parallel)

# Set up ---------------------------

args <- commandArgs(trailingOnly = TRUE)

# Define args
input_dir <- args[1]
species <- args[2]
sample_ids <- args[3:length(args)]

if (species == "human") {
  addArchRGenome("hg38")
} else if (species == "mouse") {
  addArchRGenome("mm10")
}

setwd(paste0(input_dir, "/intermediate_files/arrow_files"))

addArchRThreads(1)

#Get Input Fragment Files
inputFiles <- getInputFiles(paste0(input_dir, "/aligned_reads/", sample_ids))
print(inputFiles)
names(inputFiles) <- sample_ids

# Create Arrow Files ---------------------------
ArrowFiles <- createArrowFiles(inputFiles = inputFiles,
                               sampleNames = names(inputFiles),
                               minTSS = 1,
                               minFrags = 1000,
                               addTileMat = TRUE,
                               addGeneScoreMat = TRUE)

saveRDS(ArrowFiles, paste0(input_dir, "/intermediate_files/arrow_files/arrow_files.rds"))
