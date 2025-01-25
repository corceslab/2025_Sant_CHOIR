# ---------------------------------------------------------------------------
# Run ArchR
# ---------------------------------------------------------------------------
library(ArchR)

run_ArchR <- function(input_dir,
                      cluster_parameter_index,
                      cluster_parameters,
                      temp_dir) {

  # Set up ---------------------------

  setwd(paste0(input_dir, "/intermediate_files/arrow_files"))

  # Data type
  input_data1 <- cluster_parameters$input_data1
  input_data2 <- cluster_parameters$input_data2

  # Parameters
  batch_correction_method <- cluster_parameters$parameter1_value
  algorithm <- as.numeric(cluster_parameters$parameter2_value)
  resolution <- as.numeric(cluster_parameters$parameter3_value)
  maxClusters <- as.numeric(cluster_parameters$parameter4_value)
  seed <- as.numeric(cluster_parameters$parameter5_value)
  n_cores <- as.integer(cluster_parameters$parameter6_value)

  # Input
  modality1 <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/", input_data1, ".rds"))
  # If multi-modal
  if (!is.na(input_data2)) {
    modality2 <- readRDS(paste0(input_dir, "/intermediate_files/preprocessed/", input_data2, ".rds"))
    multi_modal <- TRUE
    # Always ATAC as modality 1 and RNA as modality 2
    if (grepl('ATAC', input_data2)) {
      swap <- modality2
      modality2 <- modality1
      modality1 <- swap
    }
  } else {
    multi_modal <- FALSE
  }

  # ArchR ---------------------------

  # Prep

  # cell_metadata
  cell_metadata <- read.csv(paste0(input_dir, "/intermediate_files/preprocessed/cell_metadata.csv"))

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
                                    outputDirectory = paste0(temp_dir, "/ArchR_", cluster_parameter_index))
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
  } else {
    object <- modality1
    # Subset to cell_metadata
    object <- subsetArchRProject(ArchRProj = object,
                                 cells = cell_metadata$CellID,
                                 dropCells = FALSE,
                                 outputDirectory = paste0(temp_dir, "/choir_", cluster_parameter_index))
  }

  # Start time
  start_time <- Sys.time()

  # Dimensionality reduction
  if (multi_modal == FALSE) {
    # ATAC
    object <- addIterativeLSI(
      ArchRProj = object,
      saveIterations = FALSE,
      useMatrix = "TileMatrix",
      depthCol = "nFrags",
      name = "IterativeLSI",
      threads = n_cores,
      seed = seed)
  } else if (multi_modal == TRUE) {
    # ATAC
    object <- addIterativeLSI(
      ArchRProj = object,
      saveIterations = FALSE,
      useMatrix = "TileMatrix",
      depthCol = "nFrags",
      name = "LSI_ATAC",
      threads = n_cores,
      seed = seed)
    # RNA
    object <- addIterativeLSI(
      ArchRProj = object,
      saveIterations = FALSE,
      useMatrix = "GeneExpressionMatrix",
      depthCol = "Gex_nUMI",
      varFeatures = 2000,
      firstSelection = "variable",
      binarize = FALSE,
      name = "LSI_RNA",
      threads = n_cores,
      seed = seed)
    # Combined
    object <- addCombinedDims(object, name = "IterativeLSI", reducedDims = c("LSI_ATAC", "LSI_RNA"))
  }

  # Batch correction
  if (batch_correction_method == "Harmony") {
    object <- addHarmony(ArchRProj = object,
                       reducedDims = "IterativeLSI",
                       name = "Harmony",
                       groupBy = "Batch")
    use_dims <- "Harmony"
  } else {
    use_dims <- "IterativeLSI"
  }

  # Cluster
  try(object <- addClusters(object,
                        reducedDims = use_dims,
                        name = "Clusters",
                        resolution = resolution,
                        maxClusters = maxClusters,
                        force = TRUE,
                        seed = seed,
                        algorithm = algorithm))
  
  if ("Clusters" %in% colnames(object@cellColData)) {
    # proceed
  } else {
    # re-try -- bug w/ multi-modal addClusters (https://github.com/GreenleafLab/ArchR/issues/2011#issue-1865763760)
    # fixed with manual renaming of dimensionality reduction columns
    resolution = resolution
    algorithm = algorithm
    
    input = object
    reducedDims = use_dims
    name = "Clusters"
    sampleCells = NULL
    seed = seed
    method = "Seurat"
    dimsToUse = NULL
    scaleDims = NULL
    corCutOff = 0.75
    knnAssign = 10
    nOutlier = 5
    maxClusters = maxClusters
    testBias = TRUE
    filterBias = FALSE
    biasClusters = 0.01
    biasCol = "nFrags"
    biasVals = NULL
    biasQuantiles = c(0.05, 0.95)
    biasEnrich = 10
    biasProportion = 0.5
    biasPval = 0.05
    nPerm = 500
    prefix = "C"
    ArchRProj = NULL
    verbose = TRUE
    tstart = NULL
    force = TRUE
    logFile = createLogFile("addClusters")
    
    ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj", "null"))
    if (is(ArchRProj, "ArchRProject")) {
      message("When running addClusters 'input' param should be used for 'ArchRProj'. Replacing 'input' param with user 'ArchRPRoj'...")
      input <- ArchRProj
      rm(ArchRProj)
      gc()
    }
    ArchR:::.validInput(input = input, name = "input", valid = c("ArchRProj", "matrix"))
    ArchR:::.validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
    ArchR:::.validInput(input = name, name = "name", valid = c("character"))
    ArchR:::.validInput(input = sampleCells, name = "sampleCells", valid = c("integer", "null"))
    ArchR:::.validInput(input = seed, name = "seed", valid = c("integer"))
    ArchR:::.validInput(input = method, name = "method", valid = c("character"))
    ArchR:::.validInput(input = dimsToUse, name = "dimsToUse", valid = c("numeric", "null"))
    ArchR:::.validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
    ArchR:::.validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
    ArchR:::.validInput(input = knnAssign, name = "knnAssign", valid = c("integer"))
    ArchR:::.validInput(input = nOutlier, name = "nOutlier", valid = c("integer"))
    ArchR:::.validInput(input = testBias, name = "testBias", valid = c("boolean"))
    ArchR:::.validInput(input = filterBias, name = "filterBias", valid = c("boolean"))
    ArchR:::.validInput(input = biasClusters, name = "biasClusters", valid = c("numeric"))
    ArchR:::.validInput(input = biasCol, name = "biasCol", valid = c("character"))
    ArchR:::.validInput(input = biasQuantiles, name = "biasQuantiles", valid = c("numeric"))
    ArchR:::.validInput(input = biasEnrich, name = "biasEnrich", valid = c("numeric"))
    ArchR:::.validInput(input = biasProportion, name = "biasProportion", valid = c("numeric"))
    ArchR:::.validInput(input = biasPval, name = "biasPval", valid = c("numeric"))
    ArchR:::.validInput(input = nPerm, name = "nPerm", valid = c("integer"))
    ArchR:::.validInput(input = prefix, name = "prefix", valid = c("character"))
    ArchR:::.validInput(input = verbose, name = "verbose", valid = c("boolean"))
    ArchR:::.validInput(input = tstart, name = "tstart", valid = c("timestamp","null"))
    ArchR:::.validInput(input = force, name = "force", valid = c("boolean"))
    ArchR:::.validInput(input = logFile, name = "logFile", valid = c("character"))
    
    ArchR:::.startLogging(logFile = logFile)
    
    if (is.null(tstart)) {
      tstart <- Sys.time()
    }
    
    if (inherits(input, "ArchRProject")) {
      #Check
      input <- addCellColData(
        ArchRProj = input, 
        data = rep(NA, nCells(input)), 
        name = name, 
        cells = getCellNames(input), 
        force = force
      )
      
      if (reducedDims %ni% names(input@reducedDims)) {
        stop("Error reducedDims not available!")
      }
      
      matDR <- getReducedDims(
        ArchRProj = input, 
        reducedDims = reducedDims, 
        dimsToUse = dimsToUse, 
        corCutOff = corCutOff, 
        scaleDims = scaleDims
      )
      
    } else if(inherits(input, "matrix")) {
      matDR <- input
    } else {
      stop("Input an ArchRProject or Cell by Reduced Dims Matrix!")
    }
    
    #Subset Matrix
    set.seed(seed)
    nr <- nrow(matDR)
    
    if (!is.null(sampleCells)) {
      if (sampleCells < nrow(matDR)) {
        ArchR:::.logDiffTime("Estimating Clusters by Sampling", tstart, verbose = verbose, logFile = logFile)
        estimatingClusters <- 1
        idx <- sample(seq_len(nrow(matDR)), sampleCells)
        matDRAll <- matDR
        matDR <- matDR[idx,,drop=FALSE]
      } else {
        estimatingClusters <- 0
      }
    } else {
      estimatingClusters <- 0
    }
    
    #################################################################################
    # Decide on which clustering setup to use
    #################################################################################
    if (grepl("seurat",tolower(method))) {
    } else if (grepl("scran",tolower(method))) {
    } else {
      stop("Clustering Method Not Recognized!")
    }
    
    clustParams <- list("resolution" = resolution, "algorithm" = algorithm)
    clustParams$verbose <- verbose
    clustParams$tstart <- tstart
    #clust <- ArchR:::.clustSeurat(mat = matDR, clustParams = clustParams, logFile = logFile)
    
    mat = matDR
    clustParams = clustParams
    logFile = logFile
    
    #Simply a wrapper on Seurats FindClusters
    ArchR:::.requirePackage("Seurat", source = "cran")
    ArchR:::.logDiffTime("Running Seurats FindClusters (Stuart et al. Cell 2019)", clustParams$tstart, verbose=clustParams$verbose, logFile = logFile)
    
    tmp <- matrix(rnorm(nrow(mat) * 3, 10), ncol = nrow(mat), nrow = 3)
    colnames(tmp) <- rownames(mat)
    rownames(tmp) <- paste0("t",seq_len(nrow(tmp)))
    
    obj <- Seurat::CreateSeuratObject(tmp, project='scATAC', min.cells=0, min.features=0)
    identical(colnames(obj),rownames(mat))
    mat2 <- mat
    colnames(mat) <- paste0("LSI", 1:ncol(mat))
    obj[['pca']] <- Seurat::CreateDimReducObject(embeddings=mat, key='PC_', assay='RNA')
    clustParams$object <- obj
    clustParams$reduction <- "pca"
    clustParams$dims <- seq_len(ncol(mat))
    
    obj <- suppressWarnings(do.call(Seurat::FindNeighbors, clustParams))
    clustParams$object <- obj
    clustParams
    
    
    clust <- tryCatch({
      
      cS <- Matrix::colSums(obj@graphs$RNA_snn)
      
      if (cS[length(cS)] == 1) {
        
        #Error Handling with Singletons
        idxSingles <- which(cS == 1)
        idxNonSingles <- which(cS != 1)
        
        rn <- rownames(mat) #original order
        mat <- mat[c(idxSingles, idxNonSingles), ,drop = FALSE]
        
        tmp <- matrix(rnorm(nrow(mat) * 3, 10), ncol = nrow(mat), nrow = 3)
        colnames(tmp) <- rownames(mat)
        rownames(tmp) <- paste0("t",seq_len(nrow(tmp)))
        
        obj <- Seurat::CreateSeuratObject(tmp, project='scATAC', min.cells=0, min.features=0)
        obj[['pca']] <- Seurat::CreateDimReducObject(embeddings=mat, key='PC_', assay='RNA')
        clustParams$object <- obj
        clustParams$reduction <- "pca"
        clustParams$dims <- seq_len(ncol(mat))
        
        obj <- ArchR:::.suppressAll(do.call(Seurat::FindNeighbors, clustParams))
        clustParams$object <- obj
        
        obj <- suppressWarnings(do.call(Seurat::FindClusters, clustParams))
        
        #Get Output
        clust <- obj@meta.data[,ncol(obj@meta.data)]
        clust <- paste0("Cluster",match(clust, unique(clust)))
        names(clust) <- rownames(mat)
        clust <- clust[rn]
        
      } else {
        
        obj <- suppressWarnings(do.call(Seurat::FindClusters, clustParams))
        
        #Get Output
        clust <- obj@meta.data[,ncol(obj@meta.data)]
        clust <- paste0("Cluster",match(clust, unique(clust)))
        names(clust) <- rownames(mat)
        
      }
      
      clust
      
    }, error = function(e) {
      
      errorList <- append(args, mget(names(formals()),sys.frame(sys.nframe())))
      ArchR:::.logError(e, fn = "FindClusters", info = "", errorList = errorList, logFile = logFile)
      
    })
    
    
    
    #################################################################################
    # If estimating clsuters we will assign to nearest neighbor cluster
    #################################################################################
    if (estimatingClusters == 1) {
      
      ArchR:::.logDiffTime("Finding Nearest Clusters", tstart, verbose = verbose, logFile = logFile)
      knnAssigni <- as.matrix(.computeKNN(matDR, matDRAll[-idx,,drop=FALSE], knnAssign))
      clustUnique <- unique(clust)
      clustMatch <- match(clust, clustUnique)
      knnAssigni <- matrix(apply(knnAssigni, 2, function(x) clustMatch[x]), ncol = knnAssign)
      
      ArchR:::.logDiffTime("Assigning Nearest Clusters", tstart, verbose = verbose, logFile = logFile)
      clustAssign <- lapply(seq_along(clustUnique), function(x){
        rowSums(knnAssigni == x)
      }) %>% Reduce("cbind", .) %>% apply(., 1, which.max)
      clustOld <- clust
      clust <- rep(NA, nr)
      clust[idx] <- clustOld
      clust[-idx] <- clustUnique[clustAssign]
      matDR <- matDRAll
      remove(matDRAll)
      gc()
      
    }
    
    #################################################################################
    # Testing Bias
    #################################################################################
    if (testBias) {
      if (inherits(input, "ArchRProject")) {
        if (is.null(biasVals)) {
          biasDF <- getCellColData(input, select = biasCol)
        } else {
          biasDF <- DataFrame(row.names = rownames(matDR), bias = biasVals)
        }
      } else {
        if (!is.null(biasVals)) {
          biasDF <- DataFrame(row.names = rownames(matDR), bias = biasVals)
        } else {
          message("No biasVals for testing bias continuing without bias detection")
          testBias <- FALSE
        }
      }
    }
    
    if (testBias) {
      clust <- tryCatch({
        biasDF$Q <- ArchR:::.getQuantiles(biasDF[,1])
        tabClust <- table(clust)
        tabClustP <- tabClust / sum(tabClust)
        idxTest <- which(tabClustP < biasClusters)
        names(clust) <- rownames(matDR)
        if(length(idxTest) > 0){
          ArchR:::.logDiffTime("Testing Biased Clusters", tstart, verbose = verbose, logFile = logFile)
          testDF <- lapply(seq_along(idxTest), function(i){
            clustTesti <- names(tabClustP)[idxTest[i]]
            biasQ <- biasDF[names(clust)[which(clust == clustTesti)], 2]
            biasBgd <- matrix(
              sample(
                x = biasDF[names(clust)[which(clust != clustTesti)], 2],
                size = nPerm * length(biasQ),
                replace = if(nPerm * length(biasQ) > nrow(biasDF[names(clust)[which(clust != clustTesti)], ])) TRUE else FALSE
              ), 
              nrow = length(biasQ), 
              ncol = nPerm
            )
            n1 <- colSums(biasBgd >= max(biasQuantiles))
            n2 <- colSums(biasBgd <= min(biasQuantiles))
            pval1 <- max(sum(sum(biasQ >= max(biasQuantiles)) < n1) * 2, 1) / length(n1)
            pval2 <- max(sum(sum(biasQ <= min(biasQuantiles)) < n2) * 2, 1) / length(n2)
            enrich1 <- sum(biasQ >= max(biasQuantiles)) / max(median(n1), 1)
            enrich2 <- sum(biasQ <= min(biasQuantiles)) / max(median(n2), 1)
            per1 <- sum(biasQ >= max(biasQuantiles)) / length(biasQ)
            per2 <- sum(biasQ <= min(biasQuantiles)) / length(biasQ)
            if(enrich1 > enrich2){
              enrichClust <- enrich1
              enrichPval <- min(pval1, 1)
              enrichPer <- per1
            }else{
              enrichClust <- enrich2
              enrichPval <- min(pval2, 1)
              enrichPer <- per2
            }
            DataFrame(Cluster = clustTesti, enrichClust = enrichClust, enrichPval = enrichPval, enrichProportion = enrichPer)
          }) %>% Reduce("rbind", .)
          
          clustAssign <- testDF[which(testDF$enrichClust > biasEnrich & testDF$enrichProportion > biasProportion & testDF$enrichPval <= biasPval),1]
          if(length(clustAssign) > 0){
            if(filterBias){
              ArchR:::.logDiffTime(sprintf("Assigning Biased Clusters (n = %s) to Neighbors", length(clustAssign)), tstart, verbose = verbose, logFile = logFile)
              for(i in seq_along(clustAssign)){
                clusti <- clustAssign[i]
                idxi <- which(clust==clusti)
                knni <- .computeKNN(matDR[-idxi,,drop=FALSE], matDR[idxi,,drop=FALSE], knnAssign)
                clustf <- unlist(lapply(seq_len(nrow(knni)), function(x) names(sort(table(clust[-idxi][knni[x,]]),decreasing=TRUE)[1])))
                clust[idxi] <- clustf
              }
            }else{
              ArchR:::.logDiffTime(sprintf("Identified Biased Clusters (n = %s), set filterBias = TRUE to re-assign these cells: ", length(clustAssign)), tstart, verbose = verbose, logFile = logFile)
              message("Biased Clusters : ", appendLF = FALSE)
              for(i in seq_along(clustAssign)){
                message(clustAssign[i], " ", appendLF = FALSE)
              }
              message("")
            }
          }
        }
        clust
      }, error = function(e) {
        
        errorList <- list(
          idxTest = if(exists("testDF", inherits = FALSE)) fragx else "Error with idxTest!",
          biasDF = if(exists("testDF", inherits = FALSE)) fragx else "Error with biasDF!",
          testDF = if(exists("testDF", inherits = FALSE)) fragx else "Error with testDF!",
          clustAssign = if(exists("idf", inherits = FALSE)) fragx else "Error with clustAssign!"
        )
        
        ArchR:::.logError(e, fn = "testBias", info = "", errorList = errorList, logFile = logFile)
        
      })
      
    }
    
    #################################################################################
    # Test if clusters are outliers identified as cells with fewer than nOutlier
    #################################################################################
    ArchR:::.logDiffTime("Testing Outlier Clusters", tstart, verbose = verbose, logFile = logFile)
    tabClust <- table(clust)
    clustAssign <- which(tabClust < nOutlier)
    if (length(clustAssign) > 0) {
      ArchR:::.logDiffTime(sprintf("Assigning Outlier Clusters (n = %s, nOutlier < %s cells) to Neighbors", length(clustAssign), nOutlier), tstart, verbose = verbose, logFile = logFile)
      for(i in seq_along(clustAssign)){
        clusti <- names(clustAssign[i])
        idxi <- which(clust==clusti)
        knni <- ArchR:::.computeKNN(matDR[-idxi,], matDR[idxi,], knnAssign)
        clustf <- unlist(lapply(seq_len(nrow(knni)), function(x) names(sort(table(clust[-idxi][knni[x,]]),decreasing=TRUE)[1])))
        clust[idxi] <- clustf
      }
    }
    
    #################################################################################
    # Merging if more than maxClusters
    #################################################################################
    if (!is.null(maxClusters)) {
      if (length(unique(clust)) > maxClusters) {
        ArchR:::.logDiffTime(sprintf("Identified more clusters than maxClusters allowed (n = %s). Merging clusters to maxClusters (n = %s).\nIf this is not desired set maxClusters = NULL!", length(clustAssign), maxClusters), tstart, verbose = verbose, logFile = logFile)
        meanDR <- t(ArchR:::.groupMeans(t(matDR), clust))
        hc <- hclust(dist(as.matrix(meanDR)))
        ct <- cutree(hc, maxClusters)
        clust <- mapLabels(
          labels = clust, 
          oldLabels = names(ct), 
          newLabels = paste0(prefix, ct)
        )
      }
    }
    
    #################################################################################
    # Renaming Clusters based on Proximity in Reduced Dimensions
    #################################################################################
    ArchR:::.logDiffTime(sprintf("Assigning Cluster Names to %s Clusters", length(unique(clust))), tstart, verbose = verbose, logFile = logFile)
    
    if (length(unique(clust)) > 1) {
      
      meanDR <- t(ArchR:::.groupMeans(t(matDR), clust))
      hc <- hclust(dist(as.matrix(meanDR)))
      out <- mapLabels(
        labels = clust, 
        oldLabels = hc$labels[hc$order], 
        newLabels = paste0(prefix, seq_along(hc$labels))
      )
      
    } else {
      out <- rep(paste0(prefix, "1"), length(clust))
    }
    
    if (inherits(input, "ArchRProject")) {
      input <- ArchR:::.suppressAll(addCellColData(
        input, 
        data = out, 
        name = name, 
        cells = rownames(matDR),
        force = TRUE
      ))
      ArchR:::.logDiffTime("Finished addClusters", t1 = tstart, verbose = verbose, logFile = logFile)
    } else if (!inherits(input, "ArchRProject")) {
      ArchR:::.logDiffTime("Finished addClusters", t1 = tstart, verbose = verbose, logFile = logFile)
    }
    
    object <- input
  }

  # End time
  end_time <- Sys.time()

  # Output ---------------------------

  # Create cluster assignments dataframe
  clusters <- data.frame(CellID = rownames(object@cellColData),
                         Clusters = object@cellColData$Clusters)

  colnames(clusters) <- c("CellID", paste0("Parameters_", cluster_parameter_index))

  # Time
  time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Return cluster assignments dataframe
  return(list("clusters" = clusters,
              "time" = time))
}

