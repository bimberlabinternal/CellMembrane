#' @title Run escape ssGSEA
#'
#' @description This will run escape to calcualte ssGSEA on Hallmark gene sets
#' @param seuratObj A Seurat object.
#' @param outputAssayBaseName The basename of the assay to store results. Escape will be run once for each gene set.
#' @param doPlot If true, a FeaturePlot will be printed for each pathway
#' @param performDimRedux If true, the standard seurat PCA/FindClusters/UMAP process will be run on the escape data. This may be most useful when using a customGeneSet or a smaller set of features/pathways
#' @param msigdbGeneSets A vector containing gene set codes specifying which gene sets should be fetched from MSigDB and calculated. Some recommendations in increasing computation time: H (hallmark, 50 gene sets), C8 (scRNASeq cell type markers, 830 gene sets), C2 (curated pathways, 6366 gene sets), GO:BP (GO biological processes, 7658). Each item will be processed into a separate assay.
#' @param customGeneSets A (preferably named) list containing gene sets to be scored by escape.
#' @param customGeneSetAssayName The name for the output assay (prefixed with outputAssayBaseName) for any customGeneSets
#' @param maxBatchSize If more than this many cells are in the object, it will be split into batches of this size and run in serial.
#' @param nCores Passed to runEscape()
#' @return The seurat object with results stored in an assay
#' @export
RunEscape <- function(seuratObj, outputAssayBaseName = "escape.", doPlot = FALSE, performDimRedux = FALSE, msigdbGeneSets = c("H", "C5" = "GO:BP", "C5" = "GO:MF"), customGeneSets = NULL, customGeneSetAssayName = 'CustomGeneSet', maxBatchSize = 100000, nCores = 1) {
  assayToGeneSets <- list()

  # NOTE: currently escape only supports RNA:
  assayName <- 'RNA'

  if (all(!is.null(customGeneSets), !(length(customGeneSets) == 0))) {
    if (!is.list(customGeneSets)){
      stop("customGeneSets is not a list. Please coerce it into a named list.")
    }

    if (length(customGeneSets) == 1) {
      stop('Due to issues with seurat5 single-row assays, customGeneSets must contain as least 2 gene sets')
    }

    # Check for feature existence. I think it's reasonable to expect someone might grab a new/non-msigdb human gene set and supply it to customGeneSets,
    # so a warning rather than stop for feature existence is sufficient.
    allGenes <- unique(unlist(customGeneSets))
    if (!(all(allGenes %in% rownames(seuratObj)))) {
      missingGenes <- allGenes[!(allGenes %in% rownames(seuratObj))]
      warning(paste0("Some genes in customGeneSets are not features in the Seurat object. These features are missing and will not contribute to scoring: ",
                     paste0(missingGenes, collapse = ', ')))
      
      if (length(allGenes) == length(missingGenes)) {
        stop('There are no genes shared between those provided in customGeneSets and the seurat object')
      }
    }

    # If the names are not provided in the custom gene sets, we should warn the user but we can label them here.
    if (is.null(names(customGeneSets))){
      names(customGeneSets) <- paste0("CustomGeneSet", seq_along(customGeneSets))
      warning("The customGeneSets list is unnamed. Naming the genesets with generic names such as: CustomGeneSet1, CustomGeneSet2 ...")
    }

    if (any(grepl(names(customGeneSets), pattern = '_'))) {
      print('Converting underscore to hyphen in set names')
      toFix <- grepl(names(customGeneSets), pattern = '_')
      names(customGeneSets)[toFix] <- gsub(names(customGeneSets)[toFix], pattern = '_', replacement = '-')
    }

    assayToGeneSets[[paste0(outputAssayBaseName, customGeneSetAssayName)]] <- customGeneSets
  }

  for (idx in seq_along(msigdbGeneSets)) {
    if (!is.null(names(msigdbGeneSets)[idx]) && names(msigdbGeneSets)[idx] != '') {
      libraryName <- names(msigdbGeneSets)[idx]
      subcategory <- msigdbGeneSets[[idx]]
      outputAssayName <- paste0(outputAssayBaseName, libraryName, '.', subcategory)
    } else {
      libraryName <- msigdbGeneSets[[idx]]
      subcategory <- NULL
      outputAssayName <- paste0(outputAssayBaseName, libraryName)
    }

    outputAssayName <- gsub(outputAssayName, pattern = ':', replacement = '.')

    print(paste0('Querying library: ', libraryName , ' / subcategory: ', subcategory))
    GS <- escape::getGeneSets(library = libraryName, subcategory = subcategory)

    assayToGeneSets[[outputAssayName]] <- GS
  }

  nBatches <- 1
  if (ncol(seuratObj) > maxBatchSize) {
    nBatches <- ceiling(ncol(seuratObj) / maxBatchSize)
    print(paste0('The object will be split into ', nBatches, ', batches'))
  }

  for (outputAssayName in names(assayToGeneSets)) {
    GS <- assayToGeneSets[[outputAssayName]]
    print(paste0('Processing ', outputAssayName, ' with ', length(GS), ' gene sets'))

    assayCounts <- NULL
    if (nBatches == 1) {
      assayCounts <- .RunEscapeOnSubset(seuratObj = seuratObj, outputAssayName = outputAssayName, GS = GS, nCores = nCores)
    }
    else {
      cellsPerBatch <- .SplitCellsIntoBatches(seuratObj, nBatches = nBatches)
      for (i in 1:nBatches) {
        toRun <- cellsPerBatch[[i]]
        print(paste0('Running escape batch ', i, ' of ', nBatches, ' with ', length(toRun), ' cells'))
        so <- subset(seuratObj, cells = toRun)
        if (ncol(so) != length(toRun)) {
          stop(paste0('Error subsetting seurat object, batch size does not match cells after subset: ', length(toRun), ' / ', ncol(seuratObj)))
        }

        mat <- .RunEscapeOnSubset(seuratObj = so, outputAssayName = outputAssayName, GS = GS, nCores = nCores)
        rm(so)

        assayCounts <- cbind(assayCounts, mat)
      }
    }

    if (ncol(assayCounts) != ncol(seuratObj)) {
      stop('The rows of the assay object are not equal to the number of cells')
    }

    assayCounts <- assayCounts[,colnames(seuratObj)]
    if (ncol(assayCounts) != ncol(seuratObj)) {
      stop('The rows of the assay object are not equal to the number of cells, after re-ordering')
    }

    if (any(colnames(assayCounts) != colnames(seuratObj))) {
      stop('The cell names did not match after batch processing')
    }

    seuratObj[[outputAssayName]] <- Seurat::CreateAssayObject(counts = assayCounts)
    seuratObj <- .NormalizeEscape(seuratObj, assayToNormalize = outputAssayName, assayForLibrarySize = assayName)

    if (doPlot) {
      pathways <- rownames(seuratObj@assays[[outputAssayName]])
      key <- seuratObj@assays[[outputAssayName]]@key
      for (fn in pathways) {
        print(suppressWarnings(Seurat::FeaturePlot(seuratObj, features = paste0(key, fn), min.cutoff = 'q02', max.cutoff = 'q98')))
      }
    }

    if (performDimRedux) {
      seuratObj <- .RunEscapePca(seuratObj, assayName = outputAssayName)
    }
  }

  return(seuratObj)
}

.RunEscapeOnSubset <- function(seuratObj, outputAssayName, GS, nCores){
  BPPARAM <- .InferBpParam(nCores, defaultValue = NULL)
  seuratObj <- escape::runEscape(seuratObj,
                                 method = "ssGSEA",
                                 gene.sets = GS,
                                 min.size = 0,
                                 BPPARAM = BPPARAM,
                                 new.assay.name = outputAssayName)

  return(SeuratObject::GetAssayData(seuratObj, assay = outputAssayName, layer = 'data'))
}

.RunEscapePca <- function(seuratObj, assayName, dimsToUse = NULL, resolutionsToUse = 0.2) {
  Seurat::VariableFeatures(seuratObj, assay = assayName) <- rownames(seuratObj@assays[[assayName]])
  seuratObj <- Seurat::ScaleData(seuratObj, assay = assayName)

  assayNameForKeys <- gsub(assayName, pattern = '\\.', replacement = '')
  pca.reduction.key <- paste0(assayNameForKeys, 'pca_')
  pca.reduction.name <- paste0('pca.', assayName)
  seuratObj <- Seurat::RunPCA(seuratObj, assay = assayName, npcs = min(50, length(Seurat::VariableFeatures(seuratObj, assay = assayName))), reduction.key = pca.reduction.key, reduction.name = pca.reduction.name)

  print(Seurat::ProjectDim(seuratObj, reduction = pca.reduction.name, assay = assayName))
  print(Seurat::VizDimLoadings(object = seuratObj, dims = 1:4, nfeatures = nrow(seuratObj@assays[[assayName]]), reduction = pca.reduction.name))

  if (all(is.null(dimsToUse))) {
    npc <- dim(seuratObj@reductions[[pca.reduction.name]]@cell.embeddings)[2]
    dimsToUse <- 1:min(npc, nrow(seuratObj@assays[[assayName]]))
  }

  graphName <- paste0(assayName, '.nn')
  seuratObj <- Seurat::FindNeighbors(seuratObj, dims = dimsToUse, reduction = pca.reduction.name, assay = assayName, graph.name = graphName)

  origIdents <- Seurat::Idents(seuratObj)
  for (resolutionToUse in resolutionsToUse) {
    seuratObj <- Seurat::FindClusters(object = seuratObj, resolution = resolutionToUse, verbose = FALSE, graph.name = graphName, seed.use = GetSeed(), cluster.name = paste0('ClusterNames.', assayName, '_', resolutionToUse))
  }
  Seurat::Idents(seuratObj) <- origIdents

  umap.reduction.name <- paste0(assayName, '.umap')
  umap.reduction.key <- paste0(assayNameForKeys, 'umap_')
  seuratObj <- Seurat::RunUMAP(seuratObj, dims = dimsToUse, assay = assayName, reduction = pca.reduction.name, seed.use = GetSeed(), reduction.name = umap.reduction.name, reduction.key = umap.reduction.key, verbose = FALSE)

  print(DimPlot(seuratObj, reduction = umap.reduction.name))

  return(seuratObj)
}

.NormalizeEscape <- function(seuratObj, assayToNormalize, assayForLibrarySize = 'RNA') {
  toNormalize <- Seurat::GetAssayData(seuratObj, assayToNormalize, layer = 'counts')
  assayForLibrarySizeData <- Seurat::GetAssayData(seuratObj, assay = assayForLibrarySize, layer = 'counts')

  if (any(colnames(toNormalize) != colnames(assayForLibrarySize))) {
    stop(paste0('The assayToNormalize and assayForLibrarySize do not have the same cell names!'))
  }

  margin <- 2
  ncells <- dim(x = toNormalize)[margin]

  for (i in seq_len(length.out = ncells)) {
    x <- toNormalize[, i]
    if (any(is.na(x))) {
      warning('NAs were found in the escape data!')
      x[is.na(x)] <- 0
    }

    librarySize <- sum(assayForLibrarySizeData[, i], na.rm = TRUE)
    if (librarySize == 0) {
      stop(paste0('librarySize was zero for: ', colnames(seuratObj)[i]))
    }

    toNormalize[, i] <- x / librarySize
  }

  seuratObj <- Seurat::SetAssayData(seuratObj, assay = assayToNormalize, layer = 'data', new.data = toNormalize)
  seuratObj <- Seurat::ScaleData(seuratObj, assay = assayToNormalize)

  return(seuratObj)
}

.SplitCellsIntoBatches <- function(seuratObj, nBatches, seed = GetSeed()) {
  if (nBatches == 1) {
    stop('It does not make sense to call this function with a single batch')
  }

  cellsPerBatch <- floor(ncol(seuratObj) / nBatches)
  remainder <- ncol(seuratObj) - (cellsPerBatch * nBatches)

  ret <- list()
  set.seed(seed)
  allCells <- colnames(seuratObj)
  ret[[1]] <- sample(allCells, (cellsPerBatch + remainder), replace=FALSE)
  allCells <- setdiff(allCells, ret[[1]])
  for (i in 2:nBatches) {
    ret[[i]] <- sample(allCells, cellsPerBatch, replace=FALSE)
    allCells <- setdiff(allCells, ret[[i]])
  }

  return(ret)
}