#' @include Utils.R
#' @include Preprocessing.R
#' @import Seurat

utils::globalVariables(
  names = c('nCount_RNA', 'nFeature_RNA', 'p.mito', 'x', 'y', 'p_val_adj', 'avg_logFC', 'groupField', 'cluster', 'pct.1', 'pct.2', 'S.Score_UCell', 'G2M.Score_UCell', 'Phase', 'S.Score_UCell', 'G2M.Score_UCell', 'dims', 'stdev', 'group1', 'group2'),
  package = 'CellMembrane',
  add = TRUE
)


#' @title Read and Filter 10X files.
#'
#' @description Reads in 10X files using Read10X and filters abberent cells using PerformEmptyDropletFiltering and returns a Seurat object.
#' @param dataDir Either the directory holding raw count data (generally the raw_feature_bc_matrix), or the parent 'outs' dir from cellranger
#' @param datasetId This will be used as a prefix for barcodes, and stored in metadata. Also used as the project name for the Seurat object.
#' @param datasetName An optional print-friendly name that will be stored in metadata
#' @param emptyDropNIters The number of iterations to use with PerformEmptyDrops()
#' @param emptyDropsFdrThreshold The FDR threshold to call cells in emptyDrops()
#' @param emptyDropsLower Passed directly to emptyDrops(). The lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets.
#' @param storeGeneIds If true, a map to translate geneId and name (by default rownames will use gene name)
#' @param useEmptyDropsCellRanger If TRUE, will use DropletUtils emptyDropsCellRanger instead of emptyDrops
#' @param nExpectedCells Only applied if emptyDropsCellRanger is selected. Passed to n.expected.cells argument
#' @param previouslyFilteredMatrix An optional filepath to a pre-filtered count matrix in h5 format. If non-null, this file will be read instead of dataDir. Empty drops and/or soupX will be skipped.
#' @param useSoupX If true, SoupX will be run against the run input data, instead of emptyDrops
#' @return A Seurat object.
#' @importFrom magrittr %>%
#' @export
#' @importFrom Seurat Read10X
ReadAndFilter10xData <- function(dataDir, datasetId, datasetName = NULL, emptyDropNIters=10000, emptyDropsFdrThreshold = 0.001, storeGeneIds=TRUE, emptyDropsLower = 100, useEmptyDropsCellRanger = FALSE, nExpectedCells = 8000, useSoupX = FALSE, previouslyFilteredMatrix = NULL) {
  if (!file.exists(dataDir)){
    stop(paste0("File does not exist: ", dataDir))
  }

  if (!dir.exists(dataDir)){
    stop(paste0("File is not a directory: ", dataDir))
  }

  if (!endsWith(dataDir, '/')) {
    dataDir <- paste0(dataDir, '/')
  }

  dirWithFeatureMatrix <- .InferMatrixDir(dataDir)
  if (!is.null(previouslyFilteredMatrix)) {
    print('Using previously filtered count matrix')

    seuratRawData <- Seurat::Read10X_h5(filename = previouslyFilteredMatrix, use.names = TRUE)
    # Drop suffixes:
    colnames(seuratRawData) <- sapply(colnames(seuratRawData), function(x){
      return(unlist(strsplit(x, split = '-'))[1])
    })

  } else if (useSoupX) {
    print('Running SoupX')

    # This expects the outs dir:
    if (dirWithFeatureMatrix == dataDir) {
      stop(paste0('When using SoupX, provide the top-level cellranger outs dir, which contains: raw_feature_bc_matrix'))
    }

    seuratRawData <- .RunSoupX(dataDir)
  } else {
    print('Loading counts and running emptyDrops')

    seuratRawData <- Read10X(data.dir = dirWithFeatureMatrix, strip.suffix = TRUE)
    seuratRawData <- PerformEmptyDropletFiltering(seuratRawData, fdrThreshold = emptyDropsFdrThreshold, emptyDropNIters=emptyDropNIters, emptyDropsLower=emptyDropsLower, useEmptyDropsCellRanger = useEmptyDropsCellRanger, nExpectedCells = nExpectedCells)
  }

  #Cannot have underscores in feature names, Seurat will replace with hyphen anyway.  Perform upfront to avoid warning
  if (sum(grepl(x = rownames(seuratRawData), pattern = '_')) > 0) {
    print('Replacing underscores with hyphens in feature names')
    rownames(seuratRawData) <- gsub(x = rownames(seuratRawData), pattern = '_', replacement = '-')
  }

  seuratObj <- CreateSeuratObj(seuratRawData, datasetId = datasetId, datasetName = datasetName)
  .PrintQcPlots(seuratObj)

  if (useSoupX || !is.null(previouslyFilteredMatrix)) {
    print('Storing unaltered raw counts in assay RNA.orig')
    rawData <- Seurat::Read10X(data.dir = dirWithFeatureMatrix, strip.suffix = TRUE)
    colnames(rawData) <- paste0(datasetId, '_', colnames(rawData))
    rawData <- rawData[,colnames(seuratObj)]
    seuratObj[['RNA.orig']] <- Seurat::CreateAssayObject(rawData)
  }

  if (storeGeneIds) {
    #store IDs in assay metadata
    geneIds <- rownames(Read10X(data.dir = dirWithFeatureMatrix, gene.column = 1, strip.suffix = TRUE))
    names(geneIds) <- rownames(seuratObj)
    assayName <- DefaultAssay(seuratObj)
    seuratObj[[assayName]] <- AddMetaData(seuratObj[[assayName]], metadata = geneIds, col.name = 'GeneId')
  }

  return(seuratObj)
}


#' @title Retrieve the gene IDs from a seuratObj created using ReadAndFilter10xData.
#'
#' @param seuratObj The seurat object
#' @param geneNames A vector of gene names to translate
#' @param throwIfGenesNotFound If true and any of the requested gene names are not found, an error will be thrown.  Otherwise, the result will contain NAs
#' @return A named vector of the gene IDs
#' @export
GetGeneIds <- function(seuratObj, geneNames, throwIfGenesNotFound = TRUE) {
  ret <- NULL

  assayData <- Seurat::GetAssay(seuratObj, assay = Seurat::DefaultAssay(seuratObj))
  featureMeta <- slot(assayData, GetAssayMetadataSlotName(assayData))
  if ('GeneId' %in% colnames(featureMeta)) {
    ret <- featureMeta$GeneId
    names(ret) <- rownames(assayData)
    ret <- ret[geneNames]
  }

  if (all(is.null(ret))) {
  	stop('Expected gene IDs to be stored under GetAssay(seuratObj)@meta.data')
  }

  if (throwIfGenesNotFound & sum(is.na(ret)) > 0) {
    notFound <- paste0(geneNames[is.na(ret)], collapse = ',')
    stop(paste0('Gene names not found: ', notFound))
  } else {
    names(ret) <- geneNames
  }

  return(ret)
}


#' @title Merge Seurat Objects
#' @description Merges a list of Seurat objects
#' @param seuratObjs A named list of seurat objects, optionally named (in which case these will be used as dataset names).
#' @param projectName The project name when creating the final seurat object
#' @param merge.data Passed directly to Seurat::merge
#' @param expectedDefaultAssay If not null, the DefaultAssay on the resulting seurat object will be set to this
#' @param excludedAssays An optional list of assay names to drop prior to merge.
#' @param enforceUniqueCells If true, all inputs must have unique cellbarcodes.
#' @param errorOnBarcodeSuffix In certain cases, software appends a digit (i.e. -1) to the end of cellbarcodes. These can be a problem when trying to make string comparisons. If true, the method will error if these are encountered.
#' @param doGC If true, in an attempt to save memory gc() will be run after each seurat object is merged
#' @param duplicateBarcodeMode This dictates the behavior when duplicate cell barcodes are encountered between merged objects. This can be either: 'exclude-all' (all duplicated dropped from all objects), or 'error'.
#' @return A modified Seurat object.
#' @export
#' @importFrom methods slot
MergeSeuratObjs <- function(seuratObjs, projectName, merge.data = FALSE, expectedDefaultAssay = 'RNA', excludedAssays = c('UCellRanks'), enforceUniqueCells = TRUE, errorOnBarcodeSuffix = FALSE, doGC = FALSE, duplicateBarcodeMode = 'exclude-all'){
  nameList <- names(seuratObjs)
  if (is.null(nameList)) {
    stop('Must provide a named list of seurat objects')
  }

  # Ensure barcodes unique
  encounteredBarcodes <- c()
  duplicates <- c()
  for (datasetId in nameList) {
    print(paste0('Adding dataset: ', datasetId))
    seuratObj <- seuratObjs[[datasetId]]
    seuratObj <- .PossiblyAddBarcodePrefix(seuratObj, datasetId = datasetId, datasetName = NULL)
    if (enforceUniqueCells && length(intersect(encounteredBarcodes, colnames(seuratObj))) > 0) {
      if (duplicateBarcodeMode == 'exclude-all') {
        duplicates <- c(duplicates, intersect(encounteredBarcodes, colnames(seuratObj)))
      } else {
        stop(paste0('Duplicate cellbarcodes found for: ', datasetId))
      }
    }

    if (errorOnBarcodeSuffix && any(grepl(colnames(seuratObj), pattern = '-[0-9]{1,2}$'))) {
      stop(paste0('Encountered barcodes with numeric suffixes (i.e. -1) for dataset: ', datasetId))
    }

    encounteredBarcodes <- c(encounteredBarcodes, colnames(seuratObj))

    if (!all(is.null(excludedAssays))) {
      for (an in excludedAssays) {
        if (an %in% names(seuratObj@assays)) {
          seuratObj[[an]] <- NULL
        }
      }
    }
    seuratObjs[[datasetId]] <- seuratObj
  }

  if (length(duplicates) > 0) {
    print(paste0('Dropping duplicated barcodes, total: ', length(duplicates)))
    for (datasetId in nameList) {
      if (any(duplicates %in% colnames(seuratObjs[[datasetId]]))) {
        print(paste0('Dropping duplicates for: ', datasetId))
        print(paste0('Original cells: ', length(colnames(seuratObjs[[datasetId]]))))
        toKeep <- colnames(seuratObjs[[datasetId]])[!(colnames(seuratObjs[[datasetId]]) %in% duplicates)]
        print(paste0('Retaining: ', length(toKeep)))
        if (length(toKeep) == 0) {
          stop(paste0('There are no cells remaining after dropping duplicates for datasetId: ', datasetId))
        }

        seuratObjs[[datasetId]] <- subset(seuratObjs[[datasetId]], cells = toKeep)

        if (length(toKeep) != length(colnames(seuratObjs[[datasetId]]))) {
          stop('Cell number does not match expected after subset')
        }
      } else {
        print(paste0('No duplicates found in: ', datasetId))
      }
    }
  }

  seuratObj <- .DoMergeSimple(seuratObjs = seuratObjs, projectName = projectName, merge.data = merge.data, expectedDefaultAssay = expectedDefaultAssay, doGC = doGC)

  return(seuratObj)
}


#' @title Run the primary seurat processing steps.
#'
#' @description This is the primary entry point for processing scRNAseq data with Seurat
#' @param seuratObj A Seurat object.
#' @param nVariableFeatures The number of variable features, passed to either FindVariableFeatures or SCTransform
#' @param block.size Passed directly to ScaleData
#' @param variableGenesWhitelist An optional vector of genes that will be included in PCA, beyond the default VariableFeatures()
#' @param variableGenesBlacklist An optional vector of genes that will be excluded from PCA, beyond the default VariableFeatures()
#' @param scaleVariableFeaturesOnly If true, ScaleData will only be performed on VariableFeatures(), which is governed by FindVariableFeatures, variableGenesWhitelist, and variableGenesBlacklist
#' @param featuresToRegress The set of features which will be passed to Seurat::ScaleData vars.to.regress
#' @param includeCellCycleGenesInScaleData If true, cell cycle genes will always be included in the features passed to ScaleData().
#' @param useSCTransform If true, SCTransform will be used in place of the standard Seurat workflow (NormalizeData, ScaleData, FindVariableFeatures)
#' @param additionalFindVariableFeatureArgList A list of arguments passed directly to FindVariableFeatures
#' @param scoreCellCycle If true, ScoreCellCycle will be run to compute Phase, which is stored in meta.data. If a field named Phase already exists, this will be skipped.
#' @param useAlternateG2M If true, this will use a smaller set of G2M genes, defined from: https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv
#' @return A modified Seurat object.
#' @export
NormalizeAndScale <- function(seuratObj, nVariableFeatures = NULL, block.size = 1000, variableGenesWhitelist = NULL, variableGenesBlacklist = NULL, featuresToRegress = c(paste0("nCount_", Seurat::DefaultAssay(seuratObj))), scaleVariableFeaturesOnly = TRUE, includeCellCycleGenesInScaleData = TRUE, useSCTransform = FALSE, additionalFindVariableFeatureArgList = NULL, scoreCellCycle = TRUE, useAlternateG2M = FALSE){
  if (!is.null(featuresToRegress)) {
    if ('p.mito' %in% featuresToRegress) {
      if ('p.mito' %in% names(seuratObj@meta.data)) {
        uniquePMito <- length(unique(seuratObj$p.mito))
      } else {
        uniquePMito <- -1
      }

      if (uniquePMito <= 1) {
        print('Either p.mito not present in the seurat object, or all cells have the same value, skipping from regression')
        featuresToRegress <- featuresToRegress[featuresToRegress != 'p.mito']
      }
    }

    if (length(featuresToRegress) == 0) {
      featuresToRegress <- NULL
    }
  }

  if (useSCTransform) {
    additionalArgs <- list()
    if (!is.null(nVariableFeatures)) {
      print(paste0('SCTransform will use the top ', nVariableFeatures, ' variableFeatures'))
      additionalArgs[['variable.features.n']] <- nVariableFeatures
    }

    seuratObj <- .NormalizeAndScaleSCTransform(seuratObj, featuresToRegress = featuresToRegress, additionalArgs = additionalArgs)
  } else {
    if (is.null(additionalFindVariableFeatureArgList)) {
      additionalFindVariableFeatureArgList <- list()
    }

    if (!is.null(nVariableFeatures)) {
      additionalFindVariableFeatureArgList[['nfeatures']] <- nVariableFeatures
    }

    seuratObj <- .NormalizeAndScaleDefault(seuratObj, featuresToRegress = featuresToRegress, scaleVariableFeaturesOnly = scaleVariableFeaturesOnly, includeCellCycleGenesInScaleData = includeCellCycleGenesInScaleData, block.size = block.size, variableGenesWhitelist = variableGenesWhitelist, variableGenesBlacklist = variableGenesBlacklist, additionalFindVariableFeatureArgList = additionalFindVariableFeatureArgList)
  }

  if (scoreCellCycle) {
    if ('Phase' %in% names(seuratObj@meta.data)) {
      print('Phase column already present, will not re-score cell cycle')
    } else {
      seuratObj <- ScoreCellCycle(seuratObj, useAlternateG2M = useAlternateG2M)
    }

  } else {
    print('Cell cycle scoring will be skipped')
  }

  return(seuratObj)
}

.NormalizeAndScaleSCTransform <- function(seuratObj, featuresToRegress, additionalArgs = NULL, verbose = FALSE) {
	print('Using SCTransform')

	toBind <- additionalArgs
	if (is.null(toBind)) {
		toBind <- list()
	}

	toBind[['vars.to.regress']] <- featuresToRegress
	toBind[['verbose']] <- verbose
	toBind[['return.only.var.genes']] <- FALSE
	print('SCTransform args: ')
	print(toBind)
    toBind[['object']] <- seuratObj

	# To avoid 'reached iteration limit' warnings
	seuratObj <- suppressWarnings(rlang::invoke(SCTransform, toBind))
    seuratObj <- .ClearSeuratCommands(seuratObj)

	return(seuratObj)
}

.NormalizeAndScaleDefault <- function(seuratObj, assayName = 'RNA', additionalFindVariableFeatureArgList = NULL, block.size = 1000, variableGenesWhitelist = NULL, variableGenesBlacklist = NULL, featuresToRegress = c(paste0("nCount_", assayName)), scaleVariableFeaturesOnly = TRUE, includeCellCycleGenesInScaleData = TRUE) {
	seuratObj <- NormalizeData(object = seuratObj, normalization.method = "LogNormalize", verbose = F)

	print('Find variable features:')
	toBind <- additionalFindVariableFeatureArgList
	if (all(is.null(toBind))) {
		toBind <- list()
	}

    if (length(toBind) > 0) {
      print('Additional FindVariableFeatures arguments: ')
      print(toBind)
    }

	toBind[['verbose']] <- FALSE
    toBind[['object']] <- seuratObj

	seuratObj <- rlang::invoke(FindVariableFeatures, toBind)
    seuratObj <- .ClearSeuratCommands(seuratObj)

	if (!all(is.null(variableGenesWhitelist))) {
      variableGenesWhitelist <- RIRA::ExpandGeneList(variableGenesWhitelist)
      missingFeats <- variableGenesWhitelist[!variableGenesWhitelist %in% rownames(seuratObj)]
      if (length(missingFeats) > 0) {
        print(paste0('Not all features in variableGenesWhitelist are present in the seuratObj, missing: ', paste0(missingFeats, collapse = ',')))
        variableGenesWhitelist <- intersect(rownames(seuratObj), variableGenesWhitelist)
      }

      preExisting <- intersect(VariableFeatures(seuratObj), variableGenesWhitelist)
      print(paste0('Adding ', length(variableGenesWhitelist), ' genes to variable gene list, of which ', length(preExisting), ' are already present in VariableFeatures'))
      VariableFeatures(seuratObj) <- unique(c(VariableFeatures(seuratObj), variableGenesWhitelist))
      print(paste0('Total after: ', length(VariableFeatures(seuratObj))))
	}

	if (!all(is.null(variableGenesBlacklist))){
      variableGenesBlacklist <- RIRA::ExpandGeneList(variableGenesBlacklist)
      preExisting <- intersect(VariableFeatures(seuratObj), variableGenesBlacklist)
      print(paste0('Excluding ', length(variableGenesBlacklist), ' gene(s) from the variable gene list, of which ', length(preExisting), ' are present in VariableFeatures'))
      VariableFeatures(seuratObj) <- unique(VariableFeatures(seuratObj)[!(VariableFeatures(seuratObj) %in% variableGenesBlacklist)])
      print(paste0('Total after: ', length(VariableFeatures(seuratObj))))
	}

  .PlotVariableFeatures(seuratObj)

  if (scaleVariableFeaturesOnly) {
		feats <- VariableFeatures(object = seuratObj)
		print(paste0('ScaleData will use the top ', length(feats), ' variableFeatures'))
	} else {
		feats <- rownames(x = seuratObj)
	}

	if (includeCellCycleGenesInScaleData) {
		cc.genes <- .GetCCGenes()
		cc.genes <- cc.genes[which(cc.genes %in% rownames(seuratObj))]
		feats <- unique(c(feats, cc.genes))
	}

  print('Scale data:')
  seuratObj <- ScaleData(object = seuratObj, features = feats, vars.to.regress = featuresToRegress, block.size = block.size, verbose = F)

  return(seuratObj)
}


#' @title Run Seurat PCA
#'
#' @param seuratObj A Seurat object.
#' @param npcs Number of PCs to use for RunPCA()
#' @param variableGeneTable If provided, a table of variable genes will be written to this file
#' @return A modified Seurat object.
#' @export
RunPcaSteps <- function(seuratObj, npcs = 50, variableGeneTable = NULL) {
  if (ncol(seuratObj) < npcs) {
    stop(paste0('More PCs were requested, ', npcs, ', than the number of cells: ', ncol(seuratObj)))
  }

  vg <- VariableFeatures(object = seuratObj)

  print(paste0('Total variable genes: ', length(vg)))
  if (!is.null(variableGeneTable)){
    write.table(sort(vg), file = variableGeneTable, sep = '\t', row.names = F, quote = F, col.names = F)
  }

  seuratObj <- RunPCA(object = seuratObj, features = vg, reduction.name = 'pca', verbose = F, npcs = npcs)
  seuratObj <- ProjectDim(object = seuratObj)

  if (!('SCT' %in% names(seuratObj@assays))) {
  	seuratObj <- JackStraw(object = seuratObj, num.replicate = 100, verbose = F)
  	seuratObj <- ScoreJackStraw(object = seuratObj, dims = 1:20)
  }

  .PrintSeuratPlots(seuratObj)

  return(seuratObj)
}


#' @title Perform Cell Filtering on the cells of the seurat object
#'
#' @description This will perform filtering on the cells of the seurat object
#' @param seuratObj A Seurat object.
#' @param nCount_RNA.high Cells with nCount_RNA greater than this value will be filtered
#' @param nCount_RNA.low Cells with nCount_RNA less than this value will be filtered
#' @param nFeature.high Cells with nFeature greater than this value will be filtered
#' @param nFeature.low Cells with nFeature less than this value will be filtered
#' @param pMito.high Cells with percent mito greater than this value will be filtered
#' @param pMito.low Cells with percent mito  less than this value will be filtered
#' @return The modified seurat object
#' @export
FilterRawCounts <- function(seuratObj, nCount_RNA.high = 20000, nCount_RNA.low = 1, nFeature.high = 3000, nFeature.low = 200, pMito.high = 0.15, pMito.low = 0) {
  print("Filtering Cells...")

  assayName <- Seurat::DefaultAssay(seuratObj)
  nCountField <- paste0('nCount_', assayName)
  nFeatureField <- paste0('nFeature_', assayName)

  print(paste0('Initial cells: ', length(colnames(x = seuratObj))))

  if ('p.mito' %in% colnames(seuratObj@meta.data)) {
	  uniquePMito <- length(unique(seuratObj$p.mito))
    if (uniquePMito > 1) {
      P1 <- FeatureScatter(object = seuratObj, feature1 = nCountField, feature2 = "p.mito")
      P1 <- P1 + geom_vline(aes(xintercept=nCount_RNA.low), color="blue", linetype="dashed", linewidth=1)
      P1 <- P1 + geom_vline(aes(xintercept=nCount_RNA.high), color="blue", linetype="dashed", linewidth=1)
      P1 <- P1 + geom_hline(aes(yintercept=pMito.low), color="blue", linetype="dashed", linewidth=1)
      P1 <- P1 + geom_hline(aes(yintercept=pMito.high), color="blue", linetype="dashed", linewidth=1)
      print(P1)
    } else {
      print("p.mito is either absent or identical across all cells")
    }
  }

	P1 <- FeatureScatter(object = seuratObj, feature1 = nCountField, feature2 = nFeatureField)
	P1 <- P1 + geom_vline(aes(xintercept=nCount_RNA.low), color="blue", linetype="dashed", linewidth=1)
	P1 <- P1 + geom_vline(aes(xintercept=nCount_RNA.high), color="blue", linetype="dashed", linewidth=1)
	P1 <- P1 + geom_hline(aes(yintercept=nFeature.low), color="blue", linetype="dashed", linewidth=1)
	P1 <- P1 + geom_hline(aes(yintercept=nFeature.high), color="blue", linetype="dashed", linewidth=1)
	print(P1)

  #See: https://github.com/satijalab/seurat/issues/1053#issuecomment-454512002
  expr <- Seurat::FetchData(object = seuratObj, vars = nCountField)
  seuratObj <- seuratObj[, which(x = expr >= nCount_RNA.low & expr <= nCount_RNA.high)]
  print(paste0('After ', nCountField, ' filter: ', length(colnames(x = seuratObj))))

  expr <- Seurat::FetchData(object = seuratObj, vars = nFeatureField)
  seuratObj <- seuratObj[, which(x = expr >= nFeature.low & expr <= nFeature.high)]
  print(paste0('After ', nFeatureField, ' filter: ', length(colnames(x = seuratObj))))

  if ('p.mito' %in% colnames(seuratObj@meta.data)) {
    expr <- Seurat::FetchData(object = seuratObj, vars = 'p.mito')
    if (!all(is.na(expr)) && max(expr) != 0) {
      seuratObj <- seuratObj[, which(x = expr >= pMito.low & expr <= pMito.high)]
      print(paste0('After p.mito filter: ', length(colnames(x = seuratObj))))
    } else {
      print('Either p.mito was NA or all values were 0')
    }
  } else {
    print('p.mito was missing, skipping')
  }

  print(paste0('Final cells: ', length(colnames(x = seuratObj))))
  
  return(seuratObj)
}


.PrintSeuratPlots <- function(seuratObj) {
  print(VizDimLoadings(object = seuratObj, dims = 1:2))

	if (!('SCT' %in% names(seuratObj@assays))) {
		suppressWarnings(print(LabelPoints(plot = VariableFeaturePlot(seuratObj), points = head(VariableFeatures(seuratObj), 20), repel = TRUE, xnudge = 0, ynudge = 0)))
	}

  print(DimPlot(object = seuratObj))
  if (('Phase' %in% names(seuratObj@meta.data))) {
		#Note: some cells might lack phase, so suppressWarnings
		suppressWarnings(print(
      DimPlot(object = seuratObj, reduction = "pca", dims = c(1, 2), group.by = 'Phase') +
      DimPlot(object = seuratObj, reduction = "pca", dims = c(2, 3), group.by = 'Phase') +
      DimPlot(object = seuratObj, reduction = "pca", dims = c(3, 4), group.by = 'Phase') +
      DimPlot(object = seuratObj, reduction = "pca", dims = c(4, 5), group.by = 'Phase') +
      patchwork::plot_layout(ncol = 2)
    ))
  }

  totalCells <- min(500, ncol(seuratObj))
  print(DimHeatmap(object = seuratObj, dims = 1, cells = totalCells, balanced = TRUE, fast = FALSE) + NoLegend())
  
  tryCatch({
    print(DimHeatmap(object = seuratObj, dims = 1:20, cells = totalCells, balanced = TRUE, fast = FALSE) + NoLegend())
  }, error = function(){
    try(print(DimHeatmap(object = seuratObj, dims = 1:6, cells = totalCells, balanced = TRUE, fast = FALSE) + NoLegend()), silent = T)
  })

  if (length(seuratObj@reductions$pca@jackstraw$empirical.p.values) == 0) {
    print('Unable to display JackStrawPlot, data not available')
  } else {
		# To avoid: "Removed 30927 rows containing missing values (geom_point)" warning
		suppressWarnings(print(JackStrawPlot(object = seuratObj, dims = 1:20)))
  }

  print(ElbowPlot(object = seuratObj))
}


#' @title Score Cell Cycle
#' @param seuratObj The seurat object
#' @param min.genes If less than min.genes are shared between the seurat object and the reference cell cycle genes, this method will abort.
#' @param useAlternateG2M If true, this will use a smaller set of G2M genes, defined from: https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv
#' @export
#' @return A modified Seurat object.
ScoreCellCycle <- function(seuratObj, min.genes = 10, useAlternateG2M = FALSE) {
	print('Scoring cell cycle:')

  # We can segregate this list into markers of G2/M phase and markers of S-phase
  s.genes <- GetSPhaseGenes()
  g2m.genes <- GetG2MGenes(useAlternateG2M)

  s.genes <- s.genes[which(s.genes %in% rownames(seuratObj))]
  g2m.genes <- g2m.genes[which(g2m.genes %in% rownames(seuratObj))]

  print(paste0("Genes present in seurat object: g2m (", length(g2m.genes), ") and s (", length(s.genes), ")"))

  if (length(g2m.genes) < min.genes || length(s.genes) < min.genes) {
    print(paste0("Error, the number of g2m and/or s genes < ", min.genes, ", aborting..."))
    return(seuratObj)
  }

  print("Running PCA with cell cycle genes")
  seuratObj <- suppressWarnings(RunPCA(object = seuratObj, reduction.name = 'cc.pca', features = c(s.genes, g2m.genes), do.print = FALSE, verbose = F))
  print(DimPlot(object = seuratObj, reduction = "cc.pca"))

  seuratObj <- CellCycleScoring(object = seuratObj,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = FALSE
  )

  print(
    DimPlot(object = seuratObj, group.by = 'Phase', reduction = "cc.pca", dims = c(1, 2)) +
    DimPlot(object = seuratObj, group.by = 'Phase', reduction = "cc.pca", dims = c(2, 3)) +
    DimPlot(object = seuratObj, group.by = 'Phase', reduction = "cc.pca", dims = c(3, 4)) +
    DimPlot(object = seuratObj, group.by = 'Phase', reduction = "cc.pca", dims = c(4, 5)) +
    patchwork::plot_layout(ncol = 2)
  )

  # Discard unnecessary data:
  seuratObj@reductions[['cc.pca']] <- NULL

  print(table(seuratObj$Phase))

  return(seuratObj)
}

#' @title Regress Cell Cycle
#' @param seuratObj The seurat object
#' @param scaleVariableFeaturesOnly If true, ScaleData will only be performed on genes specified in FindVariableFeatures()
#' @param block.size Passed directly to ScaleData
#' @export
#' @return A modified Seurat object.
RegressCellCycle <- function(seuratObj, scaleVariableFeaturesOnly = T, block.size = 1000) {
  print("Regressing out S and G2M score ...")
	if (scaleVariableFeaturesOnly) {
		feats <- VariableFeatures(object = seuratObj)
		if (length(feats) == 0) {
			stop('Must run FindVariableFeatures() upstream on this seurat object when using scaleVariableFeaturesOnly==TRUE')
		}

		print(paste0('ScaleData will use only the ', length(feats), ' variableFeatures'))
	} else {
		feats <- rownames(x = seuratObj)
	}

	usedSCTransform <- 'SCT' %in% names(seuratObj@assays)
	if (usedSCTransform) {
		seuratObj <- SCTransform(seuratObj, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE, return.only.var.genes = FALSE)
	} else {
  	    seuratObj <- ScaleData(object = seuratObj, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE, features = feats, do.scale = T, do.center = T, block.size = block.size)
	}

  return(seuratObj)
}


#' @title Find Clusters And Perform DimRedux
#' @param seuratObj A Seurat object.
#' @param dimsToUse The number of dims to use.  If null, this will be inferred using FindSeuratElbow()
#' @param minDimsToUse The minimum numbers of dims to use.  If dimsToUse is provided, this will override.
#' @param umap.method The UMAP method, either uwot or umap-learn, passed directly to Seurat::RunUMAP
#' @param umap.metric Passed directly to Seurat::RunUMAP
#' @param umap.n.neighbors Passed directly to Seurat::RunUMAP
#' @param umap.min.dist Passed directly to Seurat::RunUMAP
#' @param umap.spread Passed directly to Seurat::RunUMAP
#' @param seed.use Passed directly to Seurat::RunUMAP, FindClusters, and RunTSNE
#' @param umap.n.epochs Passed directly to Seurat::RunUMAP
#' @param umap.densmap Passed directly to Seurat::RunUMAP
#' @param max.tsne.iter The value of max_iter to provide to RunTSNE.  Increasing can help large datasets.
#' @param tsne.perplexity tSNE perplexity. Passed directly to Seurat::RunTSNE(), but CellMembrane:::.InferPerplexityFromSeuratObj() corrects it if need be based dataset dims.
#' @param clusterResolutions A vector of clustering resolutions, default is (0.2, 0.4, 0.6, 0.8, 1.2).
#' @param runTSNE If true, tSNE will be run. The default is UMAP alone.
#' @param useLeiden If true, Leiden clustering (FindClusters algorithm = 4) will be used. Otherwise it will default to algorithm = 1 (Louvain).
#' @return A modified Seurat object.
#' @export
FindClustersAndDimRedux <- function(seuratObj, dimsToUse = NULL, minDimsToUse = NULL,
                                   umap.method = 'uwot', umap.metric = NULL,
                                   umap.n.neighbors = NULL, umap.min.dist = NULL, umap.spread = NULL, seed.use = GetSeed(),
                                   umap.n.epochs = NULL, max.tsne.iter = 10000, tsne.perplexity = 30, umap.densmap = FALSE,
  clusterResolutions = c(0.2, 0.4, 0.6, 0.8, 1.2), runTSNE = FALSE,
                         useLeiden = FALSE){

  dimsToUse <- .GetDimsToUse(seuratObj, dimsToUse = dimsToUse, minDimsToUse = minDimsToUse)

  seuratObj <- FindNeighbors(object = seuratObj, dims = dimsToUse, verbose = FALSE)

  clusterMethod <- ifelse(ncol(seuratObj) > 5000, yes = 'igraph', no = 'matrix')
  algorithm <- ifelse(useLeiden, yes = 4, no = 1)

  for (resolution in clusterResolutions){
    seuratObj <- FindClusters(object = seuratObj, resolution = resolution, verbose = FALSE, random.seed = seed.use, method = clusterMethod, algorithm = algorithm, cluster.name = paste0("ClusterNames_", resolution))
  }

  if (runTSNE) {
    perplexity <- .InferPerplexityFromSeuratObj(seuratObj, perplexity = tsne.perplexity)
    seuratObj <- RunTSNE(object = seuratObj,
                                      dims.use = dimsToUse,
                                      seed.use = seed.use,
                                      check_duplicates = FALSE,
                                      perplexity = perplexity,
                                      reduction.key = 'rnaTSNE_',
                                      max_iter = max.tsne.iter)
  }

  seuratObj <- .RunUMAP(seuratObj,
        dimsToUse = dimsToUse,
        minDimsToUse = minDimsToUse,
        seed.use = seed.use,
        umap.method = umap.method,
        umap.metric = umap.metric,
        umap.n.neighbors = umap.n.neighbors,
        umap.min.dist = umap.min.dist,
        umap.spread = umap.spread,
        umap.n.epochs = umap.n.epochs,
        umap.densmap = umap.densmap
  )

  seuratObj <- .ClearSeuratCommands(seuratObj)

  toPlot <- ifelse(runTSNE, yes = c('tsne', 'umap'), no = 'umap')
  for (reduction in toPlot){
    plotLS <- list()
    i <- 0
    for (res in as.character(clusterResolutions)){
      i <- i + 1
      plotLS[[i]] <- suppressWarnings(DimPlot(object = seuratObj, reduction = reduction, group.by = paste0("ClusterNames_", res), label = TRUE)) + patchwork::plot_annotation(title = paste0("Resolution: ", res))

    }
    print(patchwork::wrap_plots(plotLS))
  }

  return(seuratObj)
}

.GetDimsToUse <- function(seuratObj, dimsToUse, minDimsToUse){
  if (is.null(dimsToUse)) {
    dimMax <- .FindSeuratElbow(seuratObj)
    print(paste0('Inferred elbow: ', dimMax))

    if (!is.null(minDimsToUse)) {
      print(paste0('Min dims to use: ', minDimsToUse))
      dimMax <- max(minDimsToUse, dimMax)
    }

    print(paste0('Selected dimsToUse: 1:', dimMax))
    dimsToUse <- 1:dimMax
  }

  return(dimsToUse)
}

.RunUMAP <- function(seuratObj, reduction = 'pca', dimsToUse = NULL, minDimsToUse = NULL,
                     umap.method = 'uwot', umap.metric = NULL,
                     umap.n.neighbors = NULL, umap.min.dist = NULL,
                     umap.spread = NULL, seed.use = GetSeed(),
                     umap.n.epochs = NULL, umap.densmap = FALSE, reduction.key = 'rnaUMAP_', reduction.name = 'umap') {

  dimsToUse <- .GetDimsToUse(seuratObj, dimsToUse = dimsToUse, minDimsToUse = minDimsToUse)

  umapArgs <- list()
  umapArgs[['verbose']] <- FALSE
  umapArgs[['object']] <- seuratObj
  umapArgs[['dims']] <- dimsToUse
  umapArgs[['reduction.key']] <- reduction.key
  umapArgs[['reduction.name']] <- reduction.name
  umapArgs[['reduction']] <- reduction

  possibleArgs <- list(
    n.neighbors = umap.n.neighbors,
    min.dist = umap.min.dist,
    metric = umap.metric,
    umap.method = umap.method,
    seed.use = seed.use,
    spread = umap.spread,
    densmap = umap.densmap,
    n.epochs = umap.n.epochs
  )

  for (param in names(possibleArgs)) {
    if (!is.null(possibleArgs[[param]])) {
      print(paste0('Setting RunUMAP arg: ', param, ' to ', possibleArgs[[param]]))
      umapArgs[[param]] <- possibleArgs[[param]]
    }
  }

  seuratObj <- rlang::invoke(RunUMAP, umapArgs)

  return(seuratObj)
}


#' @title Find_Markers
#' @param seuratObj A seurat object
#' @param identFields A vector of grouping fields for DE. Often these are the resolution, computed during FindClustersAndDimRedux()
#' @param outFile A file where a table of markers will be saved
#' @param testsToUse A vector of tests to be used.  Each will be used to run FindAllMarkers() and the results merged. Available are: wilcox, bimod, roc, t, negbinom, poisson, LR, MAST, DESeq2
#' @param numGenesToPrint The number of top markers per cluster to print in a table
#' @param onlyPos If true, only positive markers will be saved
#' @param pValThreshold Only genes with adjusted p-values below this will be reported
#' @param foldChangeThreshold Only genes with log2 fold-change above this will be reported
#' @param minPct Only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed.
#' @param minDiffPct Only test genes that show a minimum difference in the fraction of detection between the two groups.
#' @param datasetName An optional label for this dataset. If provided, this will be appended to the resulting table.
#' @param assayName The assay to use.
#' @param verbose Passed to Seurat::FindMarkers
#' @param doPairwise If true, rather than use Seurat::FindAllMarkers, the code will iterate each pair of values and generate a list of DEGs from these comparisons. This can help identify markers highly differentially expressed between some groups, but shared between others.
#' @return A DT::datatable object with the top markers, suitable for printing
#' @importFrom dplyr %>% coalesce group_by summarise filter top_n select everything
#' @import DESeq2
#' @import MAST
#' @export
Find_Markers <- function(seuratObj, identFields, outFile = NULL, testsToUse = c('wilcox', 'MAST', 'DESeq2'), numGenesToPrint = 20, onlyPos = F, pValThreshold = 0.001, foldChangeThreshold = 0.5, minPct = 0.1, minDiffPct = -Inf, datasetName = NULL, assayName = 'RNA', verbose = FALSE, doPairwise = FALSE) {
  if (is.null(assayName)) {
    assayName <- Seurat::DefaultAssay(seuratObj)
    print(paste0('Using default assay: ', assayName))
  }

  seuratObj.markers <- NULL
  fieldsToUse <- c()
  for (fieldName in identFields) {
    # Allow resolution to be passed directly:
    if (!(fieldName %in% names(seuratObj@meta.data))) {
      toTest <- paste0('ClusterNames_', fieldName)
      if (toTest %in% names(seuratObj@meta.data)) {
        fieldName <- toTest
      }
    }

    if (!(fieldName %in% names(seuratObj@meta.data))) {
      stop(paste0('Unknown field: ', fieldName))
    }

    fieldsToUse <- c(fieldsToUse, fieldName)
    print(paste0('Grouping by field: ', fieldName))
    Idents(seuratObj) <- fieldName

    for (test in testsToUse) {
      print(paste0('Running using test: ', test))
      if (doPairwise) {
        vals <- sort(unique(seuratObj@meta.data[[fieldName]]))
        for (val1 in vals) {
          for (val2 in vals) {
            if (val1 == val2) {
              next
            }

            dat <- Seurat::FindMarkers(seuratObj, assay = assayName, ident.1 = val1, ident.2 = val2, group.by = fieldName, test.use = test, min.diff.pct = minDiffPct, min.pct = minPct, random.seed = GetSeed(), verbose = verbose)
            if (nrow(dat) == 0) {
              print(paste0('No passing genes, skipping: ', fieldName, ', ', val1, ' vs. ', val2))
              next
            }

            dat$test <- test
            dat$groupField <- fieldName
            dat$gene <- rownames(dat)
            dat$group1 <- val1
            dat$group2 <- val2

            logFcField <- ifelse('avg_log2FC' %in% colnames(dat), yes = 'avg_log2FC', no = 'avg_logFC')
            if (test == 'roc') {
              toBind <- data.frame(
                groupField = dat$groupField,
                test = as.character(dat$test),
                group1 = as.character(dat$group1),
                group2 = as.character(dat$group2),
                gene = as.character(dat$gene),
                pct.1 = dat$pct.1,
                pct.2 = dat$pct.2,
                avg_logFC = NA,
                p_val_adj = NA,
                myAUC = dat$myAUC,
                power = dat$power,
                avg_diff = dat$avg_diff, stringsAsFactors=FALSE
              )
            } else {
              toBind <- data.frame(
                groupField = dat$groupField,
                test = as.character(dat$test),
                group1 = as.character(dat$group1),
                group2 = as.character(dat$group2),
                gene = as.character(dat$gene),
                pct.1 = dat$pct.1,
                pct.2 = dat$pct.2,
                avg_logFC = dat[[logFcField]],
                p_val_adj = dat$p_val_adj,
                myAUC = NA,
                power = NA,
                avg_diff = NA, stringsAsFactors=FALSE
              )
            }

            seuratObj.markers <- rbind(seuratObj.markers, toBind)
          }
        }
      } else {
        tryCatch({
          tMarkers <- FindAllMarkers(object = seuratObj, assay = assayName, only.pos = onlyPos, logfc.threshold = foldChangeThreshold, min.pct = minPct, min.diff.pct = minDiffPct, verbose = verbose, test.use = test, random.seed = GetSeed())
          if (nrow(tMarkers) == 0) {
            print(paste0('No genes returned, skipping: ', test))
          } else {
            tMarkers$test <- test
            tMarkers$groupField <- fieldName
            tMarkers$cluster <- as.character(tMarkers$cluster)

            logFcField <- ifelse('avg_log2FC' %in% colnames(tMarkers), yes = 'avg_log2FC', no = 'avg_logFC')
            if (test == 'roc') {
              toBind <- data.frame(
                groupField = tMarkers$groupField,
                test = as.character(tMarkers$test),
                cluster = as.character(tMarkers$cluster),
                gene = as.character(tMarkers$gene),
                pct.1 = tMarkers$pct.1,
                pct.2 = tMarkers$pct.2,
                avg_logFC = NA,
                p_val_adj = NA,
                myAUC = tMarkers$myAUC,
                power = tMarkers$power,
                avg_diff = tMarkers$avg_diff, stringsAsFactors=FALSE
              )
            } else {
              toBind <- data.frame(
                groupField = tMarkers$groupField,
                test = as.character(tMarkers$test),
                cluster = as.character(tMarkers$cluster),
                gene = as.character(tMarkers$gene),
                pct.1 = tMarkers$pct.1,
                pct.2 = tMarkers$pct.2,
                avg_logFC = tMarkers[[logFcField]],
                p_val_adj = tMarkers$p_val_adj,
                myAUC = NA,
                power = NA,
                avg_diff = NA, stringsAsFactors=FALSE
              )
            }

            print(paste0('Total genes returned: ', nrow(toBind)))

            if (all(is.null(seuratObj.markers))) {
              seuratObj.markers <- toBind
            } else {
              seuratObj.markers <- rbind(seuratObj.markers, toBind)
            }
          }
        }, error = function(e){
          print(paste0('Error running test: ', test))
          print(conditionMessage(e))
          traceback()
          print(utils::str(tMarkers))
          print(utils::str(seuratObj.markers))
        })
      }
    }
  }

  if (all(is.null(seuratObj.markers))) {
    print('All tests failed, no markers returned')
    return()
  }
  else if (nrow(seuratObj.markers) == 0) {
    print('No significant markers were found')
    return()
  } else {
    if (!doPairwise) {
      seuratObj.markers$cluster <- naturalsort::naturalfactor(seuratObj.markers$cluster)
    } else {
      seuratObj.markers$group1 <- naturalsort::naturalfactor(seuratObj.markers$group1)
      seuratObj.markers$group2 <- naturalsort::naturalfactor(seuratObj.markers$group2)
    }

    toWrite <- seuratObj.markers %>% filter(p_val_adj < pValThreshold) %>% filter(avg_logFC > foldChangeThreshold)
    if (nrow(toWrite) == 0) {
      print('No significant markers were found')
    } else {
      print(paste0('Total DE genes: ', length(unique(toWrite$gene))))
      if (!is.null(outFile)) {
        write.table(toWrite, file = outFile, sep = '\t', row.names = F, quote = F)
      }

      print(DimPlot(object = seuratObj))

      for (fieldName in fieldsToUse) {
        toPlot <- toWrite[toWrite$groupField == fieldName,]
        if (nrow(toPlot) == 0) {
          next
        }

        allGenes <- NULL
        if (!doPairwise) {
          print(ggplot(toPlot, aes(x = pct.1, y = pct.2, size = avg_logFC, color = cluster)) +
                  geom_point(alpha = 0.5) +
                  ggtitle(paste0('DE Genes: ', fieldName)) +
                  egg::theme_presentation(base_size = 12)
          )

          topGene <- toPlot %>% group_by(cluster, test) %>% top_n(numGenesToPrint, avg_logFC)
          if (length(unique(topGene$gene)) < 3) {
            print('Too few genes, skipping heatmap')
            next
          }
          allGenes <- rbind(allGenes, topGene)

          avgSeurat <- Seurat::AverageExpression(seuratObj, group.by = fieldName, features = unique(topGene$gene), layer = 'counts', assays = assayName, return.seurat = T)
          avgSeurat <- NormalizeData(avgSeurat, verbose = FALSE)

          # Genes as columns:
          mat <- t(as.matrix(Seurat::GetAssayData(avgSeurat, layer = 'data')))
          plot(ComplexHeatmap::Heatmap(mat %>% scale_mat(scale = 'column'),
                                       column_title = fieldName,
                                       row_names_side = "left",
                                       row_dend_side = "right",
                                       col = Seurat::BlueAndRed(10),
                                       column_names_side = "top",
                                       column_dend_side = "bottom"
          ))
        } else {
          topGene <- toPlot %>% group_by(group1, group2, test) %>% top_n(numGenesToPrint, avg_logFC)
          if (length(unique(topGene$gene)) < 3) {
            print('Too few genes, skipping heatmap')
            next
          }
          allGenes <- rbind(allGenes, topGene)
        }
      }

      #Note: return the datatable, so it will be printed correctly by Rmarkdown::render()
      return(DT::datatable(topGene,
        caption = paste0('Top DE Genes', ifelse(is.null(datasetName), yes = '', no = paste0(': ', datasetName))),
        filter = 'none',
        escape = FALSE,
        extensions = 'Buttons',
        options = list(
          dom = 'Bfrtip',
          pageLength = 25,
          scrollX = TRUE,
          buttons = c('excel', "csv")
         )
      ))
    }
  }
}

#' @import ggplot2
.FindSeuratElbow <- function(object, ndims = 25, reduction = "pca", print.plot = T, min.y = 1.3) {
  data.use <- Stdev(object = object, reduction = reduction)

  if (length(data.use) == 0) {
    stop(paste("No standard deviation info stored for", reduction))
  }
  if (ndims > length(x = data.use)) {
    warning("The object only has information for ", length(x = data.use),
    " reductions")
    ndims <- length(x = data.use)
  }

  #1 sd = 1.3
  elbowX <- try(.FindElbow(data.use[1:ndims], plot = FALSE, ignore.concavity = F, min.y = min.y))
  if (class(elbowX)=="try-error" || elbowX[1]==2) {
    if (is.null(ndims)){
      elbowX <- 2
    } else {
      elbowX <- ndims
    }
  }

  plot <- ggplot(data = data.frame(dims = 1:ndims, stdev = data.use[1:ndims])) +
    geom_point(mapping = aes(x = dims, y = stdev)) +
    labs(x = gsub(pattern = "_$", replacement = "", x = Key(object = object[[reduction]])),
    y = "Standard Deviation") + theme_bw() + geom_vline(xintercept = elbowX) + ggtitle("Elbow Identification")

  if (print.plot) {
    print(plot)
  }

  return(elbowX)
}


#' @importFrom stats coef
#' @noRd 
.FindElbow <- function(y, plot = FALSE, ignore.concavity = FALSE, min.y = NA, min.x = NA) {

  # minor modification to debug specic scenarios when fail to find elbow
  # The following helper functions were found at
  # paulbourke.net/geometry/pointlineplane/pointline.r
  # via the SO reference below.  The short segment check
  # was modified to permit vectorization.

  ##========================================================
  ##
  ##  Credits:
  ##  Theory by Paul Bourke http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/
  ##  Based in part on C code by Damian Coventry Tuesday, 16 July 2002
  ##  Based on VBA code by Brandon Crosby 9-6-05 (2 dimensions)
  ##  With grateful thanks for answering our needs!
  ##  This is an R (http://www.r-project.org) implementation by Gregoire Thomas 7/11/08
  ##
  ##========================================================
  #' @examples
  #' \dontrun{
  #' tmp <- .FindElbow(c(0.9, 1.1, 1.1, 1.9, 2.5, 2.8, 4.9, 8.5),
  #' 	plot = TRUE) # wandering
  #' tmp <- .FindElbow(c(0.9, 1.0, 1.2, 1.3, 1.5, 1.5, 10.0, 22.0),
  #' 	plot = TRUE) # late rise
  #' tmp <- .FindElbow(c(2, 4, 6, 8, 10, 12, 14, 16)^2,
  #' 	plot = TRUE) # gradual, no obvious break
  #'
  #' # Not the usual way to choose the number of PCs:
  #' library("chemometrics")
  #' data(glass)
  #' pca <- prcomp(glass)
  #' eigensum <- sum(pca$sdev * pca$sdev)
  #' vv <- 100 * (pca$sdev * pca$sdev/eigensum)
  #' cs <- cumsum(vv)
  #' tmp <- .FindElbow(vv, plot = TRUE)
  #' tmp <- .FindElbow(cs, plot = TRUE)
  #' }

  distancePointLine <- function(x, y, slope, intercept) {
    ## x, y is the point to test.
    ## slope, intercept is the line to check distance.
    ##
    ## Returns distance from the line.
    ##
    ## Returns 9999 on 0 denominator conditions.
    x1 <- x-10
    x2 <- x+10
    y1 <- x1*slope+intercept
    y2 <- x2*slope+intercept
    distancePointSegment(x,y, x1,y1, x2,y2)
  }

  distancePointSegment <- function(px, py, x1, y1, x2, y2) {
    ## px,py is the point to test.
    ## x1,y1,x2,y2 is the line to check distance.
    ##
    ## Returns distance from the line, or if the intersecting point on the line nearest
    ## the point tested is outside the endpoints of the line, the distance to the
    ## nearest endpoint.
    ##
    ## Returns 9999 on 0 denominator conditions.
    lineMagnitude <- function(x1, y1, x2, y2) sqrt((x2-x1)^2+(y2-y1)^2)
    lineMag <- lineMagnitude(x1, y1, x2, y2)
    if (any(lineMag < 0.00000001)) {
      #warning("short segment")
      #return(9999)
      warning("At least one line segment given by x1, y1, x2, y2 is very short.")
    }
    u <- (((px - x1) * (x2 - x1)) + ((py - y1) * (y2 - y1)))
    u <- u / (lineMag * lineMag)

    ans <- c()
    ix <- 0
    iy <- 0
    for (i in seq_along(u)) {
      if (u[i] < 0.00001 || u[i] > 1) {
        print('closest point does not fall within the line segment, take the shorter distance to an endpoint')
        ix <- lineMagnitude(px[i], py[i], x1[i], y1[i])
        iy <- lineMagnitude(px[i], py[i], x2[i], y2[i])
        ans <- c(ans, min(ix, iy))
      } else {
        ## Intersecting point is on the line, use the formula
        ix <- x1 + u * (x2 - x1)
        iy <- y1 + u * (y2 - y1)
        ans <- c(ans, lineMagnitude(px, py, ix, iy)[i])
      }
    }
    ans
  }

  # End of helper functions by PB

  ### Now for the actual FindElbow function!

  # Find the elbow using the method described in
  # stackoverflow.com/a/2022348/633251
  # but translated to R (see above).


  y <- sort(y, decreasing = T)

  # Add an index to argument values for easy plotting
  DF <- data.frame(x = seq_along(y), y = y)
  fit <- lm(y ~ x, DF[c(1,nrow(DF)),]) # 2 point 'fit'
  m <- coef(fit)[2]
  b <- coef(fit)[1]

  # Check to make sure the data is concave as described
  # in the documentation, as arbitrary trends could give
  # misleading answers.  The following approach simply
  # checks to make sure all values are either above or
  # below the reference line.  This allows the values
  # to vary quite a bit and still return an answer.

  concave <- FALSE
  use <- 2:(nrow(DF)-1)
  refpts <- m*DF$x[use] + b
  if (all(refpts > DF$y[use]) | all(refpts < DF$y[use])) concave <- TRUE
  if (ignore.concavity) concave <- TRUE

  if (!concave) {
    stop("Your curve doesn't appear to be concave")
  }

  # Calculate the orthogonal distances
  if (is.na(min.x)){
    if (!is.na(min.y)){
      if (!length(which(DF$y<=min.y))<1){
        min.x <- min(DF[which(DF$y<=min.y), ]$x)
      } else {
        print("min.y greater than smallest y")
        min.x <- 2
      }
    } else {
      print("min.x and min.y are NA")
      min.x <- 2
    }
  }

  use     <- min.x:(nrow(DF)-1)
  elbowd  <- distancePointLine(DF$x[use], DF$y[use], coef(fit)[2], coef(fit)[1])
  DF$dist <- rep(NA, nrow(DF))
  DF$dist[use]  <- elbowd # c(NA, elbowd, NA) # first & last points don't have a distance

  if (plot) {
    edm <- which.max(DF$dist)
    plot(DF[,1:2], type = "b", xlab = "index", ylab = "y values",
    main = "Looking for the Elbow")
    segments(DF$x[1], DF$y[1],
    DF$x[nrow(DF)], DF$y[nrow(DF)], col = "red")
    points(DF$x[edm], DF$y[edm], cex = 1.5, col = "red")
    points(DF$x[edm], DF$y[edm], pch = 20)
  }

  if (is.na(which.max(DF$dist))) {
    #if all fails return 2
    print('Max not found, defaulting to 2')
    return(2)
  } else {
    return(which.max(DF$dist))
  }
}

.PlotVariableFeatures <- function(seuratObj) {
  if (length(VariableFeatures(seuratObj)) == 0) {
    stop('There are no VariableFeatures in this seurat object')
  }

  cn <- class(seuratObj[[DefaultAssay(seuratObj)]])[1]
  df <- slot(GetAssay(seuratObj, assay = DefaultAssay(seuratObj)), GetAssayMetadataSlotName(GetAssay(seuratObj, assay = DefaultAssay(seuratObj))))
  if (all(is.null(df))) {
    stop(paste0('Unable to find assay metadata for assay: ', DefaultAssay(seuratObj), ' with class: ', cn, ', expected slot: ', GetAssayMetadataSlotName(GetAssay(seuratObj, assay = DefaultAssay(seuratObj)))))
  }

  varianceStandardizedFields <- c('vf_vst_counts_variance.standardized', 'vst.variance.standardized')
  vstVariableFields <- c('vf_vst_counts_variable', 'vst.variable')

  varianceStandardizedField <- NULL
  vstVariableField <- NULL
  for (idx in seq_along(varianceStandardizedFields)) {
    fn1 <- varianceStandardizedFields[idx]
    fn2 <- vstVariableFields[idx]
    if (fn1 %in% names(df) && fn2 %in% names(df)) {
      varianceStandardizedField <- fn1
      vstVariableField <- fn2

      break
    }
  }

  if (is.null(varianceStandardizedField)) {
    stop(paste0('Unable to find variance fields, expected one of: ', paste0(varianceStandardizedFields, collapse = ','), ' in assay metadata (assay class: ', cn, '). Found: ', paste0(names(df), collapse = ',')))
  }

  dat <- sort(df[[varianceStandardizedField]])
  dat <- dat[dat > 0]
  if (length(dat) == 0) {
    print('No features have vst.variance.standardized > 0')
  }

  countAbove <- sapply(dat, function(x){
    sum(dat >= x)
  })

  print(paste0('Total variable features: ', length(Seurat::VariableFeatures(seuratObj))))
  minVal <- min(df[[varianceStandardizedField]][df[[vstVariableField]]])

  print(ggplot(data.frame(x = countAbove, y = dat), aes(x = x, y = y)) +
          geom_point() +
          ylab("vst.variance.standardized") +
          xlab("# Features") +
          geom_hline(yintercept = minVal, color = 'red') +
          egg::theme_presentation(base_size = 12) +
          scale_x_continuous(trans = 'log10') +
          scale_y_continuous(trans = 'log10') +
          ggtitle('Top Variable Features')
  )

  # This allows a summary of the source of the top variable features, for situations like merging nimble with cellranger:
  if ('FeatureSource' %in% names(df)) {
    x <- df$FeatureSource[rownames(df) %in% Seurat::VariableFeatures(seuratObj)]
    x <- sort(table(x), decreasing = T)
    print(x)
  }
}

#' @title PerformIntegration
#'
#' @description This perform's Seurat's integration on a dataset
#' @param seuratObj A Seurat object.
#' @param splitField The metadata field on which to split the seurat object for integration
#' @param nVariableFeatures The number of variable features used for Seurat::FindVariableFeatures
#' @param nIntegrationFeatures The number of features for Seurat::SelectIntegrationFeatures
#' @param k.weight Passed to Seurat::IntegrateData()
#' @param dimsToUse Passed to Seurat::IntegrateData(), as dims = 1:dimsToUse
#' @param integrationFeaturesInclusionList An optional vector of genes that will be included in the integration genes
#' @param integrationFeaturesExclusionList An optional vector of genes that will be excluded in the integration genes
#' @param minCellsPerSubsetObject If any of the seurat objects after splitting on splitField have fewer than this many cells, they are either discarded or the function will throw an error, depending on dropGroupsBelowThreshold
#' @param dropGroupsBelowThreshold If any of the seurat objects after splitting are below minCellsPerSubsetObject and if this param is true, they will be discarded.
#' @return A modified Seurat object.
#' @export
PerformIntegration <- function(seuratObj, splitField = "SubjectId", nVariableFeatures = 4000, nIntegrationFeatures = 3500, k.weight = 20, dimsToUse = 20, integrationFeaturesInclusionList = NULL, integrationFeaturesExclusionList = NULL, minCellsPerSubsetObject = 75, dropGroupsBelowThreshold = TRUE) {
  if (!splitField %in% names(seuratObj@meta.data)) {
    stop(paste0('splitField not found: ', splitField))
  }

  print(paste0('Splitting on: ', splitField))
  print(sort(table(seuratObj@meta.data[[splitField]])))

  minGroupSize <- min(table(seuratObj@meta.data[[splitField]]))
  if (minGroupSize < minCellsPerSubsetObject) {
    if (dropGroupsBelowThreshold) {
      for (groupVal in unique(seuratObj@meta.data[[splitField]])) {
        cells <- seuratObj@meta.data[[splitField]] == groupVal
        if (sum(cells) < minGroupSize) {
          print(paste0('Dropping low count group: ', groupVal, ', with # cells: ', sum(cells)))
          seuratObj <- subset(seuratObj, cells = colnames(seuratObj)[!cells])
        }
      }
    } else {
      stop(paste0('One or more values of ', splitField, ' has fewer cells than allowed by minCellsPerSubsetObject'))
    }
  }

  Combo_LS <- Seurat::SplitObject(seuratObj, split.by = splitField)

  # normalize and identify variable features for each dataset independently
  Combo_LS <- lapply(X = Combo_LS, FUN = function(x) {
    x <- Seurat::NormalizeData(x, verbose = FALSE)
    x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = nVariableFeatures)
  })

  # select features that are repeatedly variable across datasets for integration run PCA on each
  features <- Seurat::SelectIntegrationFeatures(object.list = Combo_LS, nfeatures = nIntegrationFeatures)

  if (!all(is.null(integrationFeaturesInclusionList))) {
    integrationFeaturesInclusionList <- RIRA::ExpandGeneList(integrationFeaturesInclusionList)
    preExisting <- intersect(features, integrationFeaturesInclusionList)
    print(paste0('Adding ', length(integrationFeaturesInclusionList), ' features to the integration list, of which ', length(preExisting), ' are already present'))
    features <- unique(c(features, integrationFeaturesInclusionList))
    print(paste0('Total after: ', length(features)))
  }

  if (!all(is.null(integrationFeaturesExclusionList))){
    integrationFeaturesExclusionList <- RIRA::ExpandGeneList(integrationFeaturesExclusionList)
    preExisting <- intersect(features, integrationFeaturesExclusionList)
    print(paste0('Excluding ', length(integrationFeaturesExclusionList), ' features(s) from the integration list, of which ', length(preExisting), ' are present'))
    features <- unique(features[!(features %in% integrationFeaturesExclusionList)])
    print(paste0('Total after: ', length(features)))
  }

  print(paste0('Total features for integration: ', length(features)))

  # scale and run the PCA
  Combo_LS <- lapply(X = Combo_LS, FUN = function(x) {
    x <- Seurat::ScaleData(x, features = features, verbose = F)
    x <- Seurat::RunPCA(x, features = features, verbose = F,  npcs = 35)
  })

  seuratAnchors <- Seurat::FindIntegrationAnchors(object.list = Combo_LS, anchor.features = features, reduction = "rpca")

  # this command creates an 'integrated' data assay, dims 1:20 is pretty default but it does make a difference on how many PCA dims to anchor on...
  seuratObj <- Seurat::IntegrateData(anchorset = seuratAnchors, dims = 1:dimsToUse, k.weight = k.weight)
  Seurat::DefaultAssay(seuratObj) <- "Integrated"
  seuratObj <- Seurat::ScaleData(seuratObj, verbose = F)

  # Set this, since the PCA steps use VariableFeatures() to choose the input features
  VariableFeatures(seuratObj) <- features

  return(seuratObj)
}

#' @title CellCycleScoring_UCell
#' @description Similar to Seurat's CellCycleScoring, except using UCell AddModuleScore_UCell instead of Seurat's AddModuleScore
#' @param seuratObj The seurat object
#' @param outputFieldName The name of the meta.data field to store results
#' @param min.genes If less than min.genes are shared between the seurat object and the reference cell cycle genes, this method will abort.
#' @param s.genes A vector of features associated with S phase
#' @param g2m.genes A vector of features associated with G2M phase
#' @param assayName The name of the assay
#' @param min.score.G2M The minimum UCell score needed to consider a cell positive for G2M
#' @param min.score.S The minimum UCell score needed to consider a cell positive for S
#' @param facetField If not NA, a plot will be produced summarizing S.Core and G2M.Score, faceted by this value.
#' @param ncores Passed directly to AddModuleScore_UCell
#' @export
#' @return A modified Seurat object.
CellCycleScoring_UCell <- function(seuratObj, outputFieldName = 'Phase', min.genes = 10, s.genes = GetSPhaseGenes(), g2m.genes = GetG2MGenes(TRUE), assayName = 'RNA', min.score.G2M = 0.02, min.score.S = 0.02, facetField = 'ClusterNames_0.2', ncores = 1) {
  BPPARAM <- .InferBpParam(ncores, defaultValue = NULL)

  s.genes <- s.genes[which(s.genes %in% rownames(seuratObj))]
  g2m.genes <- g2m.genes[which(g2m.genes %in% rownames(seuratObj))]

  print(paste0("Genes present in seurat object: g2m (", length(g2m.genes), ") and s (", length(s.genes), ")"))

  if (length(g2m.genes) < min.genes || length(s.genes) < min.genes) {
    print(paste0("Error, the number of g2m and/or s genes < ", min.genes, ", aborting..."))
    return(seuratObj)
  }

  seuratObj <- UCell::AddModuleScore_UCell(seuratObj, features = list(
    S.Score = s.genes,
    G2M.Score = g2m.genes
  ), assay = assayName, BPPARAM = BPPARAM)

  seuratObj[[outputFieldName]] <- apply(
    X = seuratObj@meta.data,
    MARGIN = 1,
    FUN = function(scores, first = 'S', second = 'G2M', nullLabel = 'G1', threshS = min.score.S, threshG2M = min.score.G2M) {
      scores <- scores[c('S.Score_UCell', 'G2M.Score_UCell')]
      names(scores) <- c('S', 'G2M')

      if (scores[['S']] < threshS && scores[['G2M']] < threshG2M) {
        return(nullLabel)
      } else if (scores[['S']] >= threshS && scores[['G2M']] < threshG2M) {
        return(first)
      } else if (scores[['S']] < threshS && scores[['G2M']] >= threshG2M) {
        return(second)
      } else {
        if (length(which(x = scores == max(scores))) > 1) {
          return('Undecided')
        } else {
          # NOTE: this logic is not quite right, since we probably should not automatically take the highest score
          return(c(first, second)[which(x = scores == max(scores))])
        }
      }
    }
  )

  if (!is.null(facetField)) {
    facets <- stats::as.formula(paste0('. ~ ', facetField))
    P1 <- ggplot(seuratObj@meta.data, aes(x = S.Score_UCell, y = G2M.Score_UCell, color = !!sym(outputFieldName))) +
      geom_point() +
      geom_hline(yintercept = min.score.G2M, color = 'red', linetype = 'dashed') +
      geom_vline(xintercept = min.score.S, color = 'red', linetype = 'dashed') +
      facet_wrap(facets, ncol = 4)

    print(P1)
  }

  print(Seurat::RidgePlot(seuratObj, features = c('S.Score_UCell', 'G2M.Score_UCell'), ncol = 1, group.by = outputFieldName))
  print(Seurat::FeatureScatter(seuratObj, feature1 = 'S.Score_UCell', feature2 = 'G2M.Score_UCell', group.by = outputFieldName) +
          geom_vline(xintercept = min.score.S) +
          geom_hline(yintercept = min.score.G2M))

  return(seuratObj)
}