#' @include Utils.R
#' @include Preprocessing.R
#' @import Seurat

utils::globalVariables(
  names = c('nCount_RNA', 'nFeature_RNA', 'p.mito', 'x', 'y', 'p_val_adj', 'avg_logFC', 'groupField', 'cluster'),
  package = 'CellMembrane',
  add = TRUE
)


#' @title Read and Filter 10X files.
#'
#' @description Reads in 10X files using Read10X and filters abberent cells using PerformEmptyDropletFiltering and returns a Seurat object.
#' @param dataDir The directory holding raw count data, generally the raw_feature_bc_matrix from the cellranger outs folder
#' @param datasetId This will be used as a prefix for barcodes, and stored in metadata. Also used as the project name for the Seurat object.
#' @param datasetName An optional print-friendly name that will be stored in metadata
#' @param emptyDropNIters The number of iterations to use with PerformEmptyDrops()
#' @param emptyDropsLower Passed directly to emptyDrops(). The lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets.
#' @param storeGeneIds If true, a map to translate geneId and name (by default rownames will use gene name)
#' @param annotateMitoFromReference If true, a list of mitochondrial genes, taken from (https://www.genedx.com/wp-content/uploads/crm_docs/Mito-Gene-List.pdf) will be used to calculate p.mito
#' @return A Seurat object.
#' @export
#' @importFrom Seurat Read10X
ReadAndFilter10xData <- function(dataDir, datasetId, datasetName = NULL, emptyDropNIters=10000, storeGeneIds=TRUE, emptyDropsLower = 100, annotateMitoFromReference = TRUE) {
  if (!file.exists(dataDir)){
    stop(paste0("File does not exist: ", dataDir))
  }

  if (!dir.exists(dataDir)){
    stop(paste0("File is not a directory: ", dataDir))
  }

  seuratRawData <- Read10X(data.dir = dataDir, strip.suffix = TRUE)

  #Cannot have underscores in feature names, Seurat will replace with hyphen anyway.  Perform upfront to avoid warning
  if (sum(grepl(x = rownames(seuratRawData), pattern = '_')) > 0) {
    print('Replacing underscores with hyphens in feature names')
    rownames(seuratRawData) <- gsub(x = rownames(seuratRawData), pattern = '_', replacement = '-')
  }

  seuratRawData <- PerformEmptyDropletFiltering(seuratRawData, emptyDropNIters=emptyDropNIters, emptyDropsLower=emptyDropsLower)

  seuratObj <- CreateSeuratObj(seuratRawData, datasetId = datasetId, datasetName = datasetName, annotateMitoFromReference = TRUE)
  .PrintQcPlots(seuratObj)

  if (storeGeneIds) {
    #store IDs in assay metadata
    geneIds <- rownames(Read10X(data.dir = dataDir, gene.column = 1, strip.suffix = TRUE))
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

  featureMeta <- GetAssay(seuratObj)@meta.features
  if ('GeneId' %in% colnames(featureMeta)) {
    ret <- featureMeta$GeneId
    names(ret) <- rownames(seuratObj)
    ret <- ret[geneNames]
  }

  if (all(is.null(ret))) {
  	stop('Expected gene IDs to be stored under GetAssay(seuratObj)@meta.features')
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
#' @return A modified Seurat object.
#' @export
#' @importFrom methods slot
MergeSeuratObjs <- function(seuratObjs, projectName, merge.data = FALSE){
  nameList <- names(seuratObjs)
  if (is.null(nameList)) {
    stop('Must provide a named list of seurat objects')
  }

  # Ensure barcodes unique
  for (datasetId in nameList) {
    print(paste0('Adding dataset: ', datasetId))
    seuratObj <- seuratObjs[[datasetId]]
		seuratObjs[[datasetId]] <- .PossiblyAddBarcodePrefix(seuratObj, datasetId = datasetId, datasetName = NULL)
  }
  
  seuratObj <- .DoMergeSimple(seuratObjs = seuratObjs, projectName = projectName, merge.data = merge.data)

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
#' @return A modified Seurat object.
#' @export
NormalizeAndScale <- function(seuratObj, nVariableFeatures = NULL, block.size = 1000, variableGenesWhitelist = NULL, variableGenesBlacklist = NULL, featuresToRegress = c("nCount_RNA"), scaleVariableFeaturesOnly = TRUE, includeCellCycleGenesInScaleData = TRUE, useSCTransform = FALSE, additionalFindVariableFeatureArgList = NULL){
	if (!is.null(featuresToRegress)) {
		if ('p.mito' %in% featuresToRegress) {
			if ('p.mito' %in% names(seuratObj@meta.data)) {
				uniquePMito = length(unique(seuratObj$p.mito))
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

	seuratObj <- ScoreCellCycle(seuratObj)

	return(seuratObj)
}

.NormalizeAndScaleSCTransform <- function(seuratObj, featuresToRegress, additionalArgs = NULL, verbose = TRUE) {
	print('Using SCTransform')

	toBind <- additionalArgs
	if (is.null(toBind)) {
		toBind <- list()
	}

	toBind[['object']] <- seuratObj
	toBind[['vars.to.regress']] <- featuresToRegress
	toBind[['verbose']] <- verbose
	toBind[['return.only.var.genes']] <- FALSE

	# To avoid 'reached iteration limit' warnings
	seuratObj <- suppressWarnings(do.call(SCTransform, toBind))

	return(seuratObj)
}

.NormalizeAndScaleDefault <- function(seuratObj, additionalFindVariableFeatureArgList = NULL, block.size = 1000, variableGenesWhitelist = NULL, variableGenesBlacklist = NULL, featuresToRegress = c("nCount_RNA"), scaleVariableFeaturesOnly = TRUE, includeCellCycleGenesInScaleData = TRUE) {
	seuratObj <- NormalizeData(object = seuratObj, normalization.method = "LogNormalize", verbose = F)

	print('Find variable features:')
	toBind <- additionalFindVariableFeatureArgList
	if (is.null(toBind)) {
		toBind <- list()
	}

	toBind[['object']] <- seuratObj
	toBind[['verbose']] <- FALSE

	seuratObj <- do.call(FindVariableFeatures, toBind)

	if (!all(is.null(variableGenesWhitelist))) {
		print(paste0('Adding ', length(variableGenesWhitelist), ' genes to variable gene list'))
		VariableFeatures(seuratObj) <- unique(c(VariableFeatures(seuratObj), variableGenesWhitelist))
		print(paste0('Total after', length(VariableFeatures(seuratObj))))
	}

	if (!all(is.null(variableGenesBlacklist))){
		print(paste0('Removing ', length(variableGenesBlacklist), ' from variable gene list'))
		VariableFeatures(seuratObj) <- unique(VariableFeatures(seuratObj)[!(VariableFeatures(seuratObj) %in% variableGenesBlacklist)])
		print(paste0('Total after', length(VariableFeatures(seuratObj))))
	}

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

  print(paste0('Initial cells: ', length(colnames(x = seuratObj))))

  if ('p.mito' %in% colnames(seuratObj@meta.data)) {
	  uniquePMito = length(unique(seuratObj$p.mito))
    if (uniquePMito > 1) {
      P1 <- FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", feature2 = "p.mito")
      P1 <- P1 + geom_vline(aes(xintercept=nCount_RNA.low), color="blue", linetype="dashed", size=1)
      P1 <- P1 + geom_vline(aes(xintercept=nCount_RNA.high), color="blue", linetype="dashed", size=1)
      P1 <- P1 + geom_hline(aes(yintercept=pMito.low), color="blue", linetype="dashed", size=1)
      P1 <- P1 + geom_hline(aes(yintercept=pMito.high), color="blue", linetype="dashed", size=1)
      print(P1)
    } else {
      print("p.mito is either absent or identical across all cells")
    }
  }

	P1 <- FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
	P1 <- P1 + geom_vline(aes(xintercept=nCount_RNA.low), color="blue", linetype="dashed", size=1)
	P1 <- P1 + geom_vline(aes(xintercept=nCount_RNA.high), color="blue", linetype="dashed", size=1)
	P1 <- P1 + geom_hline(aes(yintercept=nFeature.low), color="blue", linetype="dashed", size=1)
	P1 <- P1 + geom_hline(aes(yintercept=nFeature.high), color="blue", linetype="dashed", size=1)
	print(P1)

  #See: https://github.com/satijalab/seurat/issues/1053#issuecomment-454512002
  expr <- Seurat::FetchData(object = seuratObj, vars = 'nCount_RNA')
  seuratObj <- seuratObj[, which(x = expr >= nCount_RNA.low & expr <= nCount_RNA.high)]
  print(paste0('After nCount_RNA filter: ', length(colnames(x = seuratObj))))

  expr <- Seurat::FetchData(object = seuratObj, vars = 'nFeature_RNA')
  seuratObj <- seuratObj[, which(x = expr >= nFeature.low & expr <= nFeature.high)]
  print(paste0('After nFeature_RNA filter: ', length(colnames(x = seuratObj))))

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
#' @export
#' @return A modified Seurat object.
ScoreCellCycle <- function(seuratObj, min.genes = 10) {
	print('Scoring cell cycle:')

  # We can segregate this list into markers of G2/M phase and markers of S-phase
  s.genes <- .GetSPhaseGenes()
  g2m.genes <- .GetG2MGenes()

  s.genes <- s.genes[which(s.genes %in% rownames(seuratObj))]
  g2m.genes <- g2m.genes[which(g2m.genes %in% rownames(seuratObj))]

  print(paste0("Genes present in seurat object: g2m (", length(g2m.genes), ") and s (", length(s.genes), ")"))

  if (length(g2m.genes) < min.genes || length(s.genes) < min.genes) {
    print("Error, the number of g2m and/or s genes < 5, aborting...")
    return(seuratObj)
  }

  print("Running PCA with cell cycle genes")
  seuratObj <- suppressWarnings(RunPCA(object = seuratObj, reduction.name = 'cc.pca', features = c(s.genes, g2m.genes), do.print = FALSE, verbose = F))
  print(DimPlot(object = seuratObj, reduction = "cc.pca"))

  seuratObj <- CellCycleScoring(object = seuratObj,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = TRUE
  )

  print(
    DimPlot(object = seuratObj, reduction = "cc.pca", dims = c(1, 2)) +
    DimPlot(object = seuratObj, reduction = "cc.pca", dims = c(2, 3)) +
    DimPlot(object = seuratObj, reduction = "cc.pca", dims = c(3, 4)) +
    DimPlot(object = seuratObj, reduction = "cc.pca", dims = c(4, 5)) +
    patchwork::plot_layout(ncol = 2)
  )

	# Discard un-necessary data:
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
#' @param umap.seed Passed directly to Seurat::RunUMAP
#' @param umap.n.epochs Passed directly to Seurat::RunUMAP
#' @param max.tsne.iter The value of max_iter to provide to RunTSNE.  Increasing can help large datasets.
#' @param tsne.perplexity tSNE perplexity. Passed directly to Seurat::RunTSNE(), but CellMembrane:::.InferPerplexityFromSeuratObj() corrects it if need be based dataset dims.
#' @param clusterResolutions A vector of clustering resolutions, default is (0.2, 0.4, 0.6, 0.8, 1.2).
#' @return A modified Seurat object.
#' @export
FindClustersAndDimRedux <- function(seuratObj, dimsToUse = NULL, minDimsToUse = NULL,
                                   umap.method = 'uwot', umap.metric = 'cosine',
                                   umap.n.neighbors = 30L,
																	 umap.min.dist = 0.3, umap.seed = GetSeed(),
                                   umap.n.epochs = NULL, max.tsne.iter = 10000, tsne.perplexity = 30,
                                   clusterResolutions = c(0.2, 0.4, 0.6, 0.8, 1.2) ){
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

  seuratObj <- FindNeighbors(object = seuratObj, dims = dimsToUse, verbose = FALSE)

  for (resolution in clusterResolutions){
    seuratObj <- FindClusters(object = seuratObj, resolution = resolution, verbose = FALSE)
    seuratObj[[paste0("ClusterNames_", resolution)]] <- Idents(object = seuratObj)
  }

  perplexity <- .InferPerplexityFromSeuratObj(seuratObj, perplexity = tsne.perplexity)
  seuratObj <- RunTSNE(object = seuratObj,
									dims.use = dimsToUse,
									check_duplicates = FALSE,
									perplexity = perplexity,
									reduction.key = 'rnaTSNE_',
									max_iter = max.tsne.iter)

  seuratObj <- RunUMAP(seuratObj,
                  dims = dimsToUse,
                  n.neighbors = umap.n.neighbors,
                  min.dist = umap.min.dist,
                  metric = umap.metric,
                  umap.method = umap.method,
									reduction.key = 'rnaUMAP_',
									seed.use = umap.seed,
									n.epochs = umap.n.epochs, verbose = FALSE)

  for (reduction in c('tsne', 'umap')){
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


#' @title Find_Markers
#' @param seuratObj A seurat object
#' @param identFields A vector of grouping fields for DE. Often these are the resolution, computed during FindClustersAndDimRedux()
#' @param outFile A file where a table of markers will be saved
#' @param testsToUse A vector of tests to be used.  Each will be used to run FindAllMarkers() and the results merged. Available are: wilcox, bimod, roc, t, negbinom, poisson, LR, MAST, DESeq2
#' @param numGenesToPrint The number of top markers per cluster to print in a table
#' @param onlyPos If true, only positive markers will be saved
#' @param pValThreshold Only genes with adjusted p-values below this will be reported
#' @param foldChangeThreshold Only genes with log2 fold-change above this will be reported
#' @return A DT::datatable object with the top markers, suitable for printing
#' @importFrom dplyr %>% coalesce group_by summarise filter top_n select everything
#' @import DESeq2
#' @import MAST
#' @export
Find_Markers <- function(seuratObj, identFields, outFile = NULL, testsToUse = c('wilcox', 'MAST', 'DESeq2'), numGenesToPrint = 20, onlyPos = F, pValThreshold = 0.001, foldChangeThreshold = 0.5) {
	seuratObj.markers <- NULL
	for (fieldName in identFields) {
		# Allow resolution to be passed directly:
		if (!(fieldName %in% names(seuratObj@meta.data))) {
  		toTest <- paste0('ClusterNames_', fieldName)
			if (toTest %in% names(seuratObj@meta.data)) {
  			fieldName <- toTest
			}
		}

  	print(paste0('Grouping by field: ', fieldName))
		Idents(seuratObj) <- fieldName

		for (test in testsToUse) {
			print(paste0('Running using test: ', test))
			tryCatch({
				tMarkers <- FindAllMarkers(object = seuratObj, only.pos = onlyPos, min.pct = 0.25, logfc.threshold = 0.25, verbose = F, test.use = test)
				if (nrow(tMarkers) == 0) {
					print('No genes returned, skipping')
				} else {
					tMarkers$test <- c(test)
					tMarkers$groupField <- c(fieldName)
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

	if (all(is.null(seuratObj.markers))) {
		print('All tests failed, no markers returned')
		return()
	}
  else if (nrow(seuratObj.markers) == 0) {
    print('No significant markers were found')
    return()
  } else {
		seuratObj.markers$cluster <- as.factor(seuratObj.markers$cluster)
    toWrite <- seuratObj.markers %>% filter(p_val_adj < pValThreshold) %>% filter(avg_logFC > foldChangeThreshold)
    if (nrow(toWrite) == 0) {
      print('No significant markers were found')
    } else {
      if (!is.null(outFile)) {
        write.table(toWrite, file = outFile, sep = '\t', row.names = F, quote = F)
      }

      print(DimPlot(object = seuratObj, reduction = 'tsne'))

      topGene <- toWrite %>% group_by(groupField, cluster, test) %>% top_n(numGenesToPrint, avg_logFC)
      print(DoHeatmap(object = seuratObj, features = unique(as.character(topGene$gene)), slot = 'data'))

      #Note: return the datatable, so it will be printed correctly by Rmarkdown::render()
  		return(DT::datatable(topGene,
        caption = 'Top DE Genes',
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
      elbowX = 2
    } else {
      elbowX = ndims
    }
  }

  plot <- ggplot(data = data.frame(dims = 1:ndims, stdev = data.use[1:ndims])) +
    geom_point(mapping = aes_string(x = "dims", y = "stdev")) +
    labs(x = gsub(pattern = "_$", replacement = "", x = Key(object = object[[reduction]])),
    y = "Standard Deviation") + theme_bw() + geom_vline(xintercept = elbowX) + ggtitle("Elbow Identification")

  if (print.plot) {
    print(plot)
  }

  return(elbowX)
}


#' @importFrom stats coef
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
    ans <- NULL
    ix <- iy <- 0   # intersecting point
    lineMag <- lineMagnitude(x1, y1, x2, y2)
    if (any(lineMag < 0.00000001)) {
      #warning("short segment")
      #return(9999)
      warning("At least one line segment given by x1, y1, x2, y2 is very short.")
    }
    u <- (((px - x1) * (x2 - x1)) + ((py - y1) * (y2 - y1)))
    u <- u / (lineMag * lineMag)

    ans <- c()
    for (i in 1:length(u)) {
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
  DF <- data.frame(x = 1:length(y), y = y)
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
        min.x = min(DF[which(DF$y<=min.y), ]$x)
      } else {
        print("min.y greater than smallest y")
        min.x = 2
      }
    } else {
      print("min.x and min.y are NA")
      min.x = 2
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






