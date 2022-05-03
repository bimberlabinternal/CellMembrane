utils::globalVariables(
	names = c('x', 'y'),
	package = 'CellMembrane',
	add = TRUE
)



#' @title Wrapper around Seurat::CreateSeuratObject
#'
#' @description Create Seurat object from count data (usually from Read10X()). This also sets pct.mito.
#' @param seuratData Seurat input data, usually from Read10X().
#' @param datasetId This will be used as a prefix for barcodes, and stored in metadata. Also used as the project name for the Seurat object.
#' @param datasetName An optional print-friendly name that will be stored in metadata
#' @param minFeatures Include cells where at least this many features are detected.
#' @param minCells Include features detected in at least this many cells.
#' @param mitoGenesPattern The expression to use when identifying mitochondrial genes
#' @param annotateMitoFromReference If true, a list of mitochondrial genes, taken from (https://en.wikipedia.org/wiki/Category:Human_mitochondrial_genes) will be used to calculate p.mito
#' @return A Seurat object with p.mito calculated.
#' @export
#' @importFrom Matrix colSums
CreateSeuratObj <- function(seuratData, datasetId, datasetName = NULL, minFeatures = 25, minCells = 0, mitoGenesPattern = "^MT-", annotateMitoFromReference = TRUE){
	seuratObj <- Seurat::CreateSeuratObject(counts = seuratData, min.cells = minCells, min.features = minFeatures, project = datasetId)
	seuratObj <- .PossiblyAddBarcodePrefix(seuratObj, datasetId = datasetId, datasetName = datasetName)
	seuratObj<- CalculatePercentMito(seuratObj, mitoGenesPattern = mitoGenesPattern, annotateMitoFromReference = annotateMitoFromReference)

	return(seuratObj)
}

#' @title Calculate Mitochrondrial Percentage
#'
#' @description This will identify mitochrondrial genes and calculate p.mito for each cell
#' @param seuratObj The seurat object
#' @param mitoGenesPattern The expression to use when identifying mitochondrial genes
#' @param annotateMitoFromReference If true, a list of mitochondrial genes, taken from (https://www.genedx.com/wp-content/uploads/crm_docs/Mito-Gene-List.pdf) will be used to calculate p.mito
#' @param outputColName The name of the output column to hold p.mito
#' @return A Seurat object with p.mito calculated.
#' @export
CalculatePercentMito <- function(seuratObj, mitoGenesPattern = "^MT-", annotateMitoFromReference = TRUE, outputColName = 'p.mito') {
	mito.features <- NULL
	if (!annotateMitoFromReference) {
		mito.features <- grep(pattern = mitoGenesPattern, x = rownames(x = seuratObj), value = TRUE)
	} else {
		mito.features <- c('ATP6','ATP8','COX1','COX2','COX3','CYTB','ND1','ND2','ND3','ND4','ND4L','ND5','ND6')

		mito.features.prefix <- c('MT-', mito.features)
		i1 <- length(intersect(mito.features, rownames(seuratObj)))
		i2 <- length(intersect(mito.features.prefix, rownames(seuratObj)))
		if (i2 > i1) {
			print('Selecting using reference gene set with MT- prefix')
			mito.features <- mito.features.prefix
		}
	}

	print(paste0('Total mito features: ', length(mito.features)))
	mito.features <- intersect(mito.features, rownames(seuratObj))
	print(paste0('Total intersecting with seurat rownames (total: ', length(rownames(seuratObj)),'): ', length(mito.features)))

	if (all(is.null(mito.features)) || length(mito.features) == 0) {
		print('No mito features found')
		seuratObj[[outputColName]] <- 0
	} else {
		p.mito <- Matrix::colSums(x = GetAssayData(object = seuratObj, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seuratObj, slot = 'counts'))
		seuratObj[[outputColName]] <- p.mito
	}

	return(seuratObj)
}


#' @importFrom Matrix colSums
.PrintQcPlots <- function(seuratObj) {
	if ('p.mito' %in% colnames(seuratObj@meta.data)) {
		totalPMito <- length(unique(seuratObj$p.mito))
	} else {
		totalPMito <- -1
	}

	feats <- c("nFeature_RNA", "nCount_RNA")
	if (totalPMito > 1) {
		feats <- c(feats, "p.mito")
	}

	print(VlnPlot(object = seuratObj, features = feats, ncol = length(feats)))

	if (totalPMito > 1) {
		print(FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", feature2 = "p.mito"))
	} else {
		print("p.mito absent or identical across all cells, will not plot")
	}
	print(FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))

	#10x-like plot
	nUMI <- Matrix::colSums(GetAssayData(object = seuratObj, slot = "counts"))
	nUMI <- sort(nUMI)

	countAbove <-unlist(lapply(nUMI, function(x){
		sum(nUMI >= x)
	}))

	print(ggplot(data.frame(x = log10(countAbove), y = log(nUMI)), aes(x = x, y = y)) +
		geom_point() + ylab("UMI/Cell") + xlab("log10(# Cells)") +
		egg::theme_presentation()
	)
}


#' @title PerformEmptyDropletFiltering
#'
#' @param seuratRawData Raw data
#' @param fdrThreshold FDR threshold, passed directly to PerformEmptyDrops()
#' @param emptyDropNIters Number of iterations, passed directly to PerformEmptyDrops()
#' @param emptyDropsLower Passed directly to emptyDrops(). The lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets.
#' @param useEmptyDropsCellRanger If TRUE, will use DropletUtils emptyDropsCellRanger instead of emptyDrops
#' @param nExpectedCells Only applied if emptyDropsCellRanger is selected. Passed to n.expected.cells argument
#' @return Plot
#' @importFrom DropletUtils barcodeRanks
PerformEmptyDropletFiltering <- function(seuratRawData, fdrThreshold=0.001, emptyDropNIters=10000, emptyDropsLower=200, useEmptyDropsCellRanger = FALSE, nExpectedCells = 8000) {
	br.out <- DropletUtils::barcodeRanks(seuratRawData)
	plot(br.out$rank, br.out$total+1, log="xy", xlab="Rank", ylab="Total")

	o <- order(br.out$rank)
	lines(br.out$rank[o], br.out$fitted[o], col="red")
	abline(h=br.out@metadata$knee, col="dodgerblue", lty=2)
	abline(h=br.out@metadata$inflection, col="forestgreen", lty=2)
	legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), legend=c("knee", "inflection"))

	print(paste0('Knee: ', br.out@metadata$knee))
	print(paste0('Inflection: ', br.out@metadata$inflection))
	e.out <- PerformEmptyDrops(seuratRawData, emptyDropNIters = emptyDropNIters, fdrThreshold = fdrThreshold, emptyDropsLower = emptyDropsLower, useEmptyDropsCellRanger = useEmptyDropsCellRanger, nExpectedCells = nExpectedCells)

	toPlot <- e.out[is.finite(e.out$LogProb),]
	if (nrow(toPlot) > 0) {
		plot(toPlot$Total, -toPlot$LogProb, col=ifelse(toPlot$is.cell, "red", "black"), log = "x", xlab="log(Total UMI count)", ylab="-Log Probability")
	} else {
		print('Probabilities all -Inf, unable to plot')
	}

	if (nrow(toPlot) != nrow(e.out)) {
		print(paste0('Total rows with non-finite probabilities: ', (nrow(e.out) - nrow(toPlot))))
	}

	print(paste0('Min UMI count in a droplet called a cell: ', min(e.out$Total[e.out$is.cell])))
  	print(paste0('Max UMI count in a droplet not called a cell: ', max(e.out$Total[!e.out$is.cell])))

	passingCells <- rownames(e.out)[e.out$is.cell]

	return(seuratRawData[,passingCells])
}

PerformEmptyDrops <- function(seuratRawData, emptyDropNIters, fdrThreshold=0.001, emptyDropsLower = 100, useEmptyDropsCellRanger = FALSE, nExpectedCells = 8000, seed = GetSeed()){
	print(paste0('Performing ', ifelse(useEmptyDropsCellRanger, yes = 'emptyDropsCellRanger', no = 'emptyDrops'), ' with ', emptyDropNIters, ' iterations'))

	if (!is.null(seed)) {
		set.seed(seed)
	}

	if (useEmptyDropsCellRanger) {
		e.out <- DropletUtils::emptyDropsCellRanger(seuratRawData, niters = emptyDropNIters, n.expected.cells = nExpectedCells)
	} else {
		e.out <- DropletUtils::emptyDrops(seuratRawData, niters = emptyDropNIters, lower = emptyDropsLower)
	}

	print(paste0('Input cells: ', nrow(e.out)))
	badRows <- !is.finite(e.out$FDR)
	if (sum(badRows) > 0) {
		print(paste0('Cells with non-finite FDR: ', sum(badRows)))
		e.out <- e.out[!badRows,]
	}

	e.out$is.cell <- e.out$FDR <= fdrThreshold
	print(paste0('Cells passing FDR (', fdrThreshold, '): ', sum(e.out$is.cell, na.rm=TRUE)))
	print(paste0('Cells failing FDR: ', sum(!e.out$is.cell, na.rm=TRUE)))

	totalLimited <- 0
	if (!useEmptyDropsCellRanger) {
		#If there are any entries with FDR above the desired threshold and Limited==TRUE, it indicates that npts should be increased in the emptyDrops call.
		print(table(Limited=e.out$Limited, Significant=e.out$is.cell))
		totalLimited <- sum(e.out$Limited == T & e.out$Significant == F)
	}

	if (totalLimited == 0){
		return(e.out)
	} else {
		print('Repeating emptyDrops with more iterations')
		return(PerformEmptyDrops(seuratRawData, emptyDropNIters = emptyDropNIters * 2, fdrThreshold = fdrThreshold, emptyDropsLower = emptyDropsLower, seed = seed, useEmptyDropsCellRanger = useEmptyDropsCellRanger, nExpectedCells = nExpectedCells))
	}
}


.DoMergeSimple <- function(seuratObjs, projectName, merge.data = FALSE, expectedDefaultAssay = NULL, doGC = FALSE){
	seuratObj <- NULL

	for (datasetId in names(seuratObjs)) {
		if (is.null(seuratObj)) {
			seuratObj <- seuratObjs[[datasetId]]
			if (!is.null(expectedDefaultAssay)) {
				DefaultAssay(seuratObj) <- expectedDefaultAssay
			}
		} else {
			assayName <- DefaultAssay(seuratObj)
			DefaultAssay(seuratObjs[[datasetId]]) <- assayName

			hasGeneId <- 'GeneId' %in% names(GetAssay(seuratObjs[[datasetId]], assay = assayName)@meta.features)

			if (any(rownames(seuratObj[[assayName]]) != rownames(seuratObjs[[datasetId]][[assayName]]))) {
				missing <- rownames(seuratObj[[assayName]])[!(rownames(seuratObj[[assayName]]) %in% rownames(seuratObjs[[datasetId]][[assayName]]))]
				missing <- c(missing, rownames(seuratObjs[[datasetId]][[assayName]])[!(rownames(seuratObjs[[datasetId]][[assayName]]) %in% rownames(seuratObj[[assayName]]))])
				stop(paste0('Gene names are not equal! Missing: ', paste0(missing, collapse = ',')))
			}

			if (hasGeneId) {
				# This can occur if the first object lacks GeneId but the second has it.
				if (!'GeneId' %in% names(GetAssay(seuratObj, assay = assayName)@meta.features)) {
					message('GeneId was present in the second merged object, but not the first. Adding NAs')
					seuratObj@assays[[assayName]]@meta.features$GeneId <- NA
				}

				geneIds1 <- GetAssay(seuratObj, assay = assayName)@meta.features$GeneId
				geneIds2 <- GetAssay(seuratObjs[[datasetId]], assay = assayName)@meta.features$GeneId
				names(geneIds1) <- rownames(seuratObj[[assayName]])
				names(geneIds2) <- rownames(seuratObjs[[datasetId]][[assayName]])
			}

			seuratObj <- merge(x = seuratObj, y = seuratObjs[[datasetId]], project = projectName, merge.data = merge.data)
			seuratObjs[[datasetId]] <- NULL
			if (doGC) {
				gc()
			}

			if (hasGeneId) {
				if (any(is.na(geneIds1)) & !any(is.na(geneIds2))) {
					seuratObj[[assayName]] <- AddMetaData(seuratObj[[assayName]], metadata = geneIds2, col.name = 'GeneId')
				} else if (!any(is.na(geneIds1)) & any(is.na(geneIds2))) {
					seuratObj[[assayName]] <- AddMetaData(seuratObj[[assayName]], metadata = geneIds1, col.name = 'GeneId')
				} else if (!any(is.na(geneIds1)) & !any(is.na(geneIds2))) {
					if (any(geneIds1 != geneIds2)) {
						stop('Gene IDs did not match between seurat objects!')
					}
					seuratObj[[assayName]] <- AddMetaData(seuratObj[[assayName]], metadata = geneIds1, col.name = 'GeneId')
				} else {
					seuratObj[[assayName]] <- AddMetaData(seuratObj[[assayName]], metadata = geneIds1, col.name = 'GeneId')
				}
			}
		}
	}

	return(seuratObj)
}

.RunSoupX <- function(outsDir) {
	sc <- SoupX::load10X(outsDir)
	sc <- SoupX::autoEstCont(sc, doPlot = T)
	seuratRawData <- SoupX::adjustCounts(sc)
	print(paste0('Passing cells: ', ncol(seuratRawData)))
	print(SoupX::plotMarkerMap(sc, geneSet = sc$fit$markersUsed$gene))
	print(SoupX::plotChangeMap(sc, cleanedMatrix = seuratRawData, geneSet = sc$fit$markersUsed$gene))

	# Drop suffixes:
	colnames(seuratRawData) <- sapply(colnames(seuratRawData), function(x){
		return(unlist(strsplit(x, split = '-'))[1])
	})

	return(seuratRawData)
}

.InferMatrixDir <- function(dataDir) {
	matrixFile <- paste0(dataDir, 'matrix.mtx')
	matrixFileGz <- paste0(dataDir, 'matrix.mtx.gz')
	if (file.exists(matrixFile) || file.exists(matrixFileGz)) {
		return(dataDir)
	}

	dirWithFeatureMatrix <- paste0(dataDir, 'raw_feature_bc_matrix')
	matrixFile <- paste0(dirWithFeatureMatrix, '/matrix.mtx')
	matrixFileGz <- paste0(dirWithFeatureMatrix, '/matrix.mtx.gz')
	if (file.exists(matrixFile) || file.exists(matrixFileGz)) {
		return(dirWithFeatureMatrix)
	}

	stop(paste0('Unable to find matrix file in: ', dataDir, ' or ', dirWithFeatureMatrix))
}