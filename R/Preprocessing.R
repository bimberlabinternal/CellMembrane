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
#' @param annotateMitoFromReference If true, a list of mitochondrial genes, taken from (https://www.genedx.com/wp-content/uploads/crm_docs/Mito-Gene-List.pdf) will be used to calculate p.mito
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
		mito.features <- CellMembrane::mitoGenes$Gene
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
		totalPMito = length(unique(seuratObj$p.mito))
	} else {
		totalPMito = -1
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

	print(ggplot(data.frame(x = log(countAbove), y = log(nUMI)), aes(x = x, y = y)) +
		geom_point() + ylab("UMI/Cell") + xlab("# Cells") +
		egg::theme_presentation()
	)
}


#' @title PerformEmptyDropletFiltering
#'
#' @param seuratRawData Raw data
#' @param fdrThreshold FDR threshold, passed directly to PerformEmptyDrops()
#' @param emptyDropNIters Number of iterations, passed directly to PerformEmptyDrops()
#' @param emptyDropsLower Passed directly to emptyDrops(). The lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets.
#' @return Plot
#' @importFrom DropletUtils barcodeRanks
PerformEmptyDropletFiltering <- function(seuratRawData, fdrThreshold=0.01, emptyDropNIters=10000, emptyDropsLower=100) {
	br.out <- DropletUtils::barcodeRanks(seuratRawData)

	# Making a plot.
	plot(br.out$rank, br.out$total+1, log="xy", xlab="Rank", ylab="Total")

	o <- order(br.out$rank)
	lines(br.out$rank[o], br.out$fitted[o], col="red")
	abline(h=br.out$knee, col="dodgerblue", lty=2)
	abline(h=br.out$inflection, col="forestgreen", lty=2)
	legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), legend=c("knee", "inflection"))

	e.out <- PerformEmptyDrops(seuratRawData, emptyDropNIters = emptyDropNIters, fdrThreshold = fdrThreshold, emptyDropsLower = emptyDropsLower)

	toPlot <- e.out[is.finite(e.out$LogProb),]
	if (nrow(toPlot) > 0) {
		plot(toPlot$Total, -toPlot$LogProb, col=ifelse(toPlot$is.cell, "red", "black"), xlab="Total UMI count", ylab="-Log Probability")
	} else {
		print('Probabilities all -Inf, unable to plot')
	}

	if (nrow(toPlot) != nrow(e.out)) {
		print(paste0('Total rows with non-finite probabilities: ', (nrow(e.out) - nrow(toPlot))))
	}

	passingCells <- rownames(e.out)[e.out$is.cell]

	return(seuratRawData[,passingCells])
}

PerformEmptyDrops <- function(seuratRawData, emptyDropNIters, fdrThreshold=0.01, emptyDropsLower = 100, seed = GetSeed()){
	print(paste0('Performing emptyDrops with ', emptyDropNIters, ' iterations'))

	if (!is.null(seed)) {
		set.seed(seed)
	}

	e.out <- DropletUtils::emptyDrops(seuratRawData, niters = emptyDropNIters, lower = emptyDropsLower)

	print(paste0('Input cells: ', nrow(e.out)))
	e.out <- e.out[!is.na(e.out$LogProb),]
	e.out$is.cell <- e.out$FDR <= fdrThreshold
	print(paste0('Cells passing FDR: ', sum(e.out$is.cell, na.rm=TRUE)))
	print(paste0('Cells failing FDR: ', sum(!e.out$is.cell, na.rm=TRUE)))

	#If there are any entries with FDR above the desired threshold and Limited==TRUE, it indicates that npts should be increased in the emptyDrops call.
	print(table(Limited=e.out$Limited, Significant=e.out$is.cell))
	totalLimited <- sum(e.out$Limited[e.out$Limited == T] & e.out$Significant == F)
	if (totalLimited == 0){
		return(e.out)
	} else {
		return(PerformEmptyDrops(seuratRawData, emptyDropNIters = emptyDropNIters * 2, fdrThreshold = fdrThreshold, emptyDropsLower = emptyDropsLower, seed = seed))
	}
}


.DoMergeSimple <- function(seuratObjs, projectName, merge.data = FALSE){
	seuratObj <- NULL

	for (datasetId in names(seuratObjs)) {
		if (is.null(seuratObj)) {
			seuratObj <- seuratObjs[[datasetId]]
		} else {
			assayName <- DefaultAssay(seuratObj)
			hasGeneId = ifelse(is.null(GetAssay(seuratObjs[[datasetId]])@meta.features$GeneId), F, T)

			if (any(rownames(seuratObj) != rownames(seuratObjs[[datasetId]]))) {
				missing <- rownames(seuratObj)[!(rownames(seuratObj) %in% rownames(seuratObjs[[datasetId]]))]
				missing <- c(missing, rownames(seuratObjs[[datasetId]])[!(rownames(seuratObjs[[datasetId]]) %in% rownames(seuratObj))])
				stop(paste0('Gene names are not equal! Missing: ', paste0(missing, collapse = ',')))
			}

			if (hasGeneId) {
				geneIds1 <- GetAssay(seuratObj)@meta.features$GeneId
				geneIds2 <- GetAssay(seuratObjs[[datasetId]])@meta.features$GeneId
				names(geneIds1) <- rownames(seuratObj)
				names(geneIds2) <- rownames(seuratObjs[[datasetId]])
			}

			seuratObj <- merge(x = seuratObj, y = seuratObjs[[datasetId]], project = projectName, merge.data = merge.data)

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

