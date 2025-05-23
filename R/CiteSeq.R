#' @include Utils.R
#' @include Preprocessing.R
#' @import Seurat

utils::globalVariables(
	names = c('TotalCount', 'Marker', 'sortorder'),
	package = 'CellMembrane',
	add = TRUE
)

#' @title Read/Append CITE-seq/ADT data to a Seurat Object
#'
#' @description Read/Append CITE-seq/ADT data to a Seurat Object
#' @param seuratObj The seurat object where data will be added.
#' @param unfilteredMatrixDir The directory holding raw count data, generally the raw_feature_bc_matrix from the cellranger outs folder
#' @param normalizeMethod The method for normalization. If NULL, normalization will be skipped.
#' @param datasetId If the seurat object includes cells from multiple datasets, and/or if the barcodes were prefixed with a dataset id, This will be used as a prefix for barcodes, and stored in metadata. Also used as the project name for the Seurat object.
#' @param assayName The name of the assay to store ADT data.
#' @param featureMetadata If is common for the raw barcode to use non-human friendly names (like TotalSeq-C-XXXX). This is an optional dataframe that can be used to override that, and also provide other feature metadata. This dataframe must have rownames that match the rownames of the saved matrix. If this dataframe has the column 'name', this will be used to replace the original rownames on the ADT assay.
#' @param adtWhitelist An optional whitelist of ADT names (matching the raw names from the matrix). If provided, the matrix will be subset to just these features
#' @param minRowSum If provided, any ADTs (rows) with rowSum below this value will be dropped.
#' @param failIfAdtsInWhitelistNotFound If an adtWhitelist is provided and this is TRUE, an error will be thrown if any of these features are missing in the input matrix
#' @param aggregateBarcodeFile Optional. This is the cellranger output, in antibody_analysis/aggregate_barcodes.csv, which contains barcodes marked as aggregates. These are dropped.
#' @export
#' @importFrom dplyr arrange
AppendCiteSeq <- function(seuratObj, unfilteredMatrixDir, normalizeMethod = 'dsb', datasetId = NULL, assayName = 'ADT', featureMetadata = NULL, adtWhitelist = NULL, minRowSum = NULL, failIfAdtsInWhitelistNotFound = TRUE, aggregateBarcodeFile = NULL) {
	print(paste0('Initial cell barcodes in GEX data: ', ncol(seuratObj)))
	if (!is.null(datasetId)) {
		gexCells <- colnames(seuratObj)[seuratObj$DatasetId == datasetId]
		print(paste0('Initial cell barcodes in GEX data for prefix: ', length(gexCells)))
	} else {
		gexCells <- colnames(seuratObj)
	}
	assayData <- .LoadCiteSeqData(unfilteredMatrixDir, datasetId = datasetId, featureMetadata = featureMetadata, adtWhitelist = adtWhitelist, minRowSum = minRowSum, failIfAdtsInWhitelistNotFound = failIfAdtsInWhitelistNotFound)

	if (!is.null(aggregateBarcodeFile)) {
		barcodes <- read.table(aggregateBarcodeFile, header = T, sep = ',')$barcode
		barcodes <- unname(sapply(barcodes, function(x){
			return(unlist(strsplit(x, split = '-'))[1])
		}))
		barcodes <- paste0(datasetId, '_', barcodes)
		print(paste0('Total barcodes marked as aggregates: ', length(barcodes)))
		if (sum(!barcodes %in% colnames(assayData)) > 0) {
			print(paste0('The following barcodes were not in the count matrix: ', paste0(barcodes[!barcodes %in% colnames(assayData)], collapse = ','), ', first was: ', colnames(assayData)[1]))
		}

		toKeep <- colnames(assayData)[!colnames(assayData) %in% barcodes]
		assayData <- subset(assayData, cells = toKeep)
		if (ncol(assayData) != length(toKeep)) {
			stop(paste0('Incorrect assay subset. Expected: ', length(toKeep), ', actual: ', ncol(assayData)))
		}
		print(paste0('After removing: ', ncol(assayData)))
	}

	sharedCells <- colnames(assayData)[which(colnames(assayData) %in% gexCells)]
	print(paste0('Intersect with GEX data: ', length(sharedCells)))

	if (length(sharedCells) == 0) {
		print(head(gexCells))
		print(head(colnames(assayData)))
		stop('No cells are shared')
	}

	if (is.null(normalizeMethod)){
		print('Normalization will not be performed')
		assayData <- subset(assayData, cells = sharedCells)
		if (ncol(assayData) != length(sharedCells)) {
			stop(paste0('Incorrect assay subset. Expected: ', length(sharedCells), ', actual: ', ncol(assayData)))
		}

	} else if (normalizeMethod == 'dsb') {
		assayData <- .NormalizeDsbWithEmptyDrops(seuratObj, assayData)
	} else if (normalizeMethod == 'clr') {
		assayData <- subset(assayData, cells = sharedCells)
		if (ncol(assayData) != length(sharedCells)) {
			stop(paste0('Incorrect assay subset. Expected: ', length(sharedCells), ', actual: ', ncol(assayData)))
		}

		assayData <- Seurat::NormalizeData(assayData, normalization.method = 'CLR', margin = 2, verbose = FALSE)
	} else {
		stop('Unknown normalizationMethod. Pass NULL to skip normalization')
	}

	#Note: we need to support merging with any existing data and account for datasetId:
	if (assayName %in% names(seuratObj@assays)) {
		print('Existing assay found, merging counts')
		assayData <- .MergeAdtWithExisting(assayData, GetAssay(seuratObj, assay = assayName))
	} else {
		assayData <- .PossiblyExpandAssay(seuratObj, assayData)
	}

	if (assayName %in% names(seuratObj@assays)) {
		seuratObj[[assayName]] <- NULL  # Reset this first to avoid warning
	}

	seuratObj[[assayName]] <- assayData
	.PlotMarkerQc(seuratObj, assayName = assayName)

	return(seuratObj)
}

.PossiblyExpandAssay <- function(seuratObj, assayData) {
	allCells <- unique(c(colnames(seuratObj), colnames(assayData)))
	assayData <- .EnsureCellsPresentInOrder(assayData, allCells)

	if (sum(colnames(seuratObj) != colnames(assayData)) > 0) {
		stop('The columns of the ADT matrix do not equal the seurat object')
	}

	return(assayData)
}

.EnsureFeaturesPresentInOrder <- function(assayData, featureWhitelist) {
	featuresToAdd <- featureWhitelist[!(featureWhitelist %in% rownames(assayData))]
	if (length(featuresToAdd) == 0) {
		return(.ReorderAssayFeatures(assayData, featureWhitelist))
	}

	replacementAssay <- NULL
	for (layer in c('counts', 'data')) {
		layerData <- GetAssayData(assayData, layer = layer)
		if (is.null(layerData)) {
			next
		}

		# Add any new ADTs from this dataset, if needed:
		missingMat <- matrix(rep(0, ncol(layerData) * length(featuresToAdd)), ncol = ncol(layerData), nrow = length(featuresToAdd))
		rownames(missingMat) <- featuresToAdd
		print(paste0('total ADT rows added to assay: ', length(featuresToAdd)))
		layerData <- Seurat::as.sparse(rbind(layerData, missingMat))
		layerData <- layerData[featureWhitelist, ]

		if (is.null(replacementAssay)) {
			args <- list()
			args[[layer]] <- layerData
			replacementAssay <- rlang::invoke(Seurat::CreateAssayObject, args)
		} else {
			replacementAssay <- SetAssayData(object = replacementAssay, layer = layer, new.data = layerData)
		}
	}

	return(replacementAssay)
}

.ReorderAssayFeatures <- function(assayData, featureWhitelist) {
	if (nrow(assayData) != length(featureWhitelist)) {
		stop("Feature set does not match, cannot reorder")
	}

	if (length(Seurat::GetAssayData(assayData, layer = 'counts'))) {
		assayData <- Seurat::SetAssayData(assayData, layer = 'counts', new.data = Seurat::GetAssayData(assayData, layer = 'counts')[featureWhitelist, ])
	}

	if (length(Seurat::GetAssayData(assayData, layer = 'data'))) {
		assayData <- Seurat::SetAssayData(assayData, layer = 'data', new.data = Seurat::GetAssayData(assayData, layer = 'data')[featureWhitelist, ])
	}

	if (length(Seurat::GetAssayData(assayData, layer = 'scale.data'))) {
		assayData <- Seurat::SetAssayData(assayData, layer = 'scale.data', new.data = Seurat::GetAssayData(assayData, layer = 'scale.data')[featureWhitelist, ])
	}

	slot(assayData, GetAssayMetadataSlotName(assayData)) <- slot(assayData, GetAssayMetadataSlotName(assayData))[featureWhitelist]

	return(assayData)
}

.ReorderAssayCells <- function(assayData, cellWhitelist) {
	if (ncol(assayData) != length(cellWhitelist)) {
		stop("Cell list does not match, cannot reorder")
	}

	if (length(Seurat::GetAssayData(assayData, layer = 'counts'))) {
		assayData <- Seurat::SetAssayData(assayData, layer = 'counts', new.data = Seurat::GetAssayData(assayData, layer = 'counts')[, cellWhitelist])
	}

	if (length(Seurat::GetAssayData(assayData, layer = 'data'))) {
		assayData <- Seurat::SetAssayData(assayData, layer = 'data', new.data = Seurat::GetAssayData(assayData, layer = 'data')[, cellWhitelist])
	}

	if (length(Seurat::GetAssayData(assayData, layer = 'scale.data'))) {
		assayData <- Seurat::SetAssayData(assayData, layer = 'scale.data', new.data = Seurat::GetAssayData(assayData, layer = 'scale.data')[, cellWhitelist])
	}

	return(assayData)
}

.EnsureCellsPresentInOrder <- function(assayData, cellWhitelist) {
	cellsToAdd <- cellWhitelist[!(cellWhitelist %in% colnames(assayData))]
	if (length(cellsToAdd) == 0) {
		return(.ReorderAssayCells(assayData, cellWhitelist))
	}

	replacementAssay <- NULL
	for (layer in c('counts', 'data')) {
		layerData <- GetAssayData(assayData, layer = layer)
		if (is.null(layerData)) {
			next
		}

		# Add any new cells from this dataset, if needed:
		missingMat <- matrix(rep(0, nrow(layerData) * length(cellsToAdd)), nrow = nrow(layerData), ncol = length(cellsToAdd))
		colnames(missingMat) <- cellsToAdd
		print(paste0('total cells added to assay: ', length(cellsToAdd)))
		layerData <- Seurat::as.sparse(cbind(layerData, missingMat))
		layerData <- layerData[, cellWhitelist]

		if (is.null(replacementAssay)) {
			args <- list()
			args[[layer]] <- layerData
			replacementAssay <- rlang::invoke(Seurat::CreateAssayObject, args)
		} else {
			replacementAssay <- SetAssayData(object = replacementAssay, layer = layer, new.data = layerData)
		}
	}

	# Ensure metadata preserved:
	slot(replacementAssay, GetAssayMetadataSlotName(replacementAssay)) <- slot(assayData, GetAssayMetadataSlotName(assayData))

	return(replacementAssay)
}

.MergeAdtWithExisting <- function(newAssay, existingAssay) {
	allFeatures <- unique(c(rownames(existingAssay), rownames(newAssay)))
	allCells <- unique(c(colnames(existingAssay), colnames(newAssay)))
	existingAssay <- .EnsureFeaturesPresentInOrder(existingAssay, allFeatures)
	existingAssay <- .EnsureCellsPresentInOrder(existingAssay, allCells)

	replacementAssay <- NULL
	for (layer in c('counts', 'data')) {
		data <- GetAssayData(newAssay, layer = layer)
		if (is.null(data)) {
			next
		}

		existingData <- GetAssayData(existingAssay, layer = layer)
		if (is.null(existingData)) {
			existingData <- matrix(rep(0, ncol(data)*nrow(data), nrow = nrow(data), ncol = ncol(data)))
			rownames(existingData) <- rownames(data)
			colnames(existingData) <- colnames(data)
		}

		existingData[rownames(data), colnames(data)] <- data
		if (is.null(replacementAssay)) {
			args <- list()
			args[[layer]] <- Seurat::as.sparse(existingData)
			replacementAssay <- rlang::invoke(Seurat::CreateAssayObject, args)
		} else {
			replacementAssay <- SetAssayData(object = replacementAssay, layer = layer, new.data = Seurat::as.sparse(existingData))
		}
	}

	return(replacementAssay)
}

.LoadCiteSeqData <- function(unfilteredMatrixDir, datasetId = NULL, featureMetadata = NULL, adtWhitelist = NULL, minRowSum = NULL, failIfAdtsInWhitelistNotFound = TRUE) {
	if (!dir.exists(unfilteredMatrixDir)){
		stop("Count matrix not found")
	}

	bData <- Seurat::Read10X(unfilteredMatrixDir, gene.column=1, strip.suffix = TRUE)
	bData <- bData[which(!(rownames(bData) %in% c('unmapped'))), , drop = FALSE]

	#Cannot have underscores in feature names, Seurat will replace with hyphen anyway.  Perform upfront to avoid warning
	if (sum(grepl(x = rownames(bData), pattern = '_')) > 0) {
		print('Replacing underscores with hyphens in feature names')
		rownames(bData) <- gsub(x = rownames(bData), pattern = '_', replacement = '-')
	}

	if (!is.null(datasetId)) {
		colnames(bData) <- paste0(datasetId, '_', colnames(bData))
	}
	bData <- Seurat::as.sparse(bData)

	print(paste0('Initial cells in cite-seq matrix: ', ncol(bData)))

	if (!is.null(minRowSum)) {
		toDrop <- rowSums(as.matrix(bData)) < minRowSum
		if (sum(toDrop) > 0){
			print(paste0('ADTs dropped due to low counts across cells: ', sum(toDrop)))
			print(paste0(rownames(bData)[toDrop], collapse = ','))
			bData <- bData[!toDrop, ]
			print(paste0('ADTs after filter: ', nrow(bData)))
			if (nrow(bData) == 0) {
				stop('No ADTs remain after filter')
			}
		}
	}

	# Apply this prior to rename:
	if (!is.null(adtWhitelist)) {
		print('Filtering ADTs based on whitelist')
		adtWhitelist <- gsub(x = adtWhitelist, pattern = '_', replacement = '-')

		bData <- bData[rownames(bData) %in% adtWhitelist, ]
		print(paste0('ADTs after filter: ', nrow(bData)))
		if (nrow(bData) == 0) {
			stop('No rows left after filter!')
		}

		if (failIfAdtsInWhitelistNotFound) {
			missing <- adtWhitelist[!(adtWhitelist %in% rownames(bData))]
			if (length(missing) > 0) {
				stop(paste0('The following ADTs were requested but not in adtWhitelist: ', paste0(missing, collapse = ','), '. Markers present: ', paste0(rownames(bData), collapse = ',')))
			}
		}
	}

	bData <- Seurat::CreateAssayObject(counts = bData)

	# Optionally rename features
	if (!all(is.null(featureMetadata))) {
		if (!('rowname' %in% colnames(featureMetadata))) {
			stop('featureMetadata must contain the column rowname, which matches the raw matrix')
		}

		print('Renaming ADTs')
		newRows <- data.frame(rowname = rownames(bData), sortorder = 1:nrow(bData), stringsAsFactors = F)
		featureMetadata$rowname <- gsub(x = featureMetadata$rowname, pattern = '_', replacement = '-')
		newRows <- merge(newRows, featureMetadata, by = 'rowname', all.x = T, all.y = F)
		newRows <- newRows %>% arrange(sortorder)
		newRows <- newRows[names(newRows) != 'sortorder']

		newRows$markername <- dplyr::coalesce(newRows$markername, newRows$rowname)
		newRows$markername <- gsub(x = newRows$markername, pattern = '_', replacement = '-')
		d <- duplicated(newRows$markername)
		if (sum(d) > 0) {
			stop('There were duplicate marker names after rename: ' + paste0(newRows$markername[d], collapse = ', '))
		}

		print(paste0('Total renamed: ', sum(newRows$markername != newRows$rowname)))
		bData <- .RenameAssayFeatures(bData, newRows$markername)
		featureMetadata <- newRows
		bData <- .AddFeatureMetadata(bData, featureMetadata)
	}

	return(bData)
}

.RenameAssayFeatures <- function(assayData, newnames) {
	if (nrow(assayData) != length(newnames)) {
		stop("Gene set does not match, cannot rename")
	}

	if (sum(grepl(x = newnames, pattern = '_')) > 0) {
		print('Replacing underscores with hyphens in feature names')
		newnames <- gsub(x = newnames, pattern = '_', replacement = '-')
	}

	if (class(assayData)[1] == 'Assay') {
		for (slotName in c('counts', 'data', 'scaled.data')) {
			if (slotName %in% slotNames(assayData)) {
				print(paste0('Updating layer: ', slotName))
				slot(assayData, slotName)@Dimnames[[1]] <- newnames

				if (any(rownames(Seurat::GetAssayData(assayData, layer = slotName)) != newnames)) {
					stop('Features were not updated!')
				}
			}
		}

		assayData@var.features <- character()
		rownames(assayData@meta.features) <- newnames
	} else if (class(assayData)[1] == 'Assay5') {
		rownames(assayData) <- newnames
	} else {
		stop(paste0('Unknown assay class: ', class(assayData)[1]))
	}

	if (any(rownames(assayData) != newnames)) {
		stop('Features were not updated!')
	}

	return(assayData)
}

.AddFeatureMetadata <- function(assayData, featureMetadata) {
	print('Adding feature metadata')
	for (colname in names(featureMetadata)) {
		if (colname == 'rowname') {
			next
		}

		print(paste0('adding column: ', colname))
		toAdd <- featureMetadata[[colname]]
		names(toAdd) <- featureMetadata$markername
		select <- names(toAdd) %in% rownames(assayData)
		toAdd <- toAdd[select]
		if (length(toAdd) != nrow(featureMetadata)) {
			print(paste0('Not all markers in feature table were present in found table, missing: ', paste0(featureMetadata$markername[!select])))
		}

		if (length(toAdd) > 0) {
			assayData <- Seurat::AddMetaData(object = assayData, metadata = toAdd, col.name = colname)
		}
	}

	return(assayData)
}

.GetNonZeroFeatures <- function(seuratObj, assayName) {
	assayData <- GetAssayData(seuratObj, layer = "counts", assay = assayName)
	featuresToPlot <- rownames(assayData)
	toSkip <- character()
	dat <- Seurat::GetAssayData(seuratObj, assay = assayName, layer = 'counts')
	for (feature in featuresToPlot) {
		if (all(is.na(dat[feature,])) || max(dat[feature,], na.rm = T) == 0) {
			print(paste0('Skipping feature with zero counts: ', feature))
			toSkip <- c(toSkip, feature)
		}
	}
	featuresToPlot <- featuresToPlot[!(featuresToPlot %in% toSkip)]

	return(featuresToPlot)
}

.PlotCiteSeqCountData <- function(seuratObj, assayName = 'ADT') {
	assayData <- GetAssayData(seuratObj, layer = "counts", assay = assayName)

	featuresToPlot <- .GetNonZeroFeatures(seuratObj, assayName)

	setSize <- 2
	steps <- ceiling(length(featuresToPlot) / setSize) - 1

	for (i in 0:steps) {
		start <- (i * setSize) + 1
		end <- min((start + setSize - 1), length(featuresToPlot))
		features <- featuresToPlot[start:end]

		tryCatch({
			print(RidgePlot(seuratObj, assay = assayName, features = features, ncol = 1))
		}, error = function(e){
			print(paste0('Error generating CITE-seq RidgePlot for: ', paste0(features, collapse = ',')))
			print(conditionMessage(e))
			traceback()
		})

	}

	# Also total per ADT
	countsPerAdt <- rowSums(as.matrix(assayData))
	countsPerAdt <- data.frame(Marker = names(countsPerAdt), TotalCount = countsPerAdt)

	P1 <- ggplot(countsPerAdt, aes(x = TotalCount)) +
		geom_density() +
		xlab('Total Count/ADT') +
		ylab('Density') +
		theme(
			axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
		) +
		labs(title = 'Total Counts Per ADT')

	print(P1)

	P2 <- ggplot(countsPerAdt, aes(x = Marker, y = TotalCount)) +
		geom_bar(stat = 'identity') +
		xlab('Marker') +
		ylab('Total Count') +
		theme(
			axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
		) +
		labs(title = 'Total Counts Per ADT')

	print(P2)
}

#' @import patchwork
.NormalizeDsbWithEmptyDrops <- function(seuratObj, unfilteredAdtAssay, rnaAssayName = 'RNA', fdrThreshold=0.01, emptyDropNIters=10000, emptyDropsLower=100) {
	print('Running DSB to normalize data')
	gexCountMatrix <- GetAssayData(object = seuratObj, layer = 'counts', assay = rnaAssayName)
	unfilteredAdtCountMatrix <- GetAssayData(object = unfilteredAdtAssay, layer = 'counts')

	emptyDrops <- PerformEmptyDrops(unfilteredAdtCountMatrix, fdrThreshold = fdrThreshold, emptyDropNIters = emptyDropNIters, emptyDropsLower = emptyDropsLower)
	passingCells <- rownames(emptyDrops)[emptyDrops$is.cell]
	passingCellsWithGex <- intersect(passingCells, colnames(gexCountMatrix))
	print(paste0('ADT passing: ', length(passingCells), ', intersecting with GEX: ', length(passingCellsWithGex)))

	# create a metadata dataframe of simple qc stats for each droplet
	rna_size <- log10(Matrix::colSums(gexCountMatrix))
	prot_size <- log10(Matrix::colSums(unfilteredAdtCountMatrix))

	df1 <- data.frame(rna_size = rna_size, is_cell = ifelse(colnames(gexCountMatrix) %in% passingCells, yes = 'Cell', no = 'Not Cell'))
	p1 <- ggplot(df1[df1$rna_size > 0, ], aes(x = rna_size)) + geom_density(fill = "dodgerblue") + ggtitle("RNA library size \n distribution") + facet_grid(is_cell ~ .)

	df2 <- data.frame(prot_size = prot_size, is_cell = ifelse(colnames(unfilteredAdtCountMatrix) %in% passingCells, yes = 'Cell', no = 'Not Cell'))
	p2 <- ggplot(df2[df2$prot_size> 0, ], aes(x = prot_size)) + geom_density(fill = "firebrick2") + ggtitle("Protein library size \n distribution") + facet_grid(is_cell ~ .)
	print(p1 + p2)

	# define a vector of background / empty droplet barcodes based on protein library size and mRNA content
	background_drops <- !emptyDrops$is.cell
	background_drops <- background_drops[!(background_drops %in% colnames(seuratObj))]
	negative_mtx_rawprot <- unfilteredAdtCountMatrix[ , background_drops] %>% as.matrix()

	# define a vector of cell-containing droplet barcodes based on protein library size and mRNA content
	filteredAdtCountMatrix <- unfilteredAdtCountMatrix[ , passingCellsWithGex] %>% as.matrix()

	print(paste0('ADT cells after filter: ', ncol(filteredAdtCountMatrix)))
	if (ncol(filteredAdtCountMatrix) == 0) {
		stop('No cells passed filters')
	}

	normCounts <- dsb::DSBNormalizeProtein(
		cell_protein_matrix = filteredAdtCountMatrix, # cell containing droplets
		empty_drop_matrix = negative_mtx_rawprot, # estimate ambient noise with the background drops
		denoise.counts = FALSE, # model and remove each cell's technical component
		use.isotype.control = FALSE # use isotype controls to define the technical component
	)

	a <- Seurat::CreateAssayObject(counts = Seurat::as.sparse(filteredAdtCountMatrix[ ,colnames(normCounts)]))
	a <- SetAssayData(object = a, layer = 'data', new.data = Seurat::as.sparse(normCounts))

	return(a)
}

#' @title Perform DimRedux on CiteSeq using euclidean distance
#'
#' @description Perform DimRedux on CiteSeq using euclidean distance
#' @param seuratObj The seurat object
#' @param assayName The name of the assay holding the ADT data.
#' @param dist.method The method, passed to dist()
#' @param performClrNormalization If true, Seurat's CLR normalization will be performed. Otherwise this expected data to be pre-normalized
#' @param print.plots If true, QC plots will be printed
#' @param doUMAP If true, RunUMAP will be performed on the ADT distance matrix.
#' @export
#' @importFrom dplyr arrange
#' @import Seurat
CiteSeqDimRedux.Dist <- function(seuratObj, assayName = 'ADT', dist.method = "euclidean", print.plots = TRUE, performClrNormalization = TRUE, doUMAP = TRUE){
	origAssay <- DefaultAssay(seuratObj)
	DefaultAssay(seuratObj) <- assayName
	print(paste0('Processing ADT data, features: ', paste0(rownames(seuratObj), collapse = ',')))

	# Before we recluster the data on ADT levels, we'll stash the original cluster IDs for later
	seuratObj[["origClusterID"]] <- Idents(seuratObj)

	if (performClrNormalization) {
		seuratObj <- NormalizeData(seuratObj, assay = assayName, normalization.method = 'CLR', margin = 2, verbose = FALSE)
	} else if (is.null(Seurat::GetAssayData(seuratObj, assay = assayName, layer = 'data')) || length(Seurat::GetAssayData(seuratObj, assay = assayName, layer = 'data')) == 0){
		stop('Cannot use performClrNormalization=FALSE without pre-existing ADT normalization')
	} else {
		print('Using pre-existing normalization')
	}

	#SNN:
	print("Calculating Distance Matrix")
	adt.data <- GetAssayData(seuratObj, assay = assayName, layer = "data")
	adt.dist <- dist(Matrix::t(adt.data), method = dist.method)
	seuratObj[["adt_snn.dist"]]  <- FindNeighbors(adt.dist, verbose = FALSE)$snn

	seuratObj <- FindClusters(seuratObj, resolution = 2.0, graph.name = "adt_snn.dist", verbose = FALSE)
	seuratObj[["AdtClusterNames_2.0.dist"]] <- Idents(object = seuratObj)

	#tSNE:
	# Now, we rerun tSNE using our distance matrix defined only on ADT (protein) levels.
	print("Performing tSNE on ADT")
	seuratObj[["tsne_adt.dist"]] <- RunTSNE(adt.dist, assay = assayName, reduction.name = 'adt.tsne.dist', reduction.key = "adtTSNEDist_")

	if (print.plots) {
		print(DimPlot(seuratObj, reduction = "tsne_adt.dist"))
	}

	#UMAP:
	# Now, we rerun UMAP using our distance matrix defined only on ADT (protein) levels.
	if (doUMAP) {
		print("Performing UMAP on ADT")
		seuratObj[["umap_adt.dist"]] <- RunUMAP(adt.dist, assay = assayName, reduction.name = 'adt.umap.dist', reduction.key = "adtUMAPDist_", verbose = FALSE)

		if (print.plots) {
			print(DimPlot(seuratObj, reduction = "umap_adt.dist"))
		}
	}

	#Restore original state:
	Idents(seuratObj) <- seuratObj[["origClusterID", drop = TRUE]]
	DefaultAssay(seuratObj) <- origAssay
	seuratObj[["origClusterID"]] <- NULL

	reductions <- c('tsne')
	if (doUMAP) {
		reductions <- c(reductions, 'umap')
	}

	if (print.plots) {
		for (reduction in reductions) {
			#Compare new/old:
			orig <- DimPlot(seuratObj, reduction = reduction, group.by = "ident", combine = FALSE)[[1]] + NoLegend()
			orig <- orig  +
				labs(title = 'Clustering based on RNA')  +
				theme(plot.title = element_text(hjust = 0.5))

			adt <- DimPlot(seuratObj, reduction = paste0(reduction, "_adt.dist"), group.by = 'AdtClusterNames_2.0.dist', pt.size = 0.5, combine = FALSE)[[1]] + NoLegend()
			adt <- adt  +
				labs(title = 'Clustering based on ADT signal') +
				theme(plot.title = element_text(hjust = 0.5))

			# Note: for this comparison, both the RNA and protein clustering are visualized using the ADT distance matrix.
			print(patchwork::wrap_plots(list(orig, adt), ncol = 2))
		}
	}

	return(seuratObj)
}

#' @title Perform DimRedux on CiteSeq using PCA input
#'
#' @description Perform DimRedux on CiteSeq
#' @param seuratObj The seurat object
#' @param assayName The name of the assay holding the ADT data.
#' @param performClrNormalization If true, Seurat's CLR normalization will be performed. Otherwise this expected data to be pre-normalized
#' @param print.plots If true, QC plots will be printed
#' @param doUMAP If true, RunUMAP will be performed
#' @param adtWhitelist An optional vector of ADTs that will be included in PCA. If NULL, all will be used.
#' @param adtBlacklist An optional vector of ADTs that will be excluded from PCA. If NULL, all will be used.
#' @param dimsToUse An optional vector of PCs to use in tSNE/UMAP
#' @export
#' @importFrom dplyr arrange
#' @import Seurat
CiteSeqDimRedux.PCA <- function(seuratObj, assayName = 'ADT', print.plots = TRUE, performClrNormalization = TRUE, doUMAP = TRUE, dimsToUse = NULL, adtWhitelist = NULL, adtBlacklist = NULL){
	origAssay <- DefaultAssay(seuratObj)

	# Before we recluster the data on ADT levels, we'll stash the original cluster IDs for later
	seuratObj[["origClusterID"]] <- Idents(seuratObj)
	tryCatch({
		DefaultAssay(seuratObj) <- assayName
		print(paste0('Processing ADT data, features: ', paste0(rownames(seuratObj), collapse = ',')))


		if (!is.null(adtWhitelist)) {
			sharedADTs <- intersect(adtWhitelist, rownames(seuratObj[[assayName]]))
			if (length(sharedADTs) != length(adtWhitelist)) {
				missing <- adtWhitelist[!(adtWhitelist %in% rownames(seuratObj[[assayName]]))]
				stop(paste0('The following ADTs were requested in adtWhitelist, but not found in the seurat object: ', paste0(missing, collapse = ',')))
			}
		}

		if (is.null(adtWhitelist)) {
			adtsForPca <- rownames(seuratObj[[assayName]])
		} else {
			adtsForPca <- adtWhitelist
		}

		if (!is.null(adtBlacklist)) {
			adtsForPca <- adtsForPca[!(adtsForPca %in% adtBlacklist)]
		}

		if (performClrNormalization) {
			seuratObj <- NormalizeData(seuratObj, assay = assayName, normalization.method = 'CLR', margin = 2, verbose = FALSE)
		} else if (is.null(Seurat::GetAssayData(seuratObj, assay = assayName, layer = 'data')) || length(Seurat::GetAssayData(seuratObj, assay = assayName, layer = 'data')) == 0){
			stop('Cannot use performClrNormalization=FALSE without pre-existing ADT normalization')
		} else {
			print('Using pre-existing normalization')
		}

		print("Performing PCA on ADTs")
		keyName <- paste0('pca.', tolower(assayName))
		print(paste0('ADTs used: ', paste0(adtsForPca, collapse = ',')))
		seuratObj <- ScaleData(seuratObj, verbose = FALSE, assay = assayName, features = adtsForPca)
		seuratObj <- RunPCA(seuratObj, reduction.name = keyName, assay = assayName, verbose = FALSE, reduction.key  = 'adtPCA_', npcs = length(adtsForPca)-1, features = adtsForPca)
		if (print.plots) {
			print(DimPlot(seuratObj, reduction = keyName))
			print(ElbowPlot(object = seuratObj, reduction = keyName))
		}

		if (is.null(dimsToUse)) {
			dimsToUse <- 1:ncol(seuratObj@reductions[[keyName]])
			print(paste0('Using dims 1:', max(dimsToUse)))
		}

		seuratObj <- FindNeighbors(seuratObj, verbose = FALSE, reduction = keyName, dims = dimsToUse, graph.name = "adt_snn.pca")
		seuratObj <- FindClusters(seuratObj, resolution = 2.0, graph.name = "adt_snn.pca", verbose = FALSE)
		seuratObj[["AdtClusterNames_2.0.PCA"]] <- Idents(object = seuratObj)

		#tSNE:
		print("Performing tSNE on ADT with PCA")
		seuratObj <- RunTSNE(seuratObj, assay = assayName, reduction = keyName, dims = dimsToUse, reduction.name = 'adt.tsne.pca', reduction.key = "adtTSNEPCA_", check_duplicates = FALSE)

		if (print.plots) {
			print(DimPlot(seuratObj, reduction = "adt.tsne.pca"))
		}

		#UMAP:
		# Now, we rerun UMAP using our distance matrix defined only on ADT (protein) levels.
		if (doUMAP) {
			print("Performing UMAP on ADT with PCA")
			seuratObj <- RunUMAP(seuratObj, assay = assayName, reduction = keyName, dims = dimsToUse, reduction.name = 'adt.umap.pca', reduction.key = "adtUMAPPCA_", verbose = FALSE)

			if (print.plots) {
				print(DimPlot(seuratObj, reduction = "adt.umap.pca"))
			}
		}

		reductions <- c('tsne')
		if (doUMAP) {
			reductions <- c(reductions, 'umap')
		}

		if (print.plots) {
			for (reduction in reductions) {
				#Compare new/old:
				orig <- DimPlot(seuratObj, reduction = reduction, group.by = "ident", combine = FALSE)[[1]] + NoLegend()
				orig <- orig  +
					labs(title = 'Clustering based on RNA')  +
					theme(plot.title = element_text(hjust = 0.5))

				adt <- DimPlot(seuratObj, reduction = paste0("adt.", reduction, ".pca"), group.by = 'AdtClusterNames_2.0.PCA', pt.size = 0.5, combine = FALSE)[[1]] + NoLegend()
				adt <- adt  +
					labs(title = 'Clustering based on ADT signal') +
					theme(plot.title = element_text(hjust = 0.5))

				# Note: for this comparison, both the RNA and protein clustering are visualized using the ADT distance matrix.
				print(patchwork::wrap_plots(list(orig, adt), ncol = 2))
			}
		}
	}, error = function(e){
		print('Error running CITE-seq PCA')
		print(conditionMessage(e))
		traceback()
	})

	#Restore original state:
	Idents(seuratObj) <- seuratObj[["origClusterID", drop = TRUE]]
	DefaultAssay(seuratObj) <- origAssay
	seuratObj[["origClusterID"]] <- NULL

	return(seuratObj)
}


#' @title Perform Seurat WNN on CiteSeq
#'
#' @description Perform Seurat WNN on CiteSeq
#' @param seuratObj The seurat object
#' @param dims.list Passed directly to Seurat::FindMultiModalNeighbors
#' @param reduction.list Passed directly to Seurat::FindMultiModalNeighbors
#' @param assayName The assay to use
#' @export
#' @import Seurat
RunSeuratWnn <- function(seuratObj, dims.list = list(1:30, 1:18), assayName = 'ADT', reduction.list = list("pca", paste0("pca.", tolower(assayName)))) {
	if (length(reduction.list) != length(dims.list)) {
		stop('Length of reduction.list and dims.list must be equal')
	}

	i <- 1
	for (reduction in reduction.list) {
		maxDim <- length(seuratObj@reductions[[reduction]])
		argMax <- max(dims.list[[i]])
		if (argMax > maxDim) {
			print(paste0('dims.list requested for: ', reduction, ' is greater than available dims (', maxDim, '), will reduce for FindMultiModalNeighbors'))
			dims.list[[i]] <- 1:maxDim
		}

		i <- i + 1
	}

	seuratObj <- FindMultiModalNeighbors(
		seuratObj, reduction.list = reduction.list,
		dims.list = dims.list, modality.weight.name = "RNA.weight"
	)

	seuratObj <- RunUMAP(seuratObj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", verbose = FALSE)
	seuratObj <- FindClusters(seuratObj, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)

	print(DimPlot(seuratObj, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend())

	return(seuratObj)
}

.PlotMarkerQc <- function(seuratObj, assayName = 'ADT') {
	barcodeMatrix <- Seurat::GetAssayData(seuratObj, layer = 'counts', assay = assayName)
	featuresToPlot <- .GetNonZeroFeatures(seuratObj, assayName)

	setSize <- 2
	steps <- ceiling(length(featuresToPlot) / setSize) - 1

	for (i in 0:steps) {
		start <- (i * setSize) + 1
		end <- min((start + setSize - 1), length(featuresToPlot))
		features <- featuresToPlot[start:end]

		tryCatch({
			suppressMessages(print(RidgePlot(seuratObj, assay = assayName, features = features, ncol = 1)))
		}, error = function(e){
			print(paste0('Error generating ridgeplot: ', paste0(features, collapse = ',')))
			print(conditionMessage(e))
			traceback()
		})
	}

	# Also total per ADT
	countsPerAdt <- rowSums(as.matrix(barcodeMatrix))
	countsPerAdt <- data.frame(Marker = names(countsPerAdt), TotalCount = countsPerAdt)

	P1 <- ggplot(countsPerAdt, aes(x = TotalCount)) +
		geom_density() +
		xlab('Total Count/ADT') +
		ylab('Density') +
		theme(
			axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
		) +
		labs(title = 'Total Counts Per ADT')

	print(P1)

	P2 <- ggplot(countsPerAdt, aes(x = Marker, y = TotalCount)) +
		geom_bar(stat = 'identity') +
		xlab('Marker') +
		ylab('Total Count') +
		theme(
			axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
		) +
		labs(title = 'Total Counts Per ADT')

	print(P2)
}


#' @title Plot Average ADT counts
#'
#' @description Creates a heatmap with average counts of each ADT markers based on one or more grouping fields
#' @param seuratObj The seurat object where data will be added.
#' @param groupFields The directory holding raw count data, generally the raw_feature_bc_matrix from the cellranger outs folder
#' @param assayName The name of the assay holding data
#' @param layer The assay layer to use for average expression data
#' @param normalization.method The normalization method to apply after AverageExpression(). If null, this will be skipped.
#' @param margin This is provided to Seurat::NormalizeData()
#' @param outFile If provided, the heatmap will be written to this file
#' @export
PlotAverageAdtCounts <- function(seuratObj, groupFields = c('ClusterNames_0.2', 'ClusterNames_0.4', 'ClusterNames_0.6'), assayName = 'ADT', layer = 'counts', outFile = NA, normalization.method = 'CLR', margin = 1) {
	for (fn in groupFields) {
		avgSeurat <- Seurat::AverageExpression(seuratObj, return.seurat = T, group.by = fn, assays = assayName, layer = layer)
		if (!is.null(normalization.method)) {
			avgSeurat <- NormalizeData(avgSeurat, normalization.method = normalization.method, margin = margin, verbose = FALSE)
		}

		mat <- t(as.matrix(GetAssayData(avgSeurat, layer = 'data')))
		mat <- mat[,colSums(mat) > 0]

		P1 <- ComplexHeatmap::Heatmap(mat %>% scale_mat(scale = 'column'),
			column_title = paste0('Average ADT Counts By ', fn),
			row_names_side = "left",
			row_dend_side = "right",
			col = Seurat::PurpleAndYellow(20),
			column_names_side = "top",
			column_dend_side = "bottom"
		)

		print(P1)

		if (!is.na(outFile)) {
			P1 <- ggplotify::as.ggplot(grid::grid.grabExpr(ComplexHeatmap::draw(P1)))
			ggsave(P1, filename = outFile)
		}
	}
}

#' @title CalculateUcellPerFeature
#' @description This will iterate all features of the given assay and calculate a per-feature Ucell score
#' @param seuratObj The seurat object
#' @param columnPrefix An optional prefix to be applied to the resulting column (which is otherwise the feature name)
#' @param assayName The name of the assay
#' @param ncores Passed directly to AddModuleScore_UCell
#' @param ignoreCachedRanks If true, any previously cached ranked from UCell will be ignored
#' @export
#' @return A modified Seurat object.
CalculateUcellPerFeature <- function(seuratObj, columnPrefix = NULL, assayName = 'ADT', ncores = 1, ignoreCachedRanks = TRUE) {
	print('Calculating per-feature UCell scores')
	feats <- rownames(seuratObj@assays[[assayName]])
	toCalculate <- list()
	for (feat in feats) {
		cn <- make.names(feat)
		if (!is.null(columnPrefix)) {
			cn <- paste0(columnPrefix, cn)
		}

		toCalculate[[cn]] <- feat
	}

	# Using stored ranks could be a problem if they are based on RNA and we are using a different assay here
	hasCachedRanks <- FALSE
	if (ignoreCachedRanks && 'UCellRanks' %in% names(seuratObj@assays)) {
		print('Storing cached ranks')
		cachedRanks <- seuratObj[['UCellRanks']]
		seuratObj[['UCellRanks']] <- NULL
		hasCachedRanks <- TRUE
	}

	BPPARAM <- .InferBpParam(ncores, defaultValue = NULL)
	seuratObj <- UCell::AddModuleScore_UCell(seuratObj, assay = assayName, features = toCalculate, BPPARAM = BPPARAM, storeRanks = FALSE)

	if (hasCachedRanks) {
		seuratObj[['UCellRanks']] <- cachedRanks
	}

	return(seuratObj)
}