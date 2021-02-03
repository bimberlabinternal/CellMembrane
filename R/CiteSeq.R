#' @include Utils.R
#' @include Preprocessing.R
#' @import Seurat


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
#' @export
#' @importFrom dplyr arrange
AppendCiteSeq <- function(seuratObj, unfilteredMatrixDir, normalizeMethod = 'dsb', datasetId = NULL, assayName = 'ADT', featureMetadata = NULL, adtWhitelist = NULL, minRowSum = NULL) {
	print(paste0('Initial cell barcodes in GEX data: ', ncol(seuratObj)))
	if (!is.null(datasetId)) {
		gexCells <- colnames(seuratObj)[seuratObj$DatasetId == datasetId]
		print(paste0('Initial cell barcodes in GEX data for prefix: ', length(gexCells)))
	} else {
		gexCells <- colnames(seuratObj)
	}
	assayData <- .LoadCiteSeqData(unfilteredMatrixDir, datasetId = datasetId, featureMetadata = featureMetadata, adtWhitelist = adtWhitelist, minRowSum = minRowSum)

	sharedCells <- colnames(assayData)[which(colnames(assayData) %in% gexCells)]
	print(paste0('Intersect with GEX data: ', length(sharedCells)))

	if (length(sharedCells) == 0) {
		print(head(gexCells))
		print(head(colnames(assayData)))
		stop('No cells are shared')
	}

	#TODO: normalize:
	if (normalizeMethod == 'dsb') {
		assayData <- .NormalizeDsb(seuratObj, assayData, cellWhitelist = sharedCells)
	} else {
		print('No normalization will be performed')
		assayData <- subset(assayData, cells = sharedCells)
	}

	#Note: we need to support merging with any existing data and account for datasetId:
	if (assayName %in% names(seuratObj@assays)) {
		print('Existing assay found, merging counts')
		assayData <- .MergeAdtWithExisting(assayData, GetAssayData(seuratObj, assay = assayName))
	} else {
		print('no pre-existing assay, creating new')
		assayData <- .PossiblyExpandAssay(seuratObj, assayData)
	}

	seuratObj[[assayName]] <- assayData

	# if (!skipNormalize) {
	# 	seuratObj <- NormalizeData(seuratObj, normalization.method = 'CLR', margin = 2, assay = assayName)
	# 	seuratObj <- ScaleData(seuratObj, assay = assayName)
	# }

	return(seuratObj)
}

.PossiblyExpandAssay <- function(seuratObj, assayData) {
	#now add empty cells for those lacking ADTs:
	#Append blank cells for any in GEX but missing in ADT:
	missing <- colnames(seuratObj)[!(colnames(seuratObj) %in% colnames(assayData))]
	if (length(missing) > 0) {
		missingMat <- matrix(rep(0, nrow(assayData) * length(missing)), nrow = nrow(assayData), ncol = length(missing))
		colnames(missingMat) <- missing
		print(paste0('total cells lacking data: ', length(missing)))

		replacementAssay <- NULL
		for (slot in c('counts', 'data')) {
			newData <- GetAssayData(assayData, slot = slot)
			if (is.null(newData)) {
				next
			}

			extendedData <- cbind(newData, missingMat)
			extendedData <- as.sparse(extendedData[,colnames(seuratObj), drop = F])

			if (is.null(replacementAssay)) {
				args <- list()
				args[[slot]] <- extendedData
				replacementAssay <- do.call(Seurat::CreateAssayObject, args)
			} else {
				replacementAssay <- SetAssayData(replacementAssay, slot = slot, extendedData)
			}
		}

		assayData <- replacementAssay
	}

	if (sum(colnames(seuratObj) != colnames(assayData)) > 0) {
		stop('The columns of the ADT matrix do not equal the seurat object')
	}

	return(assayData)
}

.MergeAdtWithExisting <- function(newAssay, existingAssay) {
	if (sum(colnames(seuratObj) != colnames(newAssay)) > 0) {
		stop('The columns of the ADT matrix do not equal the seurat object')
	}

	allFeatures <- unique(rownames(existingAssay), rownames(newAssay))
	newFeatures <- rownames(existingAssay)[!(rownames(existingAssay) %in% allFeatures)]

	replacementAssay <- NULL
	for (slot in c('counts', 'data')) {
		newData <- GetAssayData(newAssay, slot = slot)
		if (is.null(newData)) {
			next
		}

		existingData <- GetAssayData(existingAssay, slot = slot)
		if (is.null(existingData)) {
			existingData <- Seurat::as.sparse(matrix(rep(0, ncol(existingAssay) * nrow(existingAssay),  ncol = ncol(existingAssay), nrow = length(existingAssay))))
		}

		# Add any new ADTs from this dataset, if needed:
		if (length(newFeatures) > 0) {
			missingMat <- matrix(rep(0, ncol(existingAssay) * length(newFeatures)), ncol = ncol(assayData), nrow = length(newFeatures))
			rownames(missingMat) <- newFeatures
			print(paste0('total ADT rows added to assay: ', length(newFeatures)))
			existingData <- Seurat::as.sparse(rbind(existingData, missingMat))
		}

		existingData[rownames(newData),colnames(newData), drop = F] <- newData

		if (is.null(replacementAssay)) {
			args <- list()
			args[[slot]] <- existingData
			replacementAssay <- do.call(Seurat::CreateAssayObject, args)
		} else {
			replacementAssay <- SetAssayData(replacementAssay, slot = slot, existingData)
		}
	}

	existingAssay <- replacementAssay

	return(existingAssay)
}

.LoadCiteSeqData <- function(unfilteredMatrixDir, datasetId = NULL, featureMetadata = NULL, adtWhitelist = NULL, minRowSum = NULL) {
	if (!dir.exists(unfilteredMatrixDir)){
		stop("Count matrix not found")
	}

	bData <- Seurat::Read10X(unfilteredMatrixDir, gene.column=1, strip.suffix = TRUE)
	bData <- bData[which(!(rownames(bData) %in% c('unmapped'))), , drop = F]
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
			bData <- bData[!toDrop, , drop = FALSE]
			print(paste0('ADTs after filter: ', nrow(bData)))
			if (nrow(bData) == 0) {
				stop('No ADTs remain after filter')
			}
		}
	}

	# Apply this prior to rename:
	if (!is.null(adtWhitelist)) {
		print('Filtering ADTs based on whitelist')
		bData <- bData[rownames(bData) %in% adtWhitelist, , drop = FALSE]
		print(paste0('ADTs after filter: ', nrow(bData)))
		if (nrow(bData) == 0) {
			stop('No rows left after filter!')
		}
	}

	bData <- Seurat::CreateAssayObject(counts = bData)

	# Optionally rename features
	if (!all(is.null(featureMetadata))) {
		if (!('tagname' %in% colnames(featureMetadata))) {
			stop('featureMetadata must contain the column tagname')
		}

		print('Renaming ADTs')
		rownames(featureMetadata) <- featureMetadata$featureName
		featureMetadata$tagname <- paste0(featureMetadata$tagname, '-', featureMetadata$sequence)
		newRows <- data.frame(tagname = rownames(bData), sortorder = 1:nrow(bData), stringsAsFactors = F)
		newRows <- merge(newRows , featureMetadata, by = 'tagname', all.x = T, all.y = F)
		newRows <- newRows %>% arrange(sortorder)
		newRows <- newRows[names(newRows) != 'sortorder']

		newRows$markername <- dplyr::coalesce(newRows$markername, newRows$tagname)
		d <- duplicated(newRows$markername)
		if (sum(d) > 0) {
			stop('There were duplicate marker names after rename: ' + paste0(newRows$markername[d], collapse = ', '))
		}

		print(paste0('Total renamed: ', sum(newRows$markername != newRows$tagname)))

		rownames(bData) <- newRows$markername
		featureMetadata <- newRows

		bData <- Seurat::AddMetaData(object = bData, metadata = featureMetadata)
	}

	return(bData)
}

.AddFeatureMetadata <- function(ft, assayData) {
	print('Adding feature metadata')
	for (colname in names(ft)) {
		print(paste0('adding column: ', colname))
		toAdd <- ft[[colname]]
		names(toAdd) <- ft$markernameUnique
		select <- names(toAdd) %in% rownames(assayData)
		toAdd <- toAdd[select]
		if (length(toAdd) != nrow(ft)) {
			print(paste0('Not all markers in feature table were present in found table, missing: ', paste0(ft$markernameUnique[!select])))
		}

		if (length(toAdd) > 0) {
			assayData <- Seurat::AddMetaData(object = assayData, metadata = toAdd, col.name = colname)
		}
	}

	return(assayData)
}


.PlotCiteSeqCountData <- function(seuratObj, assayName = 'ADT') {
	assayData <- GetAssayData(seuratObj, slot = "counts", assay = assayName)

	featuresToPlot <- rownames(assayData)
	setSize <- 2
	steps <- ceiling(length(featuresToPlot) / setSize) - 1

	for (i in 0:steps) {
		start <- (i * setSize) + 1
		end <- min((start + setSize - 1), length(featuresToPlot))
		features <- featuresToPlot[start:end]

		print(RidgePlot(seuratObj, assay = assayName, features = features, ncol = 1))
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
.NormalizeDsb <- function(seuratObj, unfilteredAdtAssay, rnaAssayName = 'RNA', prot.size.min = 1.4, prot.size.max = 2.5, cellWhitelist = sharedCells) {
	gexCountMatrix <- GetAssayData(object = seuratObj, slot = 'counts', assay = rnaAssayName)
	unfilteredAdtCountMatrix <- GetAssayData(object = unfilteredAdtAssay, slot = 'counts')

	# create a metadata dataframe of simple qc stats for each droplet
	rna_size <- log10(Matrix::colSums(gexCountMatrix))
	prot_size <- log10(Matrix::colSums(unfilteredAdtCountMatrix))
	ngene <- Matrix::colSums(gexCountMatrix > 0)

	md <- as.data.frame(cbind(rna_size, ngene, prot_size))
	md$bc <- rownames(md)

	p1 <- ggplot(md[md$rna_size > 0, ], aes(x = rna_size)) + geom_density(fill = "dodgerblue") + ggtitle("RNA library size \n distribution")
	p2 <- ggplot(md[md$prot_size> 0, ], aes(x = prot_size)) + geom_density(fill = "firebrick2") + ggtitle("Protein library size \n distribution")
	print(p1 + p2)

	# define a vector of background / empty droplet barcodes based on protein library size and mRNA content
	background_drops <- md$bc[md$prot_size > prot.size.min & md$prot_size < prot.size.max]
	background_drops <- background_drops[!(background_drops %in% colnames(seuratObj))]
	negative_mtx_rawprot <- unfilteredAdtCountMatrix[ , background_drops] %>% as.matrix()

	# define a vector of cell-containing droplet barcodes based on protein library size and mRNA content
	filteredAdtCountMatrix <- unfilteredAdtCountMatrix[ , cellWhitelist]
	hasProteinCounts <- rownames(md)[md$prot_size > 2.8]
	filteredAdtCountMatrix <- unfilteredAdtCountMatrix[ , filteredAdtCountMatrix %in% hasProteinCounts] %>% as.matrix()

	print(5)

	normCounts <- dsb::DSBNormalizeProtein(
		cell_protein_matrix = filteredAdtCountMatrix, # cell containing droplets
		empty_drop_matrix = negative_mtx_rawprot, # estimate ambient noise with the background drops
		denoise.counts = TRUE, # model and remove each cell's technical component
		use.isotype.control = FALSE # use isotype controls to define the technical component
	)

	print(6)
	return(Seurat::CreateAssayObject(counts = filteredAdtCountMatrix[,colnames(normCounts)], data = normCounts))
}


CiteSeqDimRedux <- function(seuratObj, assayName = 'ADT', dist.method="euclidean", print.plots = T, rnaAssayName = 'RNA'){
	origAssay <- DefaultAssay(seuratObj)
	DefaultAssay(seuratObj) <- assayName
	print(paste0('Processing ADT data, features: ', paste0(rownames(seuratObj), collapse = ',')))

	# Before we recluster the data on ADT levels, we'll stash the original cluster IDs for later
	seuratObj[["origClusterID"]] <- Idents(seuratObj)

	#PCA:
	print("Performing PCA on ADT")
	seuratObj <- RunPCA(seuratObj, features = rownames(seuratObj), reduction.name = "pca_adt", reduction.key = "pcaadt_", verbose = FALSE)
	if (print.plots) {
		print(DimPlot(seuratObj, reduction = "pca_adt"))
	}

	#SNN:
	print("Calculating Distance Matrix")
	adt.data <- GetAssayData(seuratObj, slot = "data")
	adt.dist <- dist(t(adt.data), method = dist.method)
	seuratObj[["adt_snn"]]  <- FindNeighbors(adt.dist)$snn

	#Cluster with a few different resolutions
	seuratObj <- FindClusters(seuratObj, resolution = 0.1, graph.name = "adt_snn")
	seuratObj <- FindClusters(seuratObj, resolution = 0.2, graph.name = "adt_snn")
	seuratObj <- FindClusters(seuratObj, resolution = 0.5, graph.name = "adt_snn")
	seuratObj <- FindClusters(seuratObj, resolution = 1.0, graph.name = "adt_snn")
	
	#tSNE:
	# Now, we rerun tSNE using our distance matrix defined only on ADT (protein) levels.
	print("Performing tSNE on ADT")
	seuratObj[["tsne_adt"]] <- RunTSNE(adt.dist, assay = assayName, reduction.key = "adtTSNE_")

	if (print.plots) {
		print(DimPlot(seuratObj, reduction = "tsne_adt"))
	}

	#UMAP:
	# Now, we rerun UMAP using our distance matrix defined only on ADT (protein) levels.
	print("Performing UMAP on ADT")
	seuratObj[["umap_adt"]] <- RunUMAP(adt.dist, assay = assayName, reduction.key = "adtUMAP_")

	if (print.plots) {
		print(DimPlot(seuratObj, reduction = "umap_adt"))
	}

	#Restore original state:
	Idents(seuratObj) <- seuratObj[["origClusterID"]]
	DefaultAssay(seuratObj) <- origAssay
	seuratObj[["origClusterID"]] <- NULL

	if (print.plots) {
		#Compare new/old:
		tsne_orig <- DimPlot(seuratObj, reduction = "tsne", group.by = "origClusterID", combine = FALSE)[[1]] + NoLegend()
		tsne_orig <- tsne_orig  +
			labs(title = 'Clustering based on scRNA-seq')  +
			theme(plot.title = element_text(hjust = 0.5))
		tsne_orig <- LabelClusters(plot = tsne_orig, id = "origClusterID", size = 6)

		tsne_adt <- DimPlot(seuratObj, reduction = "tsne_adt", pt.size = 0.5, combine = FALSE)[[1]] + NoLegend()
		tsne_adt <- tsne_adt  +
			labs(title = 'Clustering based on ADT signal') + theme(plot.title = element_text(hjust = 0.5))
		tsne_adt <- LabelClusters(plot = tsne_adt, id = "ident", size = 6)

		# Note: for this comparison, both the RNA and protein clustering are visualized on a tSNE generated using the ADT distance matrix.
		print(patchwork::wrap_plots(list(tsne_orig, tsne_adt), ncol = 2))
	}

	# WNN:
	if (!is.null(rnaAssayName)) {

	}

	return(seuratObj)
}