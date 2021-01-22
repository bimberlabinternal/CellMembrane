



.CheckDuplicatedCellNames <- function(object.list, stop = TRUE){
	cell.names <- unlist(
	x = sapply(
	X = 1:length(x = object.list),
	FUN = function(x) Cells(object.list[[x]])
	)
	)

	dups <- duplicated(x = cell.names)
	if (any(dups)) {
		if (stop){
			stop(paste0('There were duplicated cell names: ',  paste0(head(unique(cell.names[dups])), collapse = ',')))
		} else {
			print('There were duplicated cell names')
			print(head(unique(cell.names[dups])))
		}
	}
}

DownsampleSeurat <- function(seuratObj, targetCells, subsetFields = NULL, seed = NULL) {
	if (!is.null(seed)) {
		set.seed(seed)
	}

	availFeats = colnames(seuratObj@meta.data)
	if (!is.null(subsetFields)){
		if (length(featureSplit) > 1) {
			warning("length of featureSplit > 1 only 1st is used")
			featureSplit = featureSplit[1]
			if (!(featureSplit %in% availFeats)) {
				stop("featureSplit not found in meta.data")
			}
		}
	} else {
		stop("featureSplit is NULL")
	}

	CellCountMat = table(seuratObj@meta.data[,featureSplit]) %>% unlist() %>% as.matrix()
	if (targetCells == 0) {
		targetCells = min(CellCountMat[,1])
		print(paste0("total cells set by min to:", targetCells))
	}

	if(!("barcode" %in% availFeats)) {
		warning("barcode not found in meta.data")
		seuratObj$barcode = paste("cell_", 1:nrow(seuratObj@meta.data))
	}

	splitLevs = levels(factor(seuratObj@meta.data[,featureSplit]))
	barcodeLS = lapply(splitLevs, function(xSL){
		availBarcodes = seuratObj@meta.data[which(seuratObj@meta.data[,featureSplit] == xSL),]$barcode

		if (length(availBarcodes)<targetCells) {
			warning(paste0(xSL, " had less cells than targetCells"))
			availBarcodes
		}else {
			sample(availBarcodes, targetCells, replace = F)
		}
	})

	return(subset(seuratObj, barcode %in% unlist(barcodeLS)))
}

SplitSeurat <- function(seuratObj, subsetField) {
	if (!(subsetField %in% names(seuratObj@meta.data))) {
		stop(paste0('Field not present in seurat object: ', subsetField))
	}

	values <- unique(seuratObj@meta.data[[subsetField]])

	ret <- list()
	for (value in values) {
		expr <- parse(text = paste0(subsetField, " == '", value, "'"))
		ret[value] <- seuratObj[, Seurat::WhiteCells(object = seuratObj, expression = expr)]
	}

	return (ret)
}

#' @title WriteSummaryMetrics
#' @export
#' @param seuratObj, A Seurat object.
#' @param file The file where metrics will be written
WriteSummaryMetrics <- function(seuratObj, file) {
	df <- data.frame(Category = "Seurat", MetricName = "TotalCells", Value = ncol(seuratObj))
	df <- rbind(df, data.frame(Category = "Seurat", MetricName = "TotalFeatures", Value = nrow(seuratObj)))

	if ('HighlyActivated.Call' %in% names(seuratObj@meta.data)) {
		val <- sum(seuratObj$HighlyActivated.Call) / ncol(seuratObj)
		df <- rbind(df, data.frame(Category = "Seurat", MetricName = "FractionActivated", Value = val))
	}

	write.table(df, file = file, quote = F, row.names = F, sep = '\t')
}

#' @title WriteCellBarcodes
#' @description Writes a table of cell barcodes to the provided file
#' @return A modified Seurat object.
#' @param seuratObj The seurat object
#' @param file The output file
#' @param The file to which barcodes will be written
#' @export
WriteCellBarcodes <- function(seuratObj, file) {
	df <- data.frame(CellBarcode = colnames(seuratObj))

	write.table(df, file = file, quote = F, row.names = F, sep = ',', col.names = F)
}

#' @title GetXYDataFromPlot
#' @description Get XY data from a seurat plot
#' @param plot The plot object
#' @param cellNames The set of cells to export
#' @export
GetXYDataFromPlot <- function(plot, cellNames) {
	xynames <- Seurat:::GetXYAesthetics(plot = plot)

	plot.data <- plot$data[cellNames, ]
	names(plot.data)[names(plot.data) == xynames$x] <- 'x'
	names(plot.data)[names(plot.data) == xynames$y] <- 'y'

	return(plot.data)
}

#' @title AddClonesToPlot
#' @description Can be used to highlight a set of cells from a seurat plot, such as overlaying specific clonotypes
#' @param seuratObj The seurat object
#' @param plot The plot object, such as the result from FeaturePlot()
# '@param fieldName The name of the field on the seurat object holding cloneName
# '@param colorField If provided, this field will be used
# '@param dotColor
# '@param pt.size
#' @export
#' @import ggplot2
AddClonesToPlot <- function(seuratObj, plot, fieldName = 'CloneNames', colorField = NA, dotColor = NA, pt.size = 1) {
	cellNames <- colnames(seuratObj)[!is.na(seuratObj[[fieldName]])]
	plot.data <- GetXYDataFromPlot(plot, cellNames)
	plot.data$Clone <- seuratObj[[fieldName]][!is.na(seuratObj[[fieldName]])]

	sel <- !is.na(seuratObj[[fieldName]])
	plot.data$CloneName <- naturalsort::naturalfactor(seuratObj[[fieldName]][sel])
	if (!is.na(colorField)) {
		plot.data$CloneColor <- naturalsort::naturalfactor(seuratObj[[colorField]][sel])

		plot <- plot + geom_point(
		mapping = aes_string(x = 'x', y = 'y', shape = 'CloneName', color = 'CloneColor'),
		data = plot.data,
		size = pt.size,
		inherit.aes = F
		)
	} else {
		plot <- plot + geom_point(
		mapping = aes_string(x = 'x', y = 'y', shape = 'CloneName'),
		data = plot.data,
		size = pt.size,
		inherit.aes = F,
		color = dotColor
		)
	}

	return(plot)
}

#' @title FilterCloneNames
#' @description Filter Clone Names
#' @param seuratObj The seurat object
#' @param minValue Filters clones not present in at least this many cells
#' @export
FilterCloneNames <- function(seuratObj, minValue) {
	ct <- table(seuratObj$CloneNames)
	ct <- ct[ct < minValue]

	seuratObj$CloneNames[seuratObj$CloneNames %in% names(ct)] <- NA

	return(seuratObj)
}
