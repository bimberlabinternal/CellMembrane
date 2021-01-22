utils::globalVariables(
	names = c('subsetField'),
	package = 'CellMembrane',
	add = TRUE
)

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


#' @title Downsample Seurat
#'
#' @description Downsample a seurat object, either globally or subset by a field
#' @param seuratObj The seurat object
#' @param targetCells The desired cell number to retain per unit of data. If a subsetField is provided, the string 'min' can also be used, in which case
#' @param subsetField If provided, data will be grouped by this field, and up to targetCells will be retained per group
#' @param seed The random seed
#' @export
DownsampleSeurat <- function(seuratObj, targetCells, subsetField = NULL, seed = NULL) {
	if (!is.null(seed)) {
		set.seed(seed)
	}

	if (!is.null(subsetField) && !(subsetField %in% rownames(seuratObj@meta.data))){
			stop(paste0('Field not found in seuratObj: ', subsetField))
	}

	cellsToRetain <- c()
	if (is.null(subsetField)) {
		cellsToRetain <- sample(colnames(seuratObj), targetCells, replace = F)
	} else {
		counts <- table(seuratObj[[subsetField]])
		for (val in unique(names(counts))) {
			availBarcodes <- rownames(seuratObj)[seuratObj[[subsetField]] == val]
			cellsToRetain <- c(cellsToRetain, sample(availBarcodes, targetCells, replace = F))
		}
	}

	return(subset(seuratObj, cells = cellsToRetain))
}

#' @title Split Seurat
#'
#' @description Split a seurat object, dividing into new objects based on the value of a field
#' @param seuratObj The seurat object
#' @param splitField The name of the field on which to split the object
#' @export
SplitSeurat <- function(seuratObj, splitField) {
	if (!(subsetField %in% names(seuratObj@meta.data))) {
		stop(paste0('Field not present in seurat object: ', splitField))
	}

	values <- unique(seuratObj@meta.data[[splitField]])

	ret <- list()
	for (value in values) {
		expr <- parse(text = paste0(splitField, " == '", value, "'"))
		ret[value] <- seuratObj[, Seurat::WhichCells(object = seuratObj, expression = expr)]
	}

	return (ret)
}

#' @title Subset Seurat
#'
#' @description Subset a seurat object, selecting cells that match each expression
#' @param seuratObj The seurat object
#' @param expressions A vector of expressions to use in selection
#' @export
SubsetSeurat <- function(seuratObj, expressions) {
	for (expression in expressions) {
		seuratObj <- subset(seuratObj, subset = expression)
	}

	return(seuratObj)
}


#' @title WriteSummaryMetrics
#' @export
#' @param seuratObj A Seurat object.
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

#' @title Write Cell Barcodes
#' @description Writes a table of cell barcodes to the provided file
#' @param seuratObj The seurat object
#' @param file The file to which barcodes will be written
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
