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
#' @param subsetFields If provided, data will be grouped by these fields, and up to targetCells will be retained per group
#' @param seed The random seed
#' @export
DownsampleSeurat <- function(seuratObj, targetCells, subsetFields = NULL, seed = NULL) {
	if (!is.null(seed)) {
		set.seed(seed)
	}

	if (!is.null(subsetFields)){
		for (subsetField in subsetFields) {
			if (!(subsetField %in% names(seuratObj@meta.data))) {
				stop(paste0('Field not found in seuratObj: ', subsetField))
			}
		}
	}

	print(paste0('Subsetting, original cells: ', ncol(seuratObj)))
	cellsToRetain <- c()
	if (is.null(subsetFields)) {
		toSample <- min(targetCells, ncol(seuratObj))
		if (toSample != targetCells) {
			print(paste0('There are only ', ncol(seuratObj), ' cells. Will not downsample.'))
			return(seuratObj)
		}

		cellsToRetain <- sample(colnames(seuratObj), toSample, replace = F)
	} else {
		groupVals <- (seuratObj@meta.data %>% tidyr::unite("x", subsetFields, remove = FALSE))$x
		names(groupVals) <- colnames(seuratObj)
		counts <- table(groupVals)
		print('Unique values: ', paste0(unique(names(counts)), collapse = ','))
		for (val in unique(names(counts))) {
			availBarcodes <- colnames(seuratObj)[groupVals == val]
			toSample <- min(targetCells, length(availBarcodes))
			if (toSample != targetCells) {
				print(paste0('There are only ', length(availBarcodes), ' cells available for group: ', val, '. Will retain all cells for this group.'))
				cellsToRetain <- c(cellsToRetain, availBarcodes)
			} else {
				cellsToRetain <- c(cellsToRetain, sample(availBarcodes, targetCells, replace = F))
			}
		}
	}

	print(paste('Total cells retained ', length(cellsToRetain), ' of ', ncol(seuratObj)))
	seuratObj <- subset(seuratObj, cells = cellsToRetain)

	if (ncol(seuratObj) != length(cellsToRetain)) {
		stop('Incorrect number of cells retained!')
	}

	print(paste0('Final cells: ', ncol(seuratObj)))

	return(seuratObj)
}

#' @title Split Seurat
#'
#' @description Split a seurat object, dividing into new objects based on the value of a field
#' @param seuratObj The seurat object
#' @param splitField The name of the field on which to split the object
#' @param minCellsToKeep If any of the resulting seurat objects have less than this many cells, they will be discarded.
#' @export
SplitSeurat <- function(seuratObj, splitField, minCellsToKeep = 0) {
	if (!(splitField %in% names(seuratObj@meta.data))) {
		stop(paste0('Field not present in seurat object: ', splitField))
	}

	values <- unique(seuratObj@meta.data[[splitField]])

	ret <- list()
	for (value in values) {
		s <- seuratObj[, colnames(seuratObj)[seuratObj@meta.data[[splitField]] == value]]
		if (ncol(s) < minCellsToKeep) {
			print(paste0('Too few cells (', ncol(s), '), discarding subset: ', value))
		} else {
			ret[[value]] <- s
		}
	}

	return (ret)
}

#' @title Subset Seurat Using String Expressions
#'
#' @description Subset a seurat object, selecting cells that match each expression
#' @param seuratObj The seurat object
#' @param expressionStrings A vector strings, each of which will be parsed into an expression, to use in selection
#' @export
SubsetSeurat <- function(seuratObj, expressionStrings = NULL) {
	toUse <- c()
	if (!is.null(expressionStrings)) {
		for (e in expressionStrings) {
			toUse <- c(toUse, parse(text = e))
		}
	}

	if (length(toUse) == 0) {
		stop('Must provide expressionStrings')
	}

	for (expression in toUse) {
		seuratObj <- subset(seuratObj, subset = expression)
	}

	return(seuratObj)
}


#' @title Write Summary Metrics
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
	xynames <- GetXYAesthetics(plot = plot)

	plot.data <- plot$data[cellNames, ]
	names(plot.data)[names(plot.data) == xynames$x] <- 'x'
	names(plot.data)[names(plot.data) == xynames$y] <- 'y'

	return(plot.data)
}

# Set a default value if an object is null
#
# @param lhs An object to set if it's null
# @param rhs The value to provide if x is null
#
# @return rhs if lhs is null, else lhs
#
# @author Hadley Wickham
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
#
`%||%` <- function(lhs, rhs) {
	if (!is.null(x = lhs)) {
		return(lhs)
	} else {
		return(rhs)
	}
}

# Get X and Y aesthetics from a plot for a certain geom
#
# @param plot A ggplot2 object
# @param geom Geom class to filter to
# @param plot.first Use plot-wide X/Y aesthetics before geom-specific aesthetics
# @author Seurat`
# @return A named list with values 'x' for the name of the x aesthetic and 'y' for the y aesthetic
#
GetXYAesthetics <- function(plot, geom = 'GeomPoint', plot.first = TRUE) {
	geoms <- sapply(
		X = plot$layers,
		FUN = function(layer) {
			return(class(x = layer$geom)[1])
		}
	)
	# handle case where raster is set to True
	if (geom == "GeomPoint" && "GeomScattermore" %in% geoms){
		geom <- "GeomScattermore"
	}
	geoms <- which(x = geoms == geom)
	if (length(x = geoms) == 0) {
		stop("Cannot find a geom of class ", geom)
	}
	geoms <- min(geoms)
	if (plot.first) {
		x <- as.character(x = plot$mapping$x %||% plot$layers[[geoms]]$mapping$x)[2]
		y <- as.character(x = plot$mapping$y %||% plot$layers[[geoms]]$mapping$y)[2]
	} else {
		x <- as.character(x = plot$layers[[geoms]]$mapping$x %||% plot$mapping$x)[2]
		y <- as.character(x = plot$layers[[geoms]]$mapping$y %||% plot$mapping$y)[2]
	}
	return(list('x' = x, 'y' = y))
}

#' @title AddClonesToPlot
#' @description Can be used to highlight a set of cells from a seurat plot, such as overlaying specific clonotypes
#' @param seuratObj The seurat object
#' @param plot The plot object, such as the result from FeaturePlot()
#' @param fieldName The name of the field on the seurat object holding cloneName
#' @param colorField If provided, this field will be used
#' @param dotColor An optional string passed to geom_point. Ignored if colorField is provided.
#' @param pt.size The size, passed to geom_point
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


#' @title Calculate Avg Expression
#' @description This is a wrapper around Seurat::AverageExpression, which computes the AverageExpression per group (all assays), then flattens this to a single data frame. It also saves library size and total cells per group.
#' @param seuratObj The seurat object
#' @param groupField The field on which to group. Must be present in seuratObj@meta.data
#' @param slot The slot on which to calculate average expression
#' @export
AvgExpression <- function(seuratObj, groupField, slot = 'counts') {
	if (!(groupField %in% names(seuratObj@meta.data))) {
		stop(paste0('Field not found in seuratObj: ', groupField))
	}

	ret <- Seurat::AverageExpression(seuratObj, assays = NULL, features = rownames(seuratObj), group.by = groupField, slot = slot, verbose = FALSE)

	df <- NULL
	for (assay in names(ret)) {
		df2 <- as.data.frame(ret[[assay]])
		df2$assay <- assay
		df2$feature <- rownames(ret[[assay]])
		df2 <- df2[unique(c('assay', 'feature', colnames(df2)))]
		if (is.null(df)) {
			df <- df2
		} else {
			df <- rbind(df, df2)
		}
	}

	libraryMetrics <- as.data.frame(t(as.matrix(table(seuratObj[[groupField]]))))
	libraryMetrics$feature <- 'TotalCells'
	libraryMetrics$assay <- 'ExperimentMetrics'
	libraryMetrics <- libraryMetrics[unique(c('assay', 'feature', colnames(libraryMetrics)))]
	df <- rbind(df, libraryMetrics)

	for (assay in names(seuratObj@assays)) {
		toAdd <- data.frame(assay = assay, feature = 'LibrarySize')
		for (val in unique(seuratObj@meta.data[[groupField]])) {
			cellsWhitelist <- colnames(seuratObj)[seuratObj@meta.data[[groupField]] == val]

			dat <- Seurat::GetAssay(seuratObj, assay = assay)
			dat <- subset(dat, cells = cellsWhitelist)

			toAdd[val] <- sum(Seurat::GetAssayData(dat, slot = slot))
		}

		df <- rbind(df, toAdd[colnames(df)])
	}

	return(df)
}

#' @title FeaturePlot Across Reductions
#' @description Create FeaturePlots for the selected reductions
#' @param seuratObj The seurat object
#' @param features A vector of features to plot
#' @param reductions The list of reductions to plot
#' @param plotsPerRow The number of plots to print per row
#' @import patchwork
#' @export
FeaturePlotAcrossReductions <- function(seuratObj, features, reductions = c('tsne', 'umap', 'wnn.umap'), plotsPerRow = 3) {
	for (feature in features) {
		reductionToPlot <- intersect(reductions, names(seuratObj@reductions))
		steps <- ceiling(length(reductionToPlot) / plotsPerRow) - 1

		for (i in 0:steps) {
			start1 <- (i * plotsPerRow) + 1
			end <- min((start1 + plotsPerRow - 1), length(reductionToPlot))
			toPlot <- reductionToPlot[start1:end]

			plots <- NULL
			for (reduction in toPlot) {
				P1 <- Seurat::FeaturePlot(seuratObj, reduction = reduction, features = c(feature))
				if (is.null(plots)) {
					plots <- P1
				} else {
					plots <- plots + P1
				}
			}

			plots <- plots + patchwork::plot_layout(ncol = plotsPerRow)

			print(plots + plot_layout(guides = 'collect'))
		}
	}
}

