#' @importFrom rlang :=

utils::globalVariables(
	names = c('subsetField', 'CountsPerCell', 'Saturation', 'Color', 'Label', 'clone_vector'),
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
DownsampleSeurat <- function(seuratObj, targetCells, subsetFields = NULL, seed = GetSeed()) {
	if (!is.null(seed)) {
		set.seed(seed)
	}

	if (!is.null(subsetFields)){
		for (subsetField in subsetFields) {
			if (!(subsetField %in% names(seuratObj@meta.data))) {
				stop(paste0('Field not found in seuratObj: [', subsetField, ']'))
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
		groupVals <- (seuratObj@meta.data %>% tidyr::unite("x", tidyr::all_of(subsetFields), remove = FALSE))$x
		names(groupVals) <- colnames(seuratObj)
		counts <- table(groupVals)
		print(paste0('Unique values: ', paste0(unique(names(counts)), collapse = ',')))
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
		stop(paste('Incorrect number of cells retained!, was: ', ncol(seuratObj), ', expected: ', length(cellsToRetain)))
	}

	print(paste0('Final cells: ', ncol(seuratObj)))

	return(seuratObj)
}

#' @title Split Seurat
#'
#' @description Split a seurat object, dividing into new objects based on the value of a field
#' @param seuratObj The seurat object
#' @param splitField The name of the field on which to split the object
#' @param minCellsToKeep If any of the resulting seurat objects have less than this many cells, they will be discarded. If this value is less than 1, it will be interpreted as a fraction of the total input cells.
#' @param naOtherLabel This string will be used to label any cells marked NA.
#' @param excludedClasses Any cells with these labels will be lumped into the NA/Other bin.
#' @param appendLowFreqToOther If true, any cells with NAs for the splitField, or terms with fewer than minCellsToKeep, will be merged into a single seurat object
#' @param alwaysRetainOtherClass If true, even if the number of cells is less than minCellsToKeep, this class with be retained.
#' @export
SplitSeurat <- function(seuratObj, splitField, minCellsToKeep = 0.02, naOtherLabel = 'Other', excludedClasses = NULL, appendLowFreqToOther = TRUE, alwaysRetainOtherClass = FALSE) {
  if (!(splitField %in% names(seuratObj@meta.data))) {
		stop(paste0('Field not present in seurat object: ', splitField))
	}

	if (minCellsToKeep > 0 && minCellsToKeep < 1) {
		minCellsToKeepOrig <- minCellsToKeep
		minCellsToKeep <- ncol(seuratObj) * minCellsToKeep
		print(paste0('Interpreting minCellsToKeep as a fraction of input cells. Converting from ', minCellsToKeepOrig, ' to: ', minCellsToKeep))
	}

	#Fix NAs in splitField
	data <- as.character(seuratObj@meta.data[[splitField]])
	data[is.na(data)] <- naOtherLabel
	if (!all(is.null(excludedClasses))) {
		toDrop <- unique(data %in% excludedClasses)
		if (length(toDrop) > 0) {
			print(paste0('Merging the following classes into ', naOtherLabel, ': ', paste0(toDrop, collapse = ',')))
			data[data %in% excludedClasses] <- naOtherLabel
		}
	}
	values <- unique(data)
	cellsForOther <- colnames(seuratObj)[data == naOtherLabel]

	ret <- list()
	for (value in values) {
		# Handle this below:
		if (value == naOtherLabel) {
			next
		}

		toKeep <- colnames(seuratObj)[data == value]
		if (length(toKeep) == 0 || length(toKeep) < minCellsToKeep) {
			print(paste0('Too few cells (', length(toKeep), '): ', value, ', ', ifelse(appendLowFreqToOther, yes = paste0('cells will be merged into: ', naOtherLabel), no = 'discarding cells')))
			if (appendLowFreqToOther) {
				cellsForOther <- c(cellsForOther, toKeep)
			}
		} else {
			s <- subset(seuratObj, cells = toKeep)
			ret[[as.character(value)]] <- s
		}
	}

	if (length(cellsForOther) > 0) {
		if (!alwaysRetainOtherClass && length(cellsForOther) < minCellsToKeep) {
			print(paste0('A total of ', length(cellsForOther), ' cells are in low-frequency groups, which is below minCellsToKeep and will be dropped'))
		} else {
			print(paste0('Adding category for other (', length(cellsForOther), ' cells):'))
			s <- subset(seuratObj, cells = cellsForOther)
			ret[[as.character(naOtherLabel)]] <- s
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

#' @title HighlightCellsOnSeuratPlot
#' @description Can be used to highlight a set of cells from a seurat plot, such as overlaying specific clonotypes. Note: this uses ggnewscale::new_scale_color to enable multiple color scales. The original aesthetics are suffixed with '_new' (i.e. color_new). This has also only been tested thoroughly with DimPlots.
#' @param seuratObj The seurat object
#' @param seuratPlot The plot object, such as the result from DimPlot()
#' @param cellSelectField The name of the field to select which cells to plot
#' @param colorLegendLabel This is passed to labs(color = XXX) to label the legend
#' @param colorField The name of the field to assign colors to the cells. This provides the option to use a different value from cellSelectField. If colorField and dotColor are NA, cellSelectField will be used.
#' @param dotColor An optional string passed to geom_point(color = XX). This will assign all cells the same color. Ignored if colorField is provided.
#' @param dotShapes An optional vector of shape values passed to scale_shape_manual(values = XX). If dotShapes is provided, but shapeField is not, all cells will receive the same shape.
#' @param pt.size The size, passed to geom_point().
#' @param shapeField If true, provided, these values will be used for shape. Note: if there are more than 6 unique values this will be ignored.
#' @param horizontalLegend If true, theme(legend.box = "horizontal") is added to the plot. This can be useful if the original plot also has a legend (such as cluster names)
#' @param resetLegendSize Using ggnewscale::new_scale_color seems to reset the dot size of the legend. If resetLegendSize=TRUE, the function will call "guides(colour_new = guide_legend(override.aes = list(size=3)))" to restore dot size to 3.
#' @param maxAllowableShapeValues If shapes are used and there are more than this many values, all cells will receive the same shape.
#' @export
#' @import ggplot2
HighlightCellsOnSeuratPlot <- function(seuratObj, seuratPlot, cellSelectField = 'CloneNames', colorLegendLabel = 'Clone', colorField = NA, dotColor = NA, pt.size = 1, shapeField = NA, dotShapes = NA, horizontalLegend = TRUE, resetLegendSize = TRUE, maxAllowableShapeValues = 6) {
	cellNames <- colnames(seuratObj)[!is.na(seuratObj[[cellSelectField, drop = TRUE]])]
	plot.data <- GetXYDataFromPlot(seuratPlot, cellNames)
	plot.data$Clone <- seuratObj[[cellSelectField, drop = TRUE]][!is.na(seuratObj[[cellSelectField, drop = TRUE]])]
	cellsWithData <- !is.na(seuratObj[[cellSelectField, drop = TRUE]])

	seuratPlot <- seuratPlot + ggnewscale::new_scale_color()

	if (is.na(dotColor) && is.na(colorField)) {
		colorField <- cellSelectField
	}

	if (!is.na(colorField)) {
		plot.data$Color <- naturalsort::naturalfactor(seuratObj[[colorField, drop = TRUE]][cellsWithData])

		seuratPlot <- seuratPlot + geom_point(
			mapping = aes(x = x, y = y, color = Color),
			data = plot.data,
			size = pt.size,
			inherit.aes = F
		)
	} else if (!is.na(dotColor)) {
		seuratPlot <- seuratPlot + geom_point(
			mapping = aes(x = x, y = y),
			data = plot.data,
			size = pt.size,
			inherit.aes = F,
			color = dotColor
		)
	} else {
		stop('Must provide either dotColor or colorField')
	}

	# If shape field is not provided, but dotShape is, add a dummy field and assume all are the same value
	if (!is.na(shapeField)) {
		plot.data$ShapeField <- as.factor(seuratObj[[shapeField, drop = TRUE]][cellsWithData])
	} else if (!is.na(dotShapes)) {
		plot.data$ShapeField <- 1
	}

	if ('ShapeField' %in% names(plot.data)) {
		if (length(unique(plot.data$ShapeField)) > maxAllowableShapeValues) {
			warning(paste0('There are more than ', maxAllowableShapeValues, ' unique values, all points will be assigned the same shape. See maxAllowableShapeValues'))
		} else {
			seuratPlot <- seuratPlot + aes(shape = 'ShapeField')
			seuratPlot <- seuratPlot + guides(shape = FALSE)

			if (!is.na(dotShapes)) {
				seuratPlot <- seuratPlot + scale_shape_manual(values = dotShapes)
			}
		}
	}

	if (!is.na(colorLegendLabel)) {
		seuratPlot <- seuratPlot + labs(color = colorLegendLabel)
	}

	if (horizontalLegend) {
		seuratPlot <- seuratPlot + theme(legend.box = 'horizontal')
	}

	# This restores the original DimPlot dot size
	if (resetLegendSize) {
		seuratPlot <- seuratPlot + guides(colour_new = guide_legend(override.aes = list(size=3)))
	}

	return(seuratPlot)
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

	libraryMetrics <- as.data.frame(t(as.matrix(table(seuratObj[[groupField, drop = TRUE]]))))
	libraryMetrics$feature <- 'TotalCells'
	libraryMetrics$assay <- 'ExperimentMetrics'
	libraryMetrics <- .CheckColnamesAreNumeric(libraryMetrics)
	libraryMetrics <- libraryMetrics[unique(c('assay', 'feature', colnames(libraryMetrics)))]
	df <- rbind(df, libraryMetrics)

	for (assay in names(seuratObj@assays)) {
		toAdd <- data.frame(assay = assay, feature = 'LibrarySize')
		for (val in unique(seuratObj@meta.data[[groupField]])) {
			cellsWhitelist <- colnames(seuratObj)[seuratObj@meta.data[[groupField]] == val]

			dat <- Seurat::GetAssay(seuratObj, assay = assay)
			dat <- subset(dat, cells = cellsWhitelist)
			if (ncol(dat) != length(cellsWhitelist)) {
				stop(paste0('Incorrect assay subset. Expected: ', length(cellsWhitelist), ', actual: ', ncol(dat)))
			}

			toAdd[val] <- sum(Seurat::GetAssayData(dat, slot = slot))
		}

		toAdd <- .CheckColnamesAreNumeric(toAdd)
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
	reductionToPlot <- intersect(reductions, names(seuratObj@reductions))
	if (length(reductionToPlot) == 0) {
		print('None of the requested reductions were present, skipping')
	}

	for (feature in features) {
		dat <- Seurat::FetchData(seuratObj, vars = feature, slot = 'data')
		if (all(is.na(dat)) || max(dat, na.rm = T) == 0) {
			print(paste0('Skipping feature with zero/NA counts: ', feature))
		}

		steps <- ceiling(length(reductionToPlot) / plotsPerRow) - 1

		for (i in 0:steps) {
			start1 <- (i * plotsPerRow) + 1
			end <- min((start1 + plotsPerRow - 1), length(reductionToPlot))
			toPlot <- reductionToPlot[start1:end]

			plots <- NULL
			for (reduction in toPlot) {
				P1 <- Seurat::FeaturePlot(seuratObj, reduction = reduction, features = feature, min.cutoff = 'q05', max.cutoff = 'q95')
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

.CalcPerCellSaturation <- function(seuratObj, molInfoFile, cellbarcodePrefix = NULL) {
	df <- DropletUtils::get10xMolInfoStats(molInfoFile)
	df$cellbarcode <- df$cell

	barcodePrefix <- NULL
	if (!is.null(cellbarcodePrefix)) {
		barcodePrefix <- cellbarcodePrefix
	} else if ('DatasetId' %in% names(seuratObj@meta.data)) {
		datasetId <- unique(seuratObj$DatasetId)
		if (length(datasetId) != 1) {
			stop('This seurat object has multiple barcode prefixes / datasetIds. Must supply cellbarcodePrefix!')
		}

		barcodePrefix <- paste0(datasetId, '_')
	} else if ('BarcodePrefix' %in% names(seuratObj@meta.data)) {
		datasetId <- unique(seuratObj$BarcodePrefix)
		if (length(datasetId) != 1) {
			stop('This seurat object has multiple barcode prefixes / datasetIds. Must supply cellbarcodePrefix!')
		}

		barcodePrefix <- paste0(datasetId, '_')
	}

	if (!is.null(barcodePrefix)) {
		df$cellbarcode <- paste0(barcodePrefix, df$cellbarcode)
	}

	if (length(intersect(colnames(seuratObj), df$cellbarcode)) == 0) {
		bc1 <- paste0(head(df$cellbarcode, n = 2), collapse = ';')
		bc2 <- paste0(head(colnames(seuratObj), n = 2), collapse = ';')
		print(paste0('No overlapping barcodes found (example: ', bc1, ' / ', bc2,'), adding gem_group'))
		df$cellbarcode <- paste0(df$cellbarcode, '-', df$gem_group)
	}

	if (length(intersect(colnames(seuratObj), df$cellbarcode)) == 0) {
		bc1 <- paste0(head(df$cellbarcode, n = 2), collapse = ';')
		bc2 <- paste0(head(colnames(seuratObj), n = 2), collapse = ';')
		stop(paste0('No overlapping barcodes found between seuratObj and molecule_info.h5 file, example: ', bc1, ' / ', bc2))
	}

	df <- df[df$cellbarcode %in% colnames(seuratObj),c('cellbarcode', 'num.umis', 'num.reads')]
	df <- data.frame(cellbarcode = df$cellbarcode, num.umis = df$num.umis, CountsPerCell = df$num.reads)
	df$Saturation <- 1 - (df$num.umis / df$CountsPerCell)

	return(df)
}

.AppendSaturation <- function(seuratObj, df, assayName) {
	print(paste0('Adding saturation to assay ', assayName, ' for ', nrow(df), ' cells'))

	if (is.null(assayName)) {
		assayName <- Seurat::DefaultAssay(seuratObj)
	}

	targetField <- paste0('Saturation.', assayName)
	targetFieldReads <- paste0('nReads_', assayName)
	if (!(targetField %in% names(seuratObj@meta.data))) {
		seuratObj[[targetField]] <- NA
	}

	if (!(targetFieldReads %in% names(seuratObj@meta.data))) {
		seuratObj[[targetFieldReads]] <- NA
	}

	toMerge <- df$Saturation
	names(toMerge) <- df$cellbarcode

	d <- seuratObj[[targetField, drop = TRUE]]
	names(d) <- colnames(seuratObj)
	d[names(toMerge)] <- toMerge
	seuratObj <- Seurat::AddMetaData(seuratObj, metadata = d, col.name = targetField)

	toMerge <- df$CountsPerCell
	names(toMerge) <- df$cellbarcode

	d <- seuratObj[[targetFieldReads, drop = TRUE]]
	names(d) <- colnames(seuratObj)
	d[names(toMerge)] <- toMerge
	seuratObj <- Seurat::AddMetaData(seuratObj, metadata = d, col.name = targetFieldReads)

	return(seuratObj)
}

#' @title Append Per Cell Saturation
#' @description Calculate and Append Per-cell Saturation to a Seurat Object
#' @param seuratObj The seurat object
#' @param molInfoFile The 10x molecule_info.h5 file
#' @param cellbarcodePrefix An optional string appended to the barcodes parsed from the molecule_info.h5 file. This is necessary if the seurat object has a prefix applied to cell barcodes. This value is directly concatenated and must include any delimiter. Note: if this is absent, but the seuratObj has the columns DatasetId or BarcodePrefix, the latter will be used.
#' @param assayName An optional string indicating the assay these data are associated with. If null, they will use the DefaultAssay
#' @param doPlot If true, plots summarizing saturation will be generated.
#' @export
AppendPerCellSaturation <- function(seuratObj, molInfoFile, cellbarcodePrefix = NULL, assayName = NULL, doPlot = TRUE) {
	df <- .CalcPerCellSaturation(seuratObj, molInfoFile = molInfoFile, cellbarcodePrefix = cellbarcodePrefix)

	if (doPlot) {
		overall <- 1 - round((sum(df$num.umis) / sum(df$CountsPerCell)), 2)
		print(ggplot(df, aes(x = CountsPerCell, y = Saturation)) +
				  labs(x = 'Counts/Cell', y = '% Saturation') +
				  egg::theme_presentation(base_size = 18) +
				  geom_point() +
				  annotate("text", x = max(df$CountsPerCell), y = min(df$Saturation), hjust = 1, vjust = -1, label = paste0(
					  'Total Counts: ', format(sum(df$CountsPerCell), big.mark=','), '\n',
					  'UMI Counts: ', format(sum(df$num.umis), big.mark=','), '\n',
					  'Saturation: ', overall
				  )) + ggtitle('Library Saturation')
		)
	}

	seuratObj <- .AppendSaturation(seuratObj, df, assayName = assayName)
	return(seuratObj)
}

#' @title Append Per Cell Saturation in Bulk
#' @description Calculate and Append Per-cell Saturation to a Seurat Object, where the seurat object contains multiple input 10x runs
#' @param seuratObj The seurat object
#' @param molInfoList A list mapping dataset/assay to the filepath of the 10x molecule_info.h5 file. The format of the names should be DatasetId-Assayname (i.e. 283729-GEX or 384729-HTO)
#' @export
AppendPerCellSaturationInBulk <- function(seuratObj, molInfoList) {
	print('Adding saturation')

	if (!('DatasetId' %in% names(seuratObj@meta.data))) {
		stop('This seuratObj must have a DatasetId column')
	}

	uniqueAssays <- c()
	dataframes <- list()
	for (i in names(molInfoList)) {
		datasetId <- unlist(strsplit(i, split = '-'))[1]
		assayName <- unlist(strsplit(i, split = '-'))[2]

		if (!(datasetId %in% unique(seuratObj$DatasetId))) {
			print(paste0('Skipping dataset: ', i))
			next
		}

		print(paste0('Calculating saturation: ', i))
		uniqueAssays <- c(uniqueAssays, assayName)

		toAppend <- .CalcPerCellSaturation(seuratObj, molInfoList[[i]], cellbarcodePrefix = paste0(datasetId, '_'))
		if (!assayName %in% names(dataframes)) {
			dataframes[[assayName]] <- toAppend
		} else {
			if (all(is.null(dataframes[[assayName]]))) {
				dataframes[[assayName]] <- toAppend
			} else {
				dataframes[[assayName]] <- rbind(dataframes[[assayName]], toAppend)
			}
		}
	}

	if (length(uniqueAssays) == 0) {
		print('No datasets had saturation added')
	}

	for (assayName in unique(uniqueAssays)) {
		if (!is.null(dataframes[[assayName]]) && nrow(dataframes[[assayName]]) > 0) {
			seuratObj <- .AppendSaturation(seuratObj, dataframes[[assayName]], assayName = assayName)
		} else {
			print('No saturation values to add. This may indicate an error matching cell barcodes.')
		}

		fieldName <- paste0('Saturation.', assayName)
		readField <- paste0('nReads_', assayName)
		umiField <- paste0('nCount_', assayName)
		if (!umiField %in% names(seuratObj@meta.data)) {
			print(paste0('seurat object lacks the assay ', assayName, ', skipping plots'))
			next
		}

		if (!readField %in% names(seuratObj@meta.data)) {
			stop(paste0('Missing field: ', readField))
		}

		if (!fieldName %in% names(seuratObj@meta.data)) {
			stop(paste0('Missing field: ', fieldName))
		}

		if (!umiField %in% names(seuratObj@meta.data)) {
			stop(paste0('Missing field: ', umiField))
		}

		df <- data.frame(CountsPerCell = seuratObj@meta.data[[readField]], Saturation = seuratObj@meta.data[[fieldName]], num.umis = seuratObj@meta.data[[umiField]])
		if ('DatasetName' %in% names(seuratObj@meta.data)) {
			df$Label <- seuratObj@meta.data$DatasetName
		} else {
			df$Label <- seuratObj@meta.data$DatasetId
		}
		df$Label <- as.character(df$Label)

		overall <- 1 - round((sum(df$num.umis) / sum(df$CountsPerCell)), 2)
		print(ggplot(df, aes(x = CountsPerCell, y = Saturation, color = Label)) +
			  labs(x = 'Counts/Cell', y = '% Saturation') +
			  egg::theme_presentation(base_size = 18) +
			  geom_point() +
			  annotate("text", x = max(df$CountsPerCell), y = min(df$Saturation), hjust = 1, vjust = -1, label = paste0(
				  'Total Counts: ', format(sum(df$CountsPerCell), big.mark=','), '\n',
				  'UMI Counts: ', format(sum(df$num.umis), big.mark=','), '\n',
				  'Saturation: ', overall
			  )) + ggtitle(paste0('Library Saturation: ', assayName))
		)
	}

	return(seuratObj)
}


#' @title Plot Seurat Variables
#' @description Create a panel of plots summarizing two variables within the seurat object
#' @param seuratObj The seurat object
#' @param xvar The first variable to summarize (i.e. ClusterNames_0.2)
#' @param yvar The second variable to summarize (i.e. SubjectId)
#' @param labelDimplot This value is passed directly to the Seurat::DimPlot label argument
#' @param reduction The reduction to use for the Seurat::DimPlot
#' @export
PlotSeuratVariables <- function(seuratObj, xvar, yvar, labelDimplot = FALSE, reduction = 'umap') {
	data <- seuratObj@meta.data[c(xvar, yvar)]
	names(data) <- c('x', 'y')

	P0 <- DimPlot(seuratObj, group.by = xvar, label = labelDimplot, reduction = reduction, shuffle = TRUE)
	P2 <- DimPlot(seuratObj, group.by = yvar, label = labelDimplot, reduction = reduction, shuffle = TRUE)

	P1 <- ggplot(data, aes(x = x, fill = y)) +
		geom_bar(position = 'fill', color = 'black') +
		xlab(xvar) +
		labs(fill = yvar) +
		egg::theme_presentation(base_size = 12) +
		theme(
			axis.text.x = element_text(angle = 45, hjust = 1)
		)

	P3 <- ggplot(data, aes(x = x, fill = y)) +
		geom_bar(color = 'black') +
		xlab(xvar) +
		labs(fill = yvar) +
		egg::theme_presentation(base_size = 12) +
		theme(
			axis.text.x = element_text(angle = 45, hjust = 1)
		)

	return(P0 + P2 + P1 + P3)
}


#' @title Inspect Seurat
#'
#' @description Print various information about the size of a seurat object. Can be useful for debugging.
#' @param seuratObj The seurat object
#' @param slotReportSize The size in bytes, above which a slot's size will be logged.
#' @param commandReportSize The size in bytes, above which a command's size will be logged.
#' @importFrom methods is slot<- slotNames
#' @export
InspectSeurat <- function(seuratObj, slotReportSize = 500000, commandReportSize = 500000) {
	print(paste0('Seurat object size: ', format(utils::object.size(seuratObj), units = 'auto')))

	for (assayName in names(seuratObj@assays)) {
		print(paste0('Assay: ', assayName))
		for (slotName in c('counts', 'data', 'scale.data')) {
			dat <- Seurat::GetAssayData(seuratObj, assay = assayName, slot = slotName)
			print(paste0(' slot: ', slotName, ', size: ', format(utils::object.size(x = dat), units = 'auto')))

			if (slotName != 'scale.data' && !is(dat, 'sparseMatrix')) {
				print(paste0('Non-sparse! Assay: ', assayName, ', slot: ', slotName))
			}
		}
	}

	print('All slots:')
	for (slotName in methods::slotNames(seuratObj)) {
		val <- utils::object.size(methods::slot(seuratObj, slotName))
		if (val > slotReportSize) {
			print(paste0(' ', slotName, ': ', format(val, units = 'auto')))
		}
	}

	# Note: we have had historic issues with commands logging an enormous amount of data, so add extra reporting:
	for (commandName in names(seuratObj@commands)) {
		val <- utils::object.size(x = seuratObj@commands[commandName])
		if (val > commandReportSize) {
			print(paste0('Command: ', commandName, ', size: ', format(val, units = 'auto')))
			for (slotName in methods::slotNames(seuratObj@commands[[commandName]])) {
				val <- utils::object.size(x = slot(seuratObj@commands[[commandName]], slotName))
				if (val > commandReportSize) {
					print(paste0(' ', slotName, ', size: ', format(val, units = 'auto')))
				}
			}
		}
	}
}

.ClearSeuratCommands <- function(seuratObj, maxSize = 500000) {
	for (commandName in names(seuratObj@commands)) {
		val <- utils::object.size(x = slot(seuratObj@commands[[commandName]], 'call.string'))
		if (val > maxSize) {
			print(paste0('Clearing call.string for: ', commandName, '. size: ', format(val, units = 'auto')))
			slot(seuratObj@commands[[commandName]], 'call.string') <- ''
		}
	}

	return(seuratObj)
}

#' @title Scale Features If Needed
#'
#' @description This accepts a list of features and appends them to the target assay, if they are not already scaled.
#' @param seuratObj The seurat object
#' @param toScale The list of features to scale.
#' @param assayName The name of the assay to use.
#' @export
ScaleFeaturesIfNeeded <- function(seuratObj, toScale, assayName = 'RNA') {
	notPresent <- toScale[!toScale %in% rownames(Seurat::GetAssayData(seuratObj, assay = assayName, slot = 'scale.data'))]
	print(paste0('Total features to scale: ', length(notPresent), ' of ', length(toScale)))

	scaled2 <- Seurat::ScaleData(seuratObj@assays[[assayName]], features = notPresent)
	scaled2 <- rbind(Seurat::GetAssayData(seuratObj, assay = assayName, slot = 'scale.data'), Seurat::GetAssayData(scaled2, slot = 'scale.data'))
	seuratObj <- Seurat::SetAssayData(seuratObj, assay = assayName, slot = 'scale.data', new.data = scaled2)
	print(dim(Seurat::GetAssayData(seuratObj, assay = assayName, slot = 'scale.data')))

	return(seuratObj)
}


#' @title Create a dataframe with chi-squared statistics of clones
#' @description This function creates a dataframe with chi-squared statistics between two metadata fields
#' @param seuratObj The seurat object
#' @param field1 The name of the first column to compare
#' @param field2 The name of the second column to compare
#' @param plot If true, the function will plot a heatmap using pheatmap
#' @return A dataframe with chi-squared statistics of clones
#' @export
# TODO: fix this
GetChiDF <- function(seuratObj, field1, field2, plot = FALSE) {
	if (!field1 %in% names(seuratObj@meta.data)) {
		stop(paste0('Missing field: ', field1))
	}
	dat1 <- seuratObj@meta.data[,field1]

	if (!field2 %in% names(seuratObj@meta.data)) {
		stop(paste0('Missing field: ', field2))
	}
	dat2 <- seuratObj@meta.data[,field2]

	tempDF <- lapply(clone_vector, function(x) {
		stats::chisq.test(table(dat1, dat2))$res[2,]
	}) %>% as.data.frame()

	colnames(tempDF) <- dat1
	if (plot) {
		pheatmap::pheatmap(asinh(tempDF))
	}

	return(tempDF)
}

#' @title Add a new Seurat object metadata column based on current metadata
#' @description This function creates a new Seurat object metadata column based on current formulas applied to current metadata columns
#' @param seuratObj The seurat object
#' @param varname The name of the new metadata column
#' @param formulavector vector of formulas, e.g.: c(Tcell_NaiveToEffector > 10 ~ "Effector", Tcell_NaiveToEffector < 5 ~ "Naive")
#' @param defaultname The default value applied when none of the formulas in formulavec are TRUE
#' @param enforceFunctionalFormulae Boolean determining whether or not to sanity check that all defined formulae in formulavector labeled at least one cell in the Seurat Object.
#' @importFrom rlang :=
#' @return Updated Seurat object
#' @export
AddNewMetaColumn <- function(seuratObj, varname, formulavector, defaultname, enforceFunctionalFormulae = TRUE) {
  seuratObj@meta.data <- seuratObj@meta.data |> mutate(
    "{varname}" := case_when(
      !!!rlang::parse_exprs(paste0(formulavector)),
      .default = defaultname)
  )
  
  if (enforceFunctionalFormulae) {
    unique_values <- list()
    #figure out the unique values supplied to the formula vector
    for (formula in formulavector){
        #pasting a RHS formula looks like: "~" ,"Tcell_NaiveToEffector > 10", "Effector"     
        value <- paste0(formula)[[3]]
      unique_values <- append(unique_values, value)
    }
    #add the default name
    unique_values <- append(unique_values, defaultname)
    #convert to vector
    unique_values <- unlist(unique_values)
    #check that all of those unique values exist in the metadata
    if (!all(unique_values %in% unique(seuratObj@meta.data[,varname]) )){
      stop("Not all of the values specified in the formulae exist in the newly defined metadata column. Please check that each of your formula in formulavector are valid.")
    }
  }
  return(seuratObj)
  
}