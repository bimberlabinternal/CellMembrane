#' @title Find Doublets using scDblFinder
#'
#' @description Find Doublets using scDblFinder
#' @param seuratObj The seurat object
#' @param assay The assay to use
#' @param rawResultFile An optional path to a file where raw results will be written
#' @param doPlot If true, a DimPlot will be generated, assuming reductions are calculated
#' @param dropDoublets If true, any cells marked as doublets will be dropped
#' @import Seurat
#' @export
FindDoublets <- function(seuratObj, assay = 'RNA', rawResultFile = NULL, doPlot = TRUE, dropDoublets = FALSE) {
	sce <- Seurat::as.SingleCellExperiment(seuratObj, assay = assay)
	df <- scDblFinder::scDblFinder(sce, returnType = 'table')
	df <- df[df$type == 'real',]

	toAdd <- df$class
	names(toAdd) <- rownames(df)

	# Sanity check:
	if (sum(names(toAdd) != colnames(seuratObj)) > 0) {
		stop('The cell barcode names do not match')
	}

	seuratObj[['scDblFinder.class']] <- toAdd

	if (!is.null(rawResultFile)) {
		df$cellbarcode <- rownames(df)
		write.table(df, file = rawResultFile, sep = '\t', row.names = FALSE)
	}

	if (doPlot) {
		if (length(seuratObj@reductions) == 0) {

		} else {
			print(DimPlot(seuratObj, group.by = 'scDblFinder.class'))
		}
	}

	if (dropDoublets) {
		nDoublet = sum(seuratObj[['scDblFinder.class']] == 'doublet')
		print('Dropping ', nDoublet, ' doublets')
		seuratObj <- subset(seuratObj, subset = scDblFinder.class != 'doublet')
	}

	return(seuratObj)
}