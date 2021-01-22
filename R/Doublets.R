#' @title Find Doublets using scDblFinder
#'
#' @description Find Doublets using scDblFinder
#' @param seuratObj The seurat object
#' @param assay The assay to use
#' @export
FindDoublets <- function(seuratObj, assay = 'RNA') {
	sce <- Seurat::as.SingleCellExperiment(seuratObj, assay = assay)
	df <- scDblFinder::scDblFinder(sce, returnType = 'table')

	#TODO
	return(df)
}