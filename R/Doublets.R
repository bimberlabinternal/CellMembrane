
FindDoublets <- function(seuratObj) {
	SerObj_CRS.sce <- as.SingleCellExperiment(SerObj_CRS)
	SerObj_CRS.sce <- scDblFinder:::scDblFinder(SerObj_CRS.sce)
}