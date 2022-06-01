#' @title Run VISION using MSigDB Gene Sets
#'
#' @description Run VISION using MSigDB Gene Sets. Returns an object that can be inspected using: VISION::viewResults(vision.obj)
#' @param seuratObj The seurat object
#' @param metadataCols A list of column names to include in VISION analysis.
#' @import Seurat
#' @export
RunVisionForMSigDB <- function(seuratObj, metadataCols = c('nCount_RNA', 'ClusterNames_0.2', 'ClusterNames_0.4', 'ClusterNames_0.6', 'ClusterNames_0.8')) {
  print('Downloading msigdb.v7.5.1')
  destFile <- tempfile(fileext = ".gmt")
  utils::download.file('https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/msigdb.v7.5.1.symbols.gmt', destfile = destFile)

  seuratObj@meta.data <- seuratObj@meta.data[,names(seuratObj@meta.data) %in% metadataCols]

  # see: https://yoseflab.github.io/VISION/articles/web_only/Seurat.html
  vision.obj <- VISION::Vision(seuratObj, signatures = destFile, projection_methods = NULL)
  vision.obj <- VISION::analyze(vision.obj)

  return(vision.obj)
}