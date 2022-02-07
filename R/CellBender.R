#' @title Run CellBender
#'
#' @description Runs CellBender using the 10x h5 file as input
#' @param rawFeatureH5File The path to 10x's raw_feature_bc_matrix.h5
#' @param expectedCells Passed to CellBender --expected-cells
#' @param totalDropletsIncluded Passed to CellBender --total-droplets-included
#' @param fpr Passed to CellBender --fpr
#'
#' @export
RunCellBender <- function(rawFeatureMatrix, expectedCells = 5000, totalDropletsIncluded = 20000, fpr = 0.01, epochs = 150) {
  if (!reticulate::py_available(initialize = TRUE)) {
    stop(paste0('Python/reticulate not configured. Run "reticulate::py_config()" to initialize python'))
  }

  if (!reticulate::py_module_available('cellbender')) {
    stop('The cellbender python package has not been installed!')
  }

  print('Running cellbender:')
  inputh5File <- tempfile(fileext = '.h5')
  DropletUtils::write10xCounts(inputh5File, x = rawFeatureMatrix)
  outH5File <- tempfile(fileext = '.h5')

  # consider: --low-count-threshold
  system2("cellbender", c(
    "remove-background",
    "--expected-cells", expectedCells,
    "--total-droplets-included", totalDropletsIncluded,
    "--fpr", fpr,
    "--epochs", epochs,
    "--output", outH5File,
    "--input", inputh5File
  ))

  if (!file.exists(outH5File)) {
    stop(paste0('Missing file: ', outH5File))
  }

  seuatRawData <- Seurat::Read10X_h5(filename = outH5File, use.names = TRUE)

  #TODO: PDF

  unlink(inputh5File)
  unlink(outH5File)

  return(seuatRawData)
}