#' @title Run CellBender
#'
#' @description Runs CellBender using the 10x h5 file as input
#' @param rawFeatureMatrix The path to 10x's raw_feature_bc_matrix.h5
#' @param expectedCells Passed to CellBender --expected-cells
#' @param totalDropletsIncluded Passed to CellBender --total-droplets-included
#' @param fpr Passed to CellBender --fpr
#' @param epochs Passed to CellBender --epochs
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
  args <- c(
    "remove-background",
    "--expected-cells", expectedCells,
    "--total-droplets-included", totalDropletsIncluded,
    "--fpr", fpr,
    "--epochs", epochs,
    "--output", outH5File,
    "--input", inputh5File
  )

  isEnabled <- system2(command = reticulate::py_exe(), shQuote(c("-c", 'import torch; print(torch.cuda.is_available())')), stdout = TRUE, stderr = TRUE)
  print(paste0('python torch.cuda.is_available(): ', isEnabled))

  if ("True" == isEnabled) {
    print('CUDA support detected, adding --cuda')
    args <- c(args, "--cuda")
  } else {
    print('CUDA support not detected, will not add --cuda')
  }

  system2("cellbender", args)

  outputFiltered <- gsub(outH5File, pattern = '.h5$', replacement = '_filtered.h5')
  if (!file.exists(outputFiltered)) {
    stop(paste0('Missing file: ', outputFiltered))
  }
  seuatRawData <- Seurat::Read10X_h5(filename = outputFiltered, use.names = TRUE)
  print(paste0('Cells in cellbender filtered matrix: ', ncol(seuatRawData)))

  outputPdf <- gsub(outH5File, pattern = '.h5$', replacement = '.pdf')
  plot(magick::image_read_pdf(outputPdf))

  unlink(inputh5File)
  unlink(outH5File)
  unlink(outputFiltered)
  unlink(gsub(outH5File, pattern = '.h5$', replacement = '.log'))
  unlink(gsub(outH5File, pattern = '.h5$', replacement = '_cell_barcodes.csv'))

  return(seuatRawData)
}