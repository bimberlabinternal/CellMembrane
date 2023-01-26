#' @include Utils.R
#' @include Preprocessing.R
#' @import Seurat

# utils::globalVariables(
#   names = c(''),
#   package = 'CellMembrane',
#   add = TRUE
# )


#' @title DoLdaParameterScan
#'
#' @description This will run LDA on the target assay
#' @param seuratObj A Seurat object.
#' @param outputFolder The path to save results. There will be subfolders for ./rawData and ./results
#' @param ntopics Passed to TITAN::runLDA.
#' @param normalizationMethod The method used for Seurat::NormalizeData()
#' @param varFeatures The number of variable features to use in the LDA model. The more features that are used, the slower the model will run and the more noise that will be introduced, but the model will be more complete in representing your entire dataset.
#' @param randomSeed Passed to TITAN::runLDA seed.number argument
#' @param assayName The name of the source assay
#' @param nCores The number of cores to use
#' @export
DoLdaParameterScan <- function(seuratObj, outputFolder, ntopics = seq(5, 50, by=5), normalizationMethod = "CLR", varFeatures = 5000, randomSeed = GetSeed(), assayName = 'RNA', nCores = 1) {
  # Perform normalization once:
  seuratObj <- Seurat::NormalizeData(seuratObj, assay = assayName, normalization.method = normalizationMethod)
  seuratObj <- Seurat::FindVariableFeatures(seuratObj, assay = assayName, nfeatures = varFeatures)

  TITAN::runLDA(seuratObj, ntopics = ntopics, normalizationMethod = normalizationMethod, seed.number = randomSeed, parallel = TRUE, outDir = outputFolder, cores = nCores, skipNormalization = TRUE)
  if (length(ntopics) > 1) {
    print(TITAN::LDAelbowPlot(outputFolder, seuratObj, skipNormalization = TRUE))
  }
}