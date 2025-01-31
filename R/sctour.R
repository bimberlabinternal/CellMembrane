
utils::globalVariables(
  names = c('GEXOutfile', 'modelBasePath', 'modelFileName', 'exclusionJsonPath', 'ptimeOutFile'),
  package = 'CellMembrane',
  add = TRUE
)

#' @title TrainSctourModel
#'
#' @description Trains an sctour model on a Seurat object
#' @param seuratObj The Seurat object containing the data. 
#' @param modelFileName The model's file name. The full name of the model will be modelFileName.pth
#' @param featureExclusionList A vector of gene names to be excluded from variable feature selection in model training. This supports RIRA's ExpandGeneList
#' @param outputBasePath A directory that will store the resulting pytorch model after model training.
#' @param assayName Assay whose data is to be written out with DropletUtils::write10xCounts. Should be "RNA".
#' @param cleanUpIntermediateFiles This boolean controls if GEXOutfile, embeddingOutFile, and exclusionJsonPath should be deleted after model training.
#' @param outputReductionName The assay name in which to store the results
#' @param metadataColName The name of the meta.data column to store the resulting pseudotime
#' @return A seurat object with pseudotime and dimensional reductions computed by scTour
#' @importFrom jsonlite write_json
#' @export
TrainSctourModel <- function(seuratObj,
                             modelFileName,
                             featureExclusionList = 'VariableGenes_Exclusion.2',
                             outputBasePath = './',
                             assayName = "RNA",
                             cleanUpIntermediateFiles = TRUE,
                             outputReductionName = 'sctour',
                             metadataColName = "pseudotime") {

  if (!reticulate::py_available(initialize = TRUE)) {
    stop(paste0('Python/reticulate not configured. Run "reticulate::py_config()" to initialize python'))
  }

  if (!reticulate::py_module_available('sctour')) {
    stop('The sctour python package has not been installed!')
  }

  # Adjust for windows paths:
  outputBasePath <- gsub(R.utils::getAbsolutePath(outputBasePath), pattern = '\\\\', replacement = '/')
  if (!endsWith(outputBasePath, '/')) {
    outputBasePath <- paste0(outputBasePath, '/')
  }

  #iterate over supplied gene sets and construct a gene set exclusion list
  exclusionList <- RIRA::ExpandGeneList(featureExclusionList)

  # write out data for scTour
  exclusionJsonPath <- R.utils::getAbsolutePath(paste0(outputBasePath, 'exclusionList.json'), mustWork = FALSE)
  jsonlite::write_json(exclusionList, exclusionJsonPath)

  GEXOutfile <- R.utils::getAbsolutePath(paste0(outputBasePath, 'gex.h5'), mustWork = FALSE)
  DropletUtils::write10xCounts(x = Seurat::GetAssayData(seuratObj, assay = assayName, layer = 'counts'),
                               path = GEXOutfile,
                               overwrite = TRUE)
  
  #copy run_scTour.py in inst/scripts and supply custom arguments 
  str <- readr::read_file(system.file("scripts/TrainScTourModel.py", package = "CellMembrane"))
  script <- tempfile()
  readr::write_file(str, script)

  ptimeOutFile <- R.utils::getAbsolutePath(paste0(outputBasePath, 'ptime.csv'), mustWork = FALSE)
  embeddingOutFile <- R.utils::getAbsolutePath(paste0(outputBasePath, 'embeddings.csv'), mustWork = FALSE)

  #train and write the model file and accessory variable genes
  newstr <- paste0("TrainScTourModel(GEXfile = '", GEXOutfile,
                   "', exclusion_json_path = '", exclusionJsonPath,
                   "', model_path_basedir = '", outputBasePath,
                   "', model_name = '", modelFileName,
                   "', embedding_out_file = '", embeddingOutFile,
                   "', ptime_out_file = '", ptimeOutFile,
                   "', random_state = ", GetSeed(),
                   ")")

  readr::write_file('\n', script, append = TRUE)
  readr::write_file(newstr, script, append = TRUE)

  system2(reticulate::py_exe(), script)

  seuratObj <- .AppendScTourAsReduction(seuratObj, embeddingOutFile, ptimeOutFile, outputReductionName, assayName, metadataColName)
  
  if (cleanUpIntermediateFiles){
    unlink(c(GEXOutfile, ptimeOutFile, embeddingOutFile, exclusionJsonPath))
  }
  
  return(seuratObj)
}

## Add ScTour Dimensional Reduction
.AppendScTourAsReduction <- function(seuratObj, embeddingOutFile, ptimeOutFile, outputReductionName, assayName, metadataColName) {
  for (fn in c(embeddingOutFile, ptimeOutFile)) {
    if (!file.exists(fn)) {
      stop('Missing file: ', fn, '. Found: ', paste0(list.files(dirname(fn)), collapse = ';'))
    }
  }

  #read the embeddings from the sctour training
  embeddings <- as.matrix(read.csv(embeddingOutFile, header = F))

  #get the cell ordering from the pseudotime output file
  pseudotimeOutputDf <- read.csv(ptimeOutFile)
  cellbarcodes <- pseudotimeOutputDf[,"X"]
  rownames(embeddings) <- cellbarcodes

  seuratObj[[outputReductionName]] <- Seurat::CreateDimReducObject(embeddings = embeddings, key = paste0(outputReductionName, "_"), assay = assayName)

  # Basic nearest neighbors graph construction to get a UMAP. Note: the embedding is five dimensaional from sctour .
  seuratObj <- Seurat::FindNeighbors(seuratObj, reduction = outputReductionName, dims = 1:5, k.param = 20)
  seuratObj <- Seurat::RunUMAP(seuratObj, dims = 1:2, reduction = outputReductionName, reduction.name = paste0(outputReductionName, "_umap"))

  #Add ptime to seurat metadata
  pseudotimeOutputVector <- pseudotimeOutputDf[,"ptime"]
  names(pseudotimeOutputVector) <- pseudotimeOutputDf[,"X"]
  seuratObj <- Seurat::AddMetaData(seuratObj, metadata = pseudotimeOutputVector, col.name = metadataColName)

  #Plot pseudotime
  print(FeaturePlot(seuratObj, features = metadataColName, reduction = paste0(outputReductionName, "_umap")))
  print(FeaturePlot(seuratObj, features = metadataColName, reduction = outputReductionName))

  return(seuratObj)
}

#' @title PredictScTourPseudotime
#'
#' @description Predicts pseudotime from a trained sctour model.
#' @param seuratObj The Seurat object containing the data. 
#' @param modelFile A path pointing to a trained scTour model.
#' @param outputBasePath The directory where the embeddings and pseudotime will be written. These are deleted unless cleanUpIntermediateFiles=TRUE
#' @param outputReductionName The assay name in which to store the results
#' @param metadataColName The name of the meta.data column to store the resulting pseudotime
#' @param assayName Assay whose data is to be written out with DropletUtils::write10xCounts. Should be "RNA".
#' @param cleanUpIntermediateFiles This boolean controls if GEXOutfile, embeddingOutFile, and exclusionJsonPath should be deleted after model training.
#' @return A Seurat Object with pseudotime from the supplied model written into its metadata.
#' @export
PredictScTourPseudotime <- function(seuratObj,
                                    modelFile,
                                    outputBasePath = tempdir(),
                                    outputReductionName = "sctour",
                                    metadataColName = "pseudotime",
                                    assayName = "RNA",
                                    cleanUpIntermediateFiles = TRUE) {

  if (!reticulate::py_available(initialize = TRUE)) {
    stop(paste0('Python/reticulate not configured. Run "reticulate::py_config()" to initialize python'))
  }

  if (!reticulate::py_module_available('sctour')) {
    stop('The sctour python package has not been installed!')
  }

  # Adjust for windows paths:
  outputBasePath <- gsub(R.utils::getAbsolutePath(outputBasePath), pattern = '\\\\', replacement = '/')
  if (!endsWith(outputBasePath, '/')) {
    outputBasePath <- paste0(outputBasePath, '/')
  }

  # write out data for scTour
  GEXOutfile <- R.utils::getAbsolutePath(paste0(outputBasePath, 'gex.h5'), mustWork = FALSE)
  DropletUtils::write10xCounts(x = Seurat::GetAssayData(seuratObj, assay = assayName, layer = 'counts'), path = GEXOutfile, overwrite = TRUE)
  
  #copy run_scTour.py in inst/scripts and supply custom arguments 
  str <- readr::read_file(system.file("scripts/PredictScTourPseudotime.py", package = "CellMembrane"))
  script <- tempfile()
  readr::write_file(str, script)

  ptimeOutFile <- R.utils::getAbsolutePath(paste0(outputBasePath, 'ptime.csv'), mustWork = FALSE)
  embeddingOutFile <- R.utils::getAbsolutePath(paste0(outputBasePath, 'embeddings.csv'), mustWork = FALSE)

  # See: https://dev-discuss.pytorch.org/t/bc-breaking-change-torch-load-is-being-flipped-to-use-weights-only-true-by-default-in-the-nightlies-after-137602/2573
  Sys.setenv('TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD', 1)

  newstr <- paste0("PredictPseudotime(GEXfile = '", GEXOutfile,
                   "', model_file = '", R.utils::getAbsolutePath(modelFile, mustWork = FALSE),
                   "', ptime_out_file = '", ptimeOutFile,
                   "', embedding_out_file = '", embeddingOutFile,
                   "')")

  #write the new arguments into the script and execute
  readr::write_file(newstr, script, append = TRUE)
  system2(reticulate::py_exe(), script)

  seuratObj <- .AppendScTourAsReduction(seuratObj, embeddingOutFile, ptimeOutFile, outputReductionName, assayName, metadataColName)

  if (cleanUpIntermediateFiles){
    unlink(c(GEXOutfile, ptimeOutFile, embeddingOutFile))
  }
  
  return(seuratObj)
}
