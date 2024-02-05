
utils::globalVariables(
  names = c('GEXOutfile', 'modelBasePath', 'modelFileName', 'exclusionJsonPath', 'ptimeOutFile'),
  package = 'CellMembrane',
  add = TRUE
)

#' @title TrainSctourModel
#'
#' @description Trains an sctour model on a Seurat object
#' @param seuratObj The Seurat object containing the data. 
#' @param featureExclusionList A vector of gene names to be excluded from variable feature selection in model training. This supports RIRA's ExpandGeneList
#' @param assayName Assay whose data is to be written out with DropletUtils::write10xCounts. Should be "RNA".
#' @param cleanUpIntermediateFiles This boolean controls if GEXOutfile, embeddingOutFile, and exclusionJsonPath should be deleted after model training.

#' @param modelBasePath A directory that will store the resulting pytorch model after model training. 
#' @param modelFileName The model's file name. The full name of the model will be modelFileName.pth 
#' @param GEXOutfile The GEX filename used for output of DropletUtils::write10xCounts (.h5 extension).
#' @param exclusionJsonPath Filename for the file containing the gene exclusion list (.json extension)
#' @param ptimeOutFile An output file that will contain a Nx2 csv where n is the number of cells in the seurat object, the first column is cell barcodes, and the second column is pseudotime as predicted by the model. 
#' @param variableGenesFile An scTour model requires both the trained model and feature space of the training data. This file stores the training data's feature space as a (N+1)x1 csv where N is the number of variable genes after exclusion list subtraction and the first row notes the column name "gene_ids".
#' @param embeddingOutFile The scTour model yields a dimensionally reduced space that a UMAP can be computed on. This file stores cell embeddings into this latent space, which are added as an "SCTOUR_" reduction into the input seurat object.
#' @param outputReductionName The assay name in which to store the results
#' @return A seurat object with pseudotime and dimensional reductions computed by scTour
#' @importFrom jsonlite write_json
#' @export
TrainSctourModel <- function(seuratObj,
                             modelFileName,
                             featureExclusionList = 'VariableGenes_Exclusion.2',
                             modelBasePath = './',
                             ptimeOutFile = NULL,
                             assayName = "RNA",
                             cleanUpIntermediateFiles = T,
                             outputReductionName = 'sctour') {

  if (!reticulate::py_available(initialize = TRUE)) {
    stop(paste0('Python/reticulate not configured. Run "reticulate::py_config()" to initialize python'))
  }

  #iterate over supplied gene sets and construct a gene set exclusion list
  exclusionList <- RIRA::ExpandGeneList(featureExclusionList)

  # write out data for scTour
  exclusionJsonPath <- R.utils::getAbsolutePath(paste0(modelBasePath, '/exclusionList.json'), mustWork = FALSE)
  jsonlite::write_json(exclusionList, exclusionJsonPath)

  GEXOutfile <- R.utils::getAbsolutePath(paste0(modelBasePath, '/gex.h5'), mustWork = FALSE)
  DropletUtils::write10xCounts(x = Seurat::GetAssayData(seuratObj, assay = assayName, layer = 'counts'),
                               path = GEXOutfile,
                               overwrite = TRUE)
  
  #copy run_scTour.py in inst/scripts and supply custom arguments 
  str <- readr::read_file(system.file("scripts/TrainScTourModel.py", package = "CellMembrane"))
  script <- tempfile()
  readr::write_file(str, script)

  ptimeOutFile <- R.utils::getAbsolutePath(paste0(modelBasePath, '/ptime.csv'), mustWork = FALSE)
  embeddingOutFile <- R.utils::getAbsolutePath(paste0(modelBasePath, '/embeddings.csv'), mustWork = FALSE)
  variableGenesFile <- R.utils::getAbsolutePath(paste0(modelBasePath, '/variableGenes.csv'), mustWork = FALSE)

  #train and write the model file and accessory variable genes
  newstr <- paste0("TrainScTourModel(GEXfile = '", GEXOutfile,
                   "', exclusion_json_path = '", exclusionJsonPath,
                   "', model_path_basedir = '", R.utils::getAbsolutePath(modelBasePath, mustWork = FALSE),
                   "', model_name = '", modelFileName,
                   "', embedding_out_file = '", embeddingOutFile,
                   "', ptime_out_file = '", ptimeOutFile,
                   "', variable_genes_out_file = '", variableGenesFile,
                   "')")
  readr::write_file(newstr, script, append = TRUE)
  readr::read_file(script)
  system2(reticulate::py_exe(), script)
  
  ## Add ScTour Dimensional Reduction
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
  seuratObj <- Seurat::AddMetaData(seuratObj, metadata = pseudotimeOutputVector, col.name = "pseudotime")

  #Plot pseudotime
  print(FeaturePlot(seuratObj, features = 'pseudotime', reduction = paste0(outputReductionName, "_umap")))
  print(FeaturePlot(seuratObj, features = 'pseudotime', reduction = outputReductionName))
  
  if (cleanUpIntermediateFiles){
    unlink(c(GEXOutfile, ptimeOutFile, embeddingOutFile, exclusionJsonPath))
  }
  
  return(seuratObj)
}

#' @title PredictScTourPseudotime
#'
#' @description Predicts pseudotime from a trained sctour model.
#' @param seuratObj The Seurat object containing the data. 
#' @param GEXOutfile The GEX filename used for output of DropletUtils::write10xCounts (.h5 extension).
#' @param assayName Assay whose data is to be written out with DropletUtils::write10xCounts. Should be "RNA".
#' @param modelFile A path pointing to a trained scTour model. 
#' @param variableGenesFile A path pointing to a trained scTour model's variable features. An scTour model requires both the trained model and feature space of the training data. This file stores the training data's feature space as a (N+1)x1 csv where N is the number of variable genes after exclusion list subtraction and the first row notes the column name "gene_ids".
#' @param embeddingOutFile The scTour model yields a dimensionally reduced space that a UMAP can be computed on. This file stores cell embeddings into this latent space, which are added as an "SCTOUR_" reduction into the input seurat object. 
#' @param ptimeOutFile An output file that will contain a Nx2 csv where n is the number of cells in the seurat object, the first column is cell barcodes, and the second column is pseudotime as predicted by the model. 
#' @param embeddingOutFile The scTour model will also compute cell embeddings and append them to the seurat object. These will be appended under the 'SCTOUR_" dimensional reduction. 
#' @param cleanUpIntermediateFiles This boolean controls if GEXOutfile, embeddingOutFile, and ptimeOutFile should be deleted after model training. 
#' @return A Seurat Object with pseudotime from the supplied model written into its metadata.  
#' @examples
#' \dontrun{
#'
#' seuratObj <- PredictScTourPseudotime(seuratObj,
#'                                      GEXOutfile = './gex_tempfile.h5',
#'                                      modelFile = "./test_model.pth",
#'                                      variableGenesFile = './variable_genes_out_file.csv',
#'                                      embeddingOutFile = "./embeddings.csv",
#'                                      ptimeOutFile ='./ptime_out_file.csv', 
#'                                      cleanUpIntermediateFiles = T)
#' }
#' @export

PredictScTourPseudotime <- function(seuratObj = NULL,
                                    GEXOutfile = NULL,
                                    modelFile = NULL,
                                    variableGenesFile = NULL, 
                                    ptimeOutFile = NULL,
                                    embeddingOutFile = NULL,
                                    assayName = "RNA", 
                                    cleanUpIntermediateFiles = T) {
  
  #TODO: Sanitize inputs
  
  # write out data for scTour 
  DropletUtils::write10xCounts(x = Seurat::GetAssayData(seuratObj, assay = assayName, layer = 'counts'), path = GEXOutfile, overwrite = TRUE)
  
  #copy run_scTour.py in inst/scripts and supply custom arguments 
  str <- readr::read_file(system.file("scripts/PredictScTourPseudotime.py", package = "CellMembrane"))
  script <- tempfile()
  readr::write_file(str, script)
  
  newstr <- paste0("PredictPseudotime(GEXfile = '", GEXOutfile,
                   "', model_file = '", R.utils::getAbsolutePath(modelFile, mustWork = FALSE),
                   "', variable_genes_file = '", R.utils::getAbsolutePath(variableGenesFile, mustWork = FALSE), 
                   "', ptime_out_file = '", R.utils::getAbsolutePath(ptimeOutFile, mustWork = FALSE),
                   "', embedding_out_file = '", R.utils::getAbsolutePath(embeddingOutFile, mustWork = FALSE),
                   "')")
  
  #write the new arguments into the script and execute
  readr::write_file(newstr, script, append = TRUE)
  system2(reticulate::py_exe(), script)
  
  ## Add ScTour Dimensional Reduction
  #read the embeddings from the sctour training
  embeddings <- as.matrix(read.csv(embeddingOutFile, header = F))
  #get the cell ordering from the pseudotime output file
  pseudotimeOutputDf <- read.csv(ptimeOutFile)
  cellbarcodes <- pseudotimeOutputDf[,"X"]
  rownames(embeddings) <- cellbarcodes
  seuratObj[['sctour']] <- CreateDimReducObject(embeddings = embeddings, key = "SCTOUR_", assay = assayName)
  print("scTour Embeddings added to Seruat object.")
  
  #basic nearest neighbors graph construction to get a UMAP. Note: the embedding is five dimensaional from sctour . 
  seuratObj <- Seurat::FindNeighbors(seuratObj, reduction = 'sctour', dims = 1:5, k.param = 20)
  seuratObj <- Seurat::RunUMAP(seuratObj, dims = 1:2, reduction = "sctour")
  
  #Add ptime to seurat metadata
  pseudotimeOutputVector <- pseudotimeOutputDf[,"ptime"]
  names(pseudotimeOutputVector) <- pseudotimeOutputDf[,"X"]
  seuratObj <- Seurat::AddMetaData(seuratObj, metadata = pseudotimeOutputVector, col.name = "pseudotime")
  #Plot pseudotime
  print(FeaturePlot(seuratObj, features = 'pseudotime', reduction = 'umap'))
  print(FeaturePlot(seuratObj, features = 'pseudotime', reduction = 'sctour'))
  
  if(cleanUpIntermediateFiles){
    unlink(c(GEXOutfile, ptimeOutFile, embeddingOutFile))
  }
  
  return(seuratObj)
}
