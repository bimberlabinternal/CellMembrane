
utils::globalVariables(
  names = c('GEXOutfile', 'modelBasePath', 'modelFileName', 'exclusionJsonPath', 'ptimeOutFile'),
  package = 'CellMembrane',
  add = TRUE
)

#' @title TrainSctourModel
#'
#' @description Trains an sctour model on a Seurat object
#' @param seuratObj The Seurat object containing the data. 
#' @param GEX_outfile The GEX filename used for output of DropletUtils::write10xCounts (.h5 extension).
#' @param meta_outfile The filename used for containing the seuratObj metadata.
#' @param assayName Assay whose data is to be written out with DropletUtils::write10xCounts. Should be "RNA".
#' @param exclusionList List of genes to be excluded from variable genes in sctour. Usually a RIRA GeneSet.
#' @param exclusion_JSON_path Filename for the file containing the gene exclusion list (.json extension)
#' @return A sctour model (saved as a pytorch model). 
#' @examples
#' \dontrun{
#' seuratObj <- TrainSctourModel(seuratObj = seuratObj, 
#'                 GEXOutfile = './gex_tempfile.h5', 
#'                 modelBasePath =  './', 
#'                 modelFileName = "test_model",
#'                 exclusionJsonPath = './exclusion_tempfile.json',
#'                 ptimeOutFile = './ptime_out_file.csv',
#'                 variablesGenesFile = './variable_genes_out_file.csv',
#'                 assayName = "RNA",
#'                 embeddingOutFile = "./embeddings.csv",
#'                 cleanUpIntermediateFiles = T)
#'
#' }
#' @export

TrainSctourModel <- function(seuratObj = NULL,
                             GEXOutfile = NULL,
                             exclusionList = NULL,
                             exclusionJsonPath = NULL,
                             modelBasePath = NULL,
                             modelFileName = NULL,
                             ptimeOutFile = NULL,
                             variablesGenesFile = NULL,
                             embeddingOutFile = NULL, 
                             assayName = "RNA", 
                             cleanUpIntermediateFiles = T) {
  #TODO: Sanitize inputs
  
  # write out data for scTour 
  jsonlite::write_json(exclusionList, exclusion_JSON_path)
  DropletUtils::write10xCounts(x = Seurat::GetAssayData(seuratObj, assay = assayName, layer = 'counts'), 
                               path = GEX_outfile,
                               overwrite = TRUE)
  
  #copy run_scTour.py in inst/scripts and supply custom arguments 
  str <- readr::read_file(system.file("scripts/TrainScTourModel.py", package = "CellMembrane"))
  script <- tempfile()
  readr::write_file(str, script)

  #train and write the model file and accessory variable genes
  newstr <- paste0("TrainScTourModel(GEXfile = '", R.utils::getAbsolutePath(GEXOutfile, mustWork = FALSE),
                   "', exclusion_json_path = '", R.utils::getAbsolutePath(exclusionJsonPath, mustWork = FALSE),
                   "', model_path_basedir = '", R.utils::getAbsolutePath(modelBasePath, mustWork = FALSE),
                   "', model_name = '", R.utils::getAbsolutePath(modelFileName, mustWork = FALSE),
                   "', embedding_out_file = '", R.utils::getAbsolutePath(embeddingOutFile, mustWork = FALSE),
                   "', ptime_out_file = '", R.utils::getAbsolutePath(ptimeOutFile, mustWork = FALSE),
                   "', variable_genes_out_file = '", R.utils::getAbsolutePath(variablesGenesFile, mustWork = FALSE),
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
  seuratObj[['sctour']] <- CreateDimReducObject(embeddings = embeddings, key = "SCTOUR_", assay = assayName)
  
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
    unlink(c(GEXOutfile, ptimeOutFile, embeddingOutFile, exclusion_JSON_path))
  }
  
  return(seuratObj)
}

#' @title PredictScTourPseudotime
#'
#' @description Predicts pseudotime from a trained sctour model.
#' @param seuratObj The Seurat object containing the data. 
#' @param GEX_outfile The GEX filename used for output of DropletUtils::write10xCounts (.h5 extension).
#' @param meta_outfile The filename used for containing the seuratObj metadata.
#' @param assayName Assay whose data is to be written out with DropletUtils::write10xCounts. Should be "RNA".
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
  DropletUtils::write10xCounts(x = Seurat::GetAssayData(seuratObj, 
                                                        assay = assayName, 
                                                        layer = 'counts'), 
                               path = GEX_outfile,
                               overwrite = TRUE)
  
  #copy run_scTour.py in inst/scripts and supply custom arguments 
  str <- readr::read_file(system.file("scripts/PredictScTourPseudotime.py", package = "CellMembrane"))
  script <- tempfile()
  readr::write_file(str, script)
  
  newstr <- paste0("PredictPseudotime(GEXfile = '", GEXOutfile,
                   "', model_file = '", R.utils::getAbsolutePath(modelFile, mustWork = FALSE),
                   "', variable_genes_file = '", R.utils::getAbsolutePath(variablesGenesFile, mustWork = FALSE), 
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
    unlink(c(GEXOutfile, ptimeOutFile, embeddingOutFile, exclusion_JSON_path))
  }
  
  return(seuratObj)
}
