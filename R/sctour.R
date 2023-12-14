
utils::globalVariables(
  names = c('GEX_outfile', 'meta_outfile', 'exclusion_JSON_path', 'ptime_file'),
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
#'   CreateMergedTcrClonotypeFile(seuratObj,
#'                                outputFile = "./tcrs.csv")
#'   sctourSeuratObj <- Run_sctour(seuratObj = seuratObj,
#'                              GEX_outfile = "./GEX.h5",
#'                              meta_outfile = ./meta.csv, 
#'                              assayName = "RNA", 
#'                              exclusionList = RIRA::GetGeneSet('VariableGenes_Exclusion.2'),
#'                              exclusion_JSON_path = ./exclusionList.json,
#'                              ptime_file = NULL)
#' }
#' @export

TrainSctourModel <- function(seuratObj = NULL,
                       modelFile = NULL,
                       variablesGenesFile = NULL, 
                       overwrite = FALSE, 
                       GEX_outfile = NULL,
                       meta_outfile = NULL, 
                       assayName = "RNA", 
                       exclusionList = NULL,
                       exclusion_JSON_path = NULL, 
                       appendPseudotimeToSeurat = T) {
  
  # write out data for scTour 
  jsonlite::write_json(exclusionList, exclusion_JSON_path)
  DropletUtils::write10xCounts(x = Seurat::GetAssayData(seuratObj, assay = assayName, layer = 'counts'), 
                               path = GEX_outfile,
                               overwrite = TRUE)
  write.csv(seuratObj@meta.data, meta_outfile)
  
  #copy run_scTour.py in inst/scripts and supply custom arguments 
  #str <- readr::read_file(system.file("scripts/TrainScTour.py", package = "CellMembrane"))
  str <- "../../inst/scripts/TrainScTourModel.py"
  script <- tempfile()
  readr::write_file(str, script)
  #train and write the model file and accessory variable genes
  newstr <- paste0("TrainScTourModel(GEXfile = '", R.utils::getAbsolutePath(GEX_outfile, mustWork = FALSE),
                   "', metafile = '", R.utils::getAbsolutePath(meta_outfile, mustWork = FALSE),
                   "', exclusion_json_path = '", R.utils::getAbsolutePath(exclusion_JSON_path, mustWork = FALSE),
                   "', model_out_file = '", R.utils::getAbsolutePath(modelFile, mustWork = FALSE),
                   "', variable_genes_out_file = '", R.utils::getAbsolutePath(variablesGenesFile, mustWork = FALSE),
                   "')")
  readr::write_file(newstr, script, append = TRUE)
  system2(reticulate::py_exe(), script)
  #if requested, predict pseudotime values to be appended to the training data in the Seurat Object. 
  if (appendPseudotimeToSeurat) { 
    seuratObj <- PredictPseudotime(seuratObj = seuratObj, 
                                   modelFile = NULL,
                                   
                                   )
    }
  return(seuratObj)
}

#' @title PredictPseudotime
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
#'   CreateMergedTcrClonotypeFile(seuratObj,
#'                                outputFile = "./tcrs.csv")
#'   sctourSeuratObj <- Run_sctour(seuratObj = seuratObj,
#'                              GEX_outfile = "./GEX.h5",
#'                              meta_outfile = ./meta.csv, 
#'                              assayName = "RNA", 
#'                              exclusionList = RIRA::GetGeneSet('VariableGenes_Exclusion.2'),
#'                              exclusion_JSON_path = ./exclusionList.json,
#'                              ptime_file = NULL)
#' }
#' @export

PredictPseudotime <- function(seuratObj = NULL,
                             modelFile = NULL,
                             overwrite = FALSE, 
                             GEX_outfile = NULL,
                             meta_outfile = NULL, 
                             assayName = "RNA", 
                             exclusionList = NULL,
                             exclusion_JSON_path = NULL,
                             ptime_file = NULL) {
  # write out data for scTour 
  jsonlite::write_json(exclusionList, exclusion_JSON_path)
  DropletUtils::write10xCounts(x = seuratObj@assays[[assayName]]@counts, 
                               path = GEX_outfile,
                               overwrite = TRUE)
  write.csv(seuratObj@meta.data, meta_outfile)
  
  #copy run_scTour.py in inst/scripts and supply custom arguments 
  str <- readr::read_file(system.file("scripts/PredictScTourPseudotime.py", package = "CellMembrane"))
  script <- tempfile()
  readr::write_file(str, script)
  
  newstr <- paste0("PredictPseudotime(GEXfile = '", GEX_outfile,
                   "', metafile = '", meta_outfile,
                   "', exclusion_json_path = '", exclusion_JSON_path,
                   "', ptime_out_file = '", R.utils::getAbsolutePath(ptime_file, mustWork = FALSE),
                   "')")
  
  #write the new arguments into the script and execute
  readr::write_file(newstr, script, append = TRUE)
  system2(reticulate::py_exe(), script)
  
  pseudotime_values <- read_csv(ptime_file) 
  pseudotime_values <- pseudotime_values %>% column_to_rownames('...1') 
  seuratObj <- AddMetaData(seuratObj, pseudotime_values)
  
  return(seuratObj)
}
