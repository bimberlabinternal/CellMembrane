
utils::globalVariables(
  names = c('GEX_outfile', 'meta_outfile', 'exclusion_JSON_path', 'ptime_file'),
  package = 'CellMembrane',
  add = TRUE
)

#' @title Run sctour
#'
#' @description Runs sctour on a seurat object
#' @param seuratObj The Seurat object containing the data to be run using sctour
#' @param GEX_outfile The GEX filename used for output of DropletUtils::write10xCounts (.h5 extension).
#' @param meta_outfile The filename used for containing the seuratObj metadata.
#' @param assayName Assay whose data is to be written out with DropletUtils::write10xCounts. Should be "RNA"
#' @param exclusionList List of genes to be excluded from variable genes in sctour. Usually a RIRA GeneSet.
#' @param exclusion_JSON_path Filename for the file containing the gene exclusion list (.json extension)
#' @param ptime_file Filename for the pseudotime results to be written to (.csv extension)
#' @return A Seurat object with the pseudotime results appended. 
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

Run_sctour <- function(seuratObj=NULL,
                     GEX_outfile = NULL,
                     meta_outfile = NULL, 
                     assayName = "RNA", 
                     exclusionList = NULL,
                     exclusion_JSON_path = NULL,
                     ptime_file = NULL) {
  
  # write out data for scTour 
  write_json(exclusionList, exclusion_JSON_path)
  DropletUtils::write10xCounts(x = seuratObj@assays[[assayName]]@counts, 
                               path = GEX_outfile,
                               overwrite = TRUE)
  write.csv(seuratObj@meta.data, meta_outfile)
  
  #copy run_scTour.py in inst/scripts and supply custom arguments 
  str <- readr::read_file(system.file("scripts/run_scTour.py", package = "CellMembrane"))
  script <- tempfile()
  readr::write_file(str, script)
  
  newstr <- paste0("run_sctour(GEXfile = '", GEX_outfile,
                   "',metafile = '", meta_outfile,
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


