context("scRNAseq")

test_that("TrainSctourModel works", {
  if (!reticulate::py_module_available('sctour')) {
    print('sctour module not found, debugging:')
    print(reticulate::py_list_packages())
    if ('sctour' %in% reticulate::py_list_packages()$package) {
      tryCatch({
        reticulate::import('sctour')
      }, error = function(e){
        print("Error with reticulate::import('sctour')")
        print(conditionMessage(e))
        traceback()
      })
    }
    
    warning('The python sctour module has not been installed!')
    return()
  }
  
  #read data
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS("../testdata/seuratOutput.rds")))
  
  
  TrainSctourModel(seuratObj = seuratObj,
                               modelFile = "./sctour_model.pth",
                               variablesGenesFile = './variableGenes.csv',
                               overwrite = TRUE, 
                               GEX_outfile = "./gex_tempfile.h5",
                               meta_outfile = "./metadata_tempfile.csv", 
                               assayName = "RNA", 
                               exclusionList = RIRA::GetGeneSet('VariableGenes_Exclusion.2'),
                               exclusion_JSON_path = "./exclusion_tempfile.json", 
                               appendPseudotimeToSeurat = F)
  
  #test that a model file was written by sctour
  ##testthat::expect_true(file.exists("./sctour_model.pth"))
  #test the model's associated gene space was also written
  ##testthat::expect_true(file.exists("./variableGenes.csv"))
  #file cleanup
  ##unlink(c("./sctour_model.pth", "./variableGenes.csv"))
})
