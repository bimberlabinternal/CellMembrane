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
                   GEXOutfile = './gex_tempfile.h5', 
                   modelBasePath =  './', 
                   modelFileName = "test_model",
                   exclusionJsonPath = './exclusion_tempfile.json',
                   ptimeOutFile = './ptime_out_file.csv',
                   variablesGenesFile = './variable_genes_out_file.csv',
                   assayName = "RNA",
                   embeddingOutFile = "./embeddings.csv",
                   cleanUpIntermediateFiles = T
  )
  
  
  #test that a model file was written by sctour
  testthat::expect_true(file.exists("./test_model.pth"))
  #test the model's associated gene space was also written
  testthat::expect_true(file.exists("./variable_genes_out_file.csv"))
  
  
  seuratObj <- PredictScTourPseudotime(seuratObj,
                                       GEXOutfile = './gex_tempfile.h5',
                                       modelFile = "./test_model.pth",
                                       variableGenesFile = './variable_genes_out_file.csv',
                                       embeddingOutFile = "./embeddings.csv",
                                       ptimeOutFile ='./ptime_out_file.csv', 
                                       cleanUpIntermediateFiles = T)
  testthat::expect_true("pseudotime" %in% colnames(seuratObj@meta.data))
  #file cleanup
  ##unlink(c("./sctour_model.pth", "./variableGenes.csv"))
})
