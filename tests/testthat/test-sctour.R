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

    stop('The python sctour module has not been installed!')
  }
  
  #read data
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS("../testdata/seuratOutput.rds")))

  TrainSctourModel(seuratObj = seuratObj, 
                   modelFileName = "test_model",
                   outputBasePath =  './'
  )
  
  
  #test that a model file was written by sctour
  testthat::expect_true(file.exists("./test_model.pth"))
  
  seuratObj <- PredictScTourPseudotime(seuratObj, modelFile = "./test_model.pth")

  testthat::expect_true("pseudotime" %in% colnames(seuratObj@meta.data))
  testthat::expect_true("sctour" %in% names(seuratObj@reductions))

  #file cleanup
  unlink("./sctour_model.pth")
})
