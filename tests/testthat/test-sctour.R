context("scRNAseq")
library(Seurat)

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
  
  #read data, subset for speed:
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS("../testdata/seuratOutput.rds")))
  seuratObj <- subset(seuratObj, cells = colnames(seuratObj)[1:500])

  seuratObj<- TrainSctourModel(seuratObj = seuratObj, 
                   modelFileName = "test_model",
                   outputBasePath =  './'
  )
  
  #test that a model file was written by sctour
  testthat::expect_true(file.exists("./test_model.pth"))
  testthat::expect_true("pseudotime" %in% colnames(seuratObj@meta.data))
  testthat::expect_true("sctour" %in% names(seuratObj@reductions))

  seuratObj <- PredictScTourPseudotime(seuratObj, modelFile = "./test_model.pth", outputReductionName = 'sctour2', metadataColName = 'pseudotime2')

  testthat::expect_true("pseudotime2" %in% colnames(seuratObj@meta.data))
  testthat::expect_true("sctour2" %in% names(seuratObj@reductions))
  testthat::expect_equal(max(seuratObj$pseudotime2), 0.961, tolerance = 0.65)

  #file cleanup
  unlink("./sctour_model.pth")
})
