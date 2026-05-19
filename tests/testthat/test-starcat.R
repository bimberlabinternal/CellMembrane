context("scRNAseq")
library(Seurat)

test_that("RunStarCAT works", {
  if (!reticulate::py_module_available('starcat')) {
    print('starcat module not found, debugging:')
    print(reticulate::py_list_packages())
    if ('starcatpy' %in% reticulate::py_list_packages()$package) {
      tryCatch({
        reticulate::import('starcat')
      }, error = function(e) {
        print("Error with reticulate::import('starcat')")
        print(conditionMessage(e))
        traceback()
      })
    }
    warning('The python starcat (starcatpy) module has not been installed!')
    return()
  }

  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS("../testdata/seuratOutput.rds")))
  seuratObj <- subset(seuratObj, cells = colnames(seuratObj)[1:500])

  seuratObj <- RunStarCAT(seuratObj,
                          reference                = "TCAT.V1",
                          outputDirectory          = "./starcat_test_output",
                          name                     = "starcat_test",
                          cleanUpIntermediateFiles = FALSE)

  testthat::expect_true(file.exists("./starcat_test_output/starcat_test.scores.txt"))
  testthat::expect_true(file.exists("./starcat_test_output/starcat_test.rf_usage_normalized.txt"))
  testthat::expect_true(any(grepl("^starcat_score_", colnames(seuratObj@meta.data))))
  testthat::expect_true(any(grepl("^starcat_usage_", colnames(seuratObj@meta.data))))

  unlink("./starcat_test_output", recursive = TRUE)

  seuratObj <- RunStarCAT(seuratObj,
                          reference                = "TCAT.V1",
                          outputDirectory          = "./starcat_test_output_clean",
                          name                     = "starcat_test",
                          cleanUpIntermediateFiles = TRUE)

  testthat::expect_false(dir.exists("./starcat_test_output_clean"))
  testthat::expect_true(any(grepl("^starcat_score_", colnames(seuratObj@meta.data))))
  testthat::expect_true(any(grepl("^starcat_usage_", colnames(seuratObj@meta.data))))
})
