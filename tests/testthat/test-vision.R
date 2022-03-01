context("scRNAseq")

test_that("Doublet detection works as expected", {
  seuratObj <- readRDS('../testdata/seuratOutput.rds')
  vision.out <- RunVisionForMSigDB(seuratObj)

  #expect_equal(1540, ncol(seuratObj), tolerance = 3)
})

