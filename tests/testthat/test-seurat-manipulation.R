context("scRNAseq")

test_that("SingleR works as expected", {
  set.seed(1234)

  seuratObj <- readRDS('../testdata/seuratOutput.rds')

  seuratObj2 <- DownsampleSeurat(seuratObj, targetCells = 500)
  expect_equal(500, ncol(seuratObj2))
  
  seuratObj2 <- DownsampleSeurat(seuratObj, targetCells = 100, subsetField = 'ClusterNames_0.2')
  expect_equal(493, ncol(seuratObj2))
  
  #seuratObj2 <- SubsetSeurat(seuratObj, expressionStrings = c('ClusterNames_0.2 == 0'))

  #seuratList <- SplitSeurat(seuratObj, splitField = 'ClusterNames_0.2')
  #expect_equal(6, length(seuratList))
})

