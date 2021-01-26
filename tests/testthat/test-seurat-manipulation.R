context("scRNAseq")

test_that("SingleR works as expected", {
  set.seed(1234)

  seuratObj <- readRDS('../testdata/seuratOutput.rds')

  seuratObj2 <- DownsampleSeurat(seuratObj, targetCells = 500)
  expect_equal(500, ncol(seuratObj2))
  
  seuratObj2 <- DownsampleSeurat(seuratObj, targetCells = 100, subsetField = 'ClusterNames_0.2')
  expect_equal(493, ncol(seuratObj2))

  seuratList <- SplitSeurat(seuratObj, splitField = 'ClusterNames_0.2')
  expect_equal(5, length(seuratList))

  seuratList <- SplitSeurat(seuratObj, splitField = 'ClusterNames_0.2', minCellsToKeep = 200)
  expect_equal(3, length(seuratList))

  #seuratObj2 <- SubsetSeurat(seuratObj, expressionStrings = c('ClusterNames_0.2 == 0'))
  #e <- expression(ClusterNames_0.2 == 0)
  #seuratObj2 <- subset(seuratObj, subset = e)
})

