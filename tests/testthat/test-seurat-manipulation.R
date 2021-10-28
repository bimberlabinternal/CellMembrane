context("scRNAseq")

test_that("Seurat-manipulation works as expected", {
  set.seed(CellMembrane::GetSeed())

  seuratObj <- readRDS('../testdata/seuratOutput.rds')

  seuratObj2 <- DownsampleSeurat(seuratObj, targetCells = 500)
  expect_equal(500, ncol(seuratObj2))
  
  seuratObj2 <- DownsampleSeurat(seuratObj, targetCells = 100, subsetFields = 'ClusterNames_0.2')
  expect_equal(493, ncol(seuratObj2))

  seuratObj2 <- DownsampleSeurat(seuratObj, targetCells = 100, subsetFields = c('ClusterNames_0.2'))
  expect_equal(493, ncol(seuratObj2))

  seuratObj2 <- DownsampleSeurat(seuratObj, targetCells = 100, subsetFields = c('ClusterNames_0.2', 'Phase'))
  expect_equal(1062, ncol(seuratObj2))

  seuratList <- SplitSeurat(seuratObj, splitField = 'ClusterNames_0.2')
  expect_equal(5, length(seuratList))

  seuratList <- SplitSeurat(seuratObj, splitField = 'ClusterNames_0.2', minCellsToKeep = 200)
  expect_equal(3, length(seuratList))

  df <- AvgExpression(seuratObj, groupField = 'ClusterNames_0.4')
  expect_equal(ncol(df), 8)
  expect_equal(nrow(df), nrow(seuratObj) + 2)
  expect_equal(max(df[df$feature == 'LibrarySize', colnames(df) %in% 0:5]), 939113)

  #seuratObj2 <- SubsetSeurat(seuratObj, expressionStrings = c('ClusterNames_0.2 == 0'))
  #e <- expression(ClusterNames_0.2 == 0)
  #seuratObj2 <- subset(seuratObj, subset = e)
})

context("scRNAseq")

test_that("Seurat-saturation works as expected", {
	molInfoFile <- '../testdata/512-5-molecule_info.h5'
	df <- DropletUtils::get10xMolInfoStats(molInfoFile)
	df <- data.frame(cellbarcode = df$cell, num.umis = df$num.umis, CountsPerCell = df$num.reads)
	df <- df[df$CountsPerCell > 100,]
	dat <- matrix(df$CountsPerCell, nrow = 1)
	colnames(dat) <- df$cellbarcode
	rownames(dat) <- c('Feat')
	
	seuratObj <- suppressWarnings(Seurat::CreateSeuratObject(dat))
	seuratObj <- AppendPerCellSaturation(seuratObj, molInfoFile)
	expect_equal(max(seuratObj$Saturation.RNA), 0.9765625)
	expect_equal(length(unique((seuratObj$Saturation.RNA))), 7126)

	colnames(dat) <- paste0('1234_', colnames(dat))
  seuratObj <- suppressWarnings(Seurat::CreateSeuratObject(dat))
  seuratObj$DatasetId <- 1234
  
  molInfoFileList <- list()
  molInfoFileList[[paste0(unique(seuratObj$DatasetId), '-RNA')]] <- molInfoFile

  seuratObj <- AppendPerCellSaturationInBulk(seuratObj, molInfoFileList)
  expect_equal(max(seuratObj$Saturation.RNA), 0.9765625)
  expect_equal(length(unique((seuratObj$Saturation.RNA))), 7126)
})