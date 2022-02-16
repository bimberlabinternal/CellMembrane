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

  seuratList <- SplitSeurat(seuratObj, splitField = 'ClusterNames_0.2', minCellsToKeep = 200, appendLowFreqToOther = FALSE)
  expect_equal(3, length(seuratList))

  df <- AvgExpression(seuratObj, groupField = 'ClusterNames_0.4')
  expect_equal(ncol(df), 8)
  expect_equal(nrow(df), nrow(seuratObj) + 2)
  expect_equal(max(df[df$feature == 'LibrarySize', colnames(df) %in% 0:5]), 939113)

  seuratList <- SplitSeurat(seuratObj, splitField = 'ClusterNames_0.2', minCellsToKeep = 200)
  expect_equal(4, length(seuratList))

  seuratObj$Phase2 <- seuratObj$Phase
  seuratObj$Phase2[seuratObj$Phase == 'G1'] <- NA
  seuratList <- SplitSeurat(seuratObj, splitField = 'Phase2', minCellsToKeep = 200)
  expect_equal(3, length(seuratList))
  expect_equal(629, ncol(seuratList[['Other']]))
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
	
	seuratObj <- Seurat::CreateSeuratObject(dat)
	seuratObj <- AppendPerCellSaturation(seuratObj, molInfoFile)
	expect_equal(max(seuratObj$Saturation.RNA), 0.9765625)
	expect_equal(length(unique((seuratObj$Saturation.RNA))), 7126)
	expect_equal(sum(is.na(seuratObj$Saturation.RNA)), 0)
	expect_equal(sum(is.na(seuratObj$nReads_RNA)), 0)

	dat1 <- dat
	dat2 <- dat
	colnames(dat1) <- paste0('1234_', colnames(dat1))
  seuratObj <- Seurat::CreateSeuratObject(dat1)
  seuratObj$DatasetId <- 1234
  seuratObj$BarcodePrefix <- 1234
  
  colnames(dat2) <- paste0('12345_', colnames(dat2))
  seuratObj2 <- Seurat::CreateSeuratObject(dat2)
  seuratObj2$DatasetId <- 12345
  seuratObj2$BarcodePrefix <- 12345
  
  seuratObj3 <- CellMembrane::MergeSeuratObjs(seuratObjs = list('1234' = seuratObj, '12345' = seuratObj2), projectName = 'test')
  
  molInfoFileList <- list()
  molInfoFileList[['1234-RNA']] <- molInfoFile
  molInfoFileList[['12345-RNA']] <- molInfoFile

  seuratObj3 <- AppendPerCellSaturationInBulk(seuratObj3, molInfoFileList)
  expect_equal(max(seuratObj3$Saturation.RNA), 0.9765625)
  expect_equal(length(unique((seuratObj3$Saturation.RNA))), 7126)
  expect_equal(sum(is.na(seuratObj3$Saturation.RNA)), 0)
  expect_equal(sum(is.na(seuratObj3$nReads_RNA)), 0)
})
