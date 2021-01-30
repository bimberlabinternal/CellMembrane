context("scRNAseq")

test_that("Serat processing works as expected", {
  set.seed(1234)

  outDir <- './'
  outPrefix <- paste0(outDir, 'testData')
  resolutionToUse <- 0.6

  seuratObj <- ReadAndFilter10xData('../testdata/CellRanger2/raw_gene_bc_matrices/cellRanger-3204293', datasetId = 'Set1', datasetName = 'datasetName', emptyDropNIters=5000)
  expect_true('BarcodePrefix' %in% colnames(seuratObj@meta.data))
  expect_true('DatasetId' %in% colnames(seuratObj@meta.data))
  expect_true('DatasetName' %in% colnames(seuratObj@meta.data))

  expect_equal(ncol(seuratObj), 3353, tolerance = 5)

  expect_equal(nrow(seuratObj), length(seuratObj@assays$RNA@meta.features$GeneId))
  geneIds <- GetGeneIds(seuratObj, c('HES4', 'CALML6'))
  names(geneIds) <- NULL
  expect_equal(geneIds, c('ENSMMUG00000001817', 'ENSMMUG00000012392'))

  gn <- c('HES4', 'CALML6', 'FAKE')
  geneIds <- GetGeneIds(seuratObj, gn, throwIfGenesNotFound = FALSE)
  expect_equal(names(geneIds), gn)
  names(geneIds) <- NULL
  expect_equal(geneIds, c('ENSMMUG00000001817', 'ENSMMUG00000012392', NA))

  #for speed, subset:
  cellsToUse <- sort(colnames(seuratObj))[1:500]
  seuratObj <- seuratObj[,cellsToUse]

  seuratObj <- FilterRawCounts(seuratObj)
  expect_equal(ncol(seuratObj), 485)
  
  seuratObj <- NormalizeAndScale(seuratObj)
  
  seuratObj <- RemoveCellCycle(seuratObj)
  tbl <- table(seuratObj$Phase)
  expect_equal(tbl[['G1']], 250)
  expect_equal(tbl[['G2M']], 92)
  expect_equal(tbl[['S']], 143)
  expect_equal(ncol(seuratObj), 485)
  
  vgFile <- 'variableGenes.txt'
  seuratObj <- RunPcaSteps(seuratObj, variableGeneTable = vgFile)
  expect_equal(ncol(seuratObj), 485)

  expect_equal(file.exists(vgFile), T)
  expect_equal(nrow(utils::read.table(vgFile, sep = '\t', header = F)), 2000)
  unlink(vgFile)

  seuratObj <- FindClustersAndDimRedux(seuratObj)
  expect_equal(ncol(seuratObj), 485)
  expect_equal(length(unique(seuratObj$ClusterNames_0.6)), 7)

  expect_equal(length(rownames(seuratObj@assays$RNA@scale.data)), length(rownames(seuratObj@assays$RNA@counts)))

  #Note: Seurat::PercentageFeatureSet returns 0-100.  our code is currently a fraction (0-1.0)
  expect_true(max(seuratObj$p.mito) < 1.0)
  expect_true(max(seuratObj$p.mito) > 0)

  seuratObj0 <- FindClustersAndDimRedux(seuratObj, minDimsToUse = 12)
  expect_equal(length(unique(seuratObj$ClusterNames_0.6)), 7)
  rm(seuratObj0)

  mf <- paste0(outPrefix, '.markers.txt')
  Find_Markers(seuratObj, identFields = c(resolutionToUse), outFile = mf, testsToUse = c('wilcox', 't'))

  expect_equal(nrow(utils::read.table(mf, sep = '\t', header = T)), 176, tolerance = 0)

  unlink(mf)

  sf <- paste0(outPrefix, '.summary.txt')
  WriteSummaryMetrics(seuratObj, file = sf)

  expect_equal(nrow(utils::read.table(sf, sep = '\t', header = T)), 2)

  unlink(sf)

  #At least execute this code, so over errors are caught
  PlotImmuneMarkers(seuratObj)

  #Note: if the expectations change, save this output as a reference:
  #seuratObjSS <- seuratObj[1:100]
  #saveRDS(seuratObjSS, file = '../testdata/seuratOutputSS.rds')
})