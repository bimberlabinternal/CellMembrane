library(DropletUtils)

context("scRNAseq")

test_that("Seurat-merge using emptyDropsCellRanger works", {
  if (!'emptyDropsCellRanger' %in% ls("package:DropletUtils")) {
    print('The installed DropletUtils lacks emptyDropsCellRanger, skipping test')
    return()
  }
	
  set.seed(CellMembrane::GetSeed())

  seuratObj <- ReadAndFilter10xData('../testdata/CellRanger2/raw_gene_bc_matrices/cellRanger-3204293', datasetId = 'Set1', datasetName = 'datasetName', emptyDropNIters=5000, useEmptyDropsCellRanger = T)

  expect_true('BarcodePrefix' %in% colnames(seuratObj@meta.data))
  expect_true('DatasetId' %in% colnames(seuratObj@meta.data))
  expect_true('DatasetName' %in% colnames(seuratObj@meta.data))

  expect_equal(ncol(seuratObj), 3244, tolerance = 5)
})

test_that("Seurat processing works as expected", {
  set.seed(CellMembrane::GetSeed())

  outDir <- './'
  outPrefix <- paste0(outDir, 'testData')
  resolutionToUse <- 0.6

  seuratObj <- ReadAndFilter10xData('../testdata/CellRanger2/raw_gene_bc_matrices/cellRanger-3204293', datasetId = 'Set1', datasetName = 'datasetName', emptyDropNIters=5000)
  expect_true('BarcodePrefix' %in% colnames(seuratObj@meta.data))
  expect_true('DatasetId' %in% colnames(seuratObj@meta.data))
  expect_true('DatasetName' %in% colnames(seuratObj@meta.data))

  expect_equal(ncol(seuratObj), 3353, tolerance = 5)

  expect_equal(nrow(seuratObj), length(slot(seuratObj@assays$RNA, GetAssayMetadataSlotName(seuratObj@assays$RNA))$GeneId))
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
  expect_equal(ncol(seuratObj), 487)
  
  seuratObj <- NormalizeAndScale(seuratObj)
  tbl <- table(seuratObj$Phase)
  expect_equal(tbl[['G1']], 190)
  expect_equal(tbl[['G2M']], 92)
  expect_equal(tbl[['S']], 205)
  expect_equal(ncol(seuratObj), 487)

  seuratObj <- RegressCellCycle(seuratObj)

  vgFile <- 'variableGenes.txt'
  seuratObj <- RunPcaSteps(seuratObj, variableGeneTable = vgFile)
  expect_equal(ncol(seuratObj), 487)

  expect_equal(file.exists(vgFile), T)
  expect_equal(nrow(utils::read.table(vgFile, sep = '\t', header = F)), 2000)
  unlink(vgFile)

  seuratObj <- FindClustersAndDimRedux(seuratObj)
  expect_equal(ncol(seuratObj), 487)
  expect_equal(length(unique(seuratObj$ClusterNames_0.6)), 6)

  #Note: Seurat::PercentageFeatureSet returns 0-100.  our code is currently a fraction (0-1.0)
  expect_true(max(seuratObj$p.mito) < 1.0)
  expect_true(max(seuratObj$p.mito) > 0)

  seuratObj0 <- FindClustersAndDimRedux(seuratObj, minDimsToUse = 12)
  expect_equal(length(unique(seuratObj$ClusterNames_0.6)), 6)
  rm(seuratObj0)

  mf <- paste0(outPrefix, '.markers.txt')
  dt <- Find_Markers(seuratObj, identFields = c(resolutionToUse), outFile = mf, testsToUse = c('wilcox', 't'), datasetName = 'Label')
  expect_equal(dt$x$caption, '<caption>Top DE Genes: Label</caption>')
  dt

  df <- utils::read.table(mf, sep = '\t', header = T)
  expect_equal(nrow(df), 856)
  expect_equal(sum(df$avg_logFC > 0.5), nrow(df))

  unlink(mf)

  sf <- paste0(outPrefix, '.summary.txt')
  WriteSummaryMetrics(seuratObj, file = sf)

  expect_equal(nrow(utils::read.table(sf, sep = '\t', header = T)), 2)

  unlink(sf)

  P1 <- Seurat::DimPlot(seuratObj)
  seuratObj$CloneName <- NA

  seuratObj$CloneName[rep(c( rep(FALSE, 5), TRUE ), ncol(seuratObj))] <- 'Clone1'
  seuratObj$CloneName[rep(c( rep(FALSE, 6), TRUE ), ncol(seuratObj))] <- 'Clone2'

  HighlightCellsOnSeuratPlot(seuratPlot = P1, seuratObj = seuratObj, cellSelectField = 'CloneName')

  PlotSeuratVariables(seuratObj, xvar = 'ClusterNames_0.6', yvar = 'ClusterNames_0.2')

  #Note: if the expectations change, save this output as a reference:
  #seuratObjSS <- seuratObj[1:100]
  #saveRDS(seuratObjSS, file = '../testdata/seuratOutputSS.rds')
})

test_that("Seurat SCTransform works as expected", {
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
  seuratObjSCT <- CreateSeuratObj(seuratData = Seurat::GetAssayData(seuratObj, assay = 'RNA', slot = 'counts'), datasetId = '1234', datasetName = 'Set1')

  seuratObjSCT <- NormalizeAndScale(seuratObjSCT, useSCTransform = T)
  expect_equal(length(rownames(Seurat::GetAssayData(seuratObjSCT, assay = 'SCT', slot = 'scale.data'))), length(rownames(Seurat::GetAssayData(seuratObjSCT, assay = 'SCT', slot = 'counts'))))
  expect_equal(ncol(seuratObjSCT), ncol(seuratObj))
})


test_that("Seurat CellCycleScoring_UCell", {
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
  seuratObj <- CellCycleScoring_UCell(seuratObj, outputFieldName = 'PhaseUCell')
  expect_equal(sum(seuratObj$PhaseUCell == 'S'), 619)
  expect_equal(sum(seuratObj$PhaseUCell == 'G2M'), 660)
  expect_equal(sum(seuratObj$PhaseUCell == 'G1'), 278)
})


test_that("LogNormalizeUsingAlternateAssay works as expected", {
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))

  assayToAdd <- Seurat::GetAssayData(seuratObj, assay = 'RNA', layer = 'counts')
  assayToAdd <- floor(assayToAdd[1:10,] / 5)
  
  rownames(assayToAdd) <- paste0('Feature', LETTERS[1:10])

  seuratObj[['Norm']] <- Seurat::CreateAssayObject(assayToAdd)

  seuratObj <- LogNormalizeUsingAlternateAssay(seuratObj, assayToNormalize = 'Norm', assayForLibrarySize = 'RNA')

  nd <- Seurat::GetAssayData(seuratObj, assay = 'Norm', layer = 'data')
  expect_equal(max(nd[,4]), 3.442982, tolerance = 0.000001)
  expect_equal(max(nd[,101]), 2.823479, tolerance = 0.000001)
})