context("scRNAseq")

test_that("ClrNormalizeByGroup works as expected", {
  set.seed(CellMembrane::GetSeed())

  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))

  # This is primarily to ensure it runs w/o error:
  seuratObj <- ClrNormalizeByGroup(seuratObj, groupingVar = 'ClusterNames_0.2', assayName = 'RNA', targetAssayName = 'ADT2')
  expect_equal(7.258214, max(Seurat::GetAssayData(seuratObj, assay = 'ADT2', slot = 'data')), tolerance = 0.001)

  seuratObj <- ClrNormalizeByGroup(seuratObj, groupingVar = 'ClusterNames_0.2', assayName = 'RNA', targetAssayName = 'ADT2', featureInclusionList = c(rownames(seuratObj@assays$RNA)[1:20]))
  expect_equal(3.135941, max(Seurat::GetAssayData(seuratObj, assay = 'ADT2', slot = 'data')), tolerance = 0.001)
})

test_that("AddNewMetaColumn works as expected", {
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
  
  seuratObj <- AddNewMetaColumn(seuratObj, varname = "binnedCounts", 
                                formulavector = c(nCount_RNA > 2000 ~ "High", 
                                                  nCount_RNA < 1000 ~ "Low"), 
                                defaultname = "Mid")
  expect_equal(1190, table(seuratObj$binnedCounts)[["High"]])
  expect_equal(97, table(seuratObj$binnedCounts)[["Low"]])
  expect_equal(270, table(seuratObj$binnedCounts)[["Mid"]])
  
  
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutputWithTCR.rds')))
  seuratObj <- AddNewMetaColumn(seuratObj, varname = "TRAVs", 
                                formulavector = c(TRA_V %in% c("TRAV1-2", "TRAV13-1") ~ "OfInterest"),
                                defaultname = "Meh")
  expect_equal(26, table(seuratObj$TRAVs)[["Meh"]])
  expect_equal(2, table(seuratObj$TRAVs)[["OfInterest"]])
  
  expect_error(AddNewMetaColumn(seuratObj, varname = "binnedCounts", 
                                 formulavector = c(nCount_RNA > 1e10 ~ "High", 
                                                   nCount_RNA < 1000 ~ "Low"), 
                                 defaultname = "Mid"))
})

test_that("GetMsigdbGeneSet works as expected", {
  #I think we might add GO:MF in the future, so this serves as a gotcha to check the codebase more fully to ensure compatibility
  testthat::expect_error(GetMsigdbGeneSet(msigdbGeneSets = "GO:MF"))

  #if this fails, then MsigDB added a "C9" category, and the Utils function GetMsigdbGeneSet needs to be updated to include C9
  testthat::expect_error(GetMsigdbGeneSet(msigdbGeneSets = "C9"))
  
  x <- GetMsigdbGeneSet(msigdbGeneSets = "GO:BP")
  testthat::expect_equal(length(x), 7658)
  
  y <- GetMsigdbGeneSet(msigdbGeneSets = "H")
  testthat::expect_equal(length(y), 50)
})

