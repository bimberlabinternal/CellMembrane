context("scRNAseq")

test_that("ClrNormalizeByGroup works as expected", {
    set.seed(CellMembrane::GetSeed())

    seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))

    # This is primarily to ensure it runs w/o error:
    seuratObj <- ClrNormalizeByGroup(seuratObj, groupingVar = 'ClusterNames_0.2', assayName = 'RNA', targetAssayName = 'ADT2')
    expect_equal(7.258214, max(Seurat::GetAssayData(seuratObj, assay = 'ADT2', slot = 'data')), tolerance = 0.001)

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

