context("spatialRNAseq")

test_that("Q3 normalization works", {
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
  
  testthat::expect_no_error(seuratObj[["RNA"]])
  seuratObj<- RUVg_Housekeeping_Normalization(seuratObj)
  testthat::expect_no_error(seuratObj[["RUVg"]])
  expect_equal(max(seuratObj@assays$RUVg@counts), 2335)
  testthat::expect_error(object = NanoString_Housekeeping_Normalization(seuratObj),"Geometric mean of housekeeping counts for one or more samples is zero")
  testthat::expect_error(object = Q3_Normalization(seuratObj),"Error: 75th percentile is zero for some samples.")
  seuratObj@assays$RNA@counts <- SeuratObject::as.sparse(as.matrix( seuratObj@assays$RNA@counts) +1 )
  
  seuratObj<- Q3_Normalization(seuratObj)
  testthat::expect_no_error(seuratObj[["Q3"]])
  testthat::expect_equal(max(seuratObj@assays$Q3@counts), 1663)

  seuratObj<- NanoString_Housekeeping_Normalization(seuratObj)
  testthat::expect_no_error(seuratObj[["Housekeeping"]])
  testthat::expect_equal(round(max(seuratObj@assays$Housekeeping@counts)), 76)

})