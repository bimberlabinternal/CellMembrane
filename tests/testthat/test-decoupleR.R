context("scRNAseq")

test_that("RunDecoupleR works as expected", {
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
  seuratObj <- RunDecoupleR(seuratObj)

  expect_true('tfsulm' %in% names(seuratObj@assays))
  expect_equal(nrow(seuratObj@assays$tfsulm), 726)
  expect_equal(ncol(seuratObj@assays$tfsulm), ncol(seuratObj))
})
