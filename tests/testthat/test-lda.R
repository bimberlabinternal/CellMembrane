context("scRNAseq")
library(Seurat)

test_that("LDA works as expected", {
    seuratObj <- readRDS('../testdata/seuratOutput.rds')

    results <- CellMembrane::DoLdaParameterScan(seuratObj, ntopics = c(5, 10))
})
