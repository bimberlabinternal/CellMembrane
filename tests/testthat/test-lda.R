context("scRNAseq")
library(Seurat)

test_that("LDA works as expected", {
    seuratObj <- readRDS('../testdata/seuratOutput.rds')

    results <- CellMembrane::DoLdaParameterScan(seuratObj, ntopics = c(5, 10, 15))
    expect_equal(length(results), 3)
    expect_equal(length(results[[1]]$topic_sums), 5)
    expect_equal(length(results[[2]]$topic_sums), 10)
    expect_equal(length(results[[3]]$topic_sums), 15)
})
