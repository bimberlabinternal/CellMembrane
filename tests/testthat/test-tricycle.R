context("scRNAseq")

test_that("tricycle works as expected", {
    seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))

    seuratObj <- RunTricycle(seuratObj = seuratObj)

    expect_equal(6.28, max(seuratObj$tricyclePosition), tolerance = 0.05)
    expect_equal(169, sum(seuratObj$CCStage == 'S', na.rm = TRUE))
    expect_equal(133, sum(seuratObj$CCStage == 'G1.S', na.rm = TRUE))
})

