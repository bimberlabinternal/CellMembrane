context("scRNAseq")

test_that("escape works as expected", {
    seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
    seuratObj <- RunEscape(seuratObj)

    expect_equal(length(rownames(seuratObj@assays$escape.ssGSEA)), 50)
    expect_equal(max(seuratObj@assays$escape.ssGSEA$data[1]), 626)
})
