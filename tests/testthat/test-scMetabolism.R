context("scRNAseq")

test_that("scMetabolism works as expected", {
    seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
    seuratObj <- RunScMetabolism(seuratObj)
    
    expect_equal(length(rownames(seuratObj@assays$METABOLISM$score)), 85)
    expect_equal(max(seuratObj@assays$METABOLISM$score$AAACCTGAGCCAGGAT), 0.155963)

    expect_equal(length(rownames(seuratObj@misc$METABOLISM.KEGG)), 85)
    expect_equal(max(seuratObj@misc$METABOLISM.KEGG$AAACCTGAGCCAGGAT), 0.155963)
})
