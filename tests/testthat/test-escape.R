context("scRNAseq")

test_that("escape works as expected", {
    seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
    
    #test minimal functionality
    testthat::expect_no_error(RunEscape(seuratObj, msigdbGeneSets = NULL, doPca = FALSE, customGeneSets = list("CD3E" = c("CD3E"))))
    #if no genes match any feature, escape will error on its own, and also the genes missing will warn the user. 
    testthat::expect_warning(testthat::expect_error(RunEscape(seuratObj, msigdbGeneSets = NULL, doPca = FALSE, customGeneSets = list("FakeGene" = c("FakeGene")))))
    
    #test full functionality
    seuratObj <- RunEscape(seuratObj, msigdbGeneSets = "H", customGeneSets = list("CD3E" = c("CD3E")))
    expect_equal(length(rownames(seuratObj@assays$escape.ssGSEA)), 51) #50 hallmark gene sets, 1 custom gene set
    expect_equal(max(seuratObj@assays$escape.ssGSEA$data[1]), 626, tolerance = 0.5)
})
