context("scRNAseq")

test_that("escape works as expected", {
    seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))

    # if no genes match any feature, escape will error on its own, and also the genes missing will warn the user.
    testthat::expect_warning(testthat::expect_error(RunEscape(seuratObj, msigdbGeneSets = NULL, customGeneSets = list("FakeGene" = c("FakeGene")))))
    
    #test full functionality
    seuratObj <- RunEscape(seuratObj, msigdbGeneSets = "H", customGeneSets = list("CD3E" = c("CD3E")), performDimRedux = TRUE)
    expect_equal(length(rownames(seuratObj@assays$escape.ssGSEA)), 51) #50 hallmark gene sets, 1 custom gene set
    expect_equal(max(seuratObj@assays$escape.ssGSEA$counts[1]), 263, tolerance = 0.5)
    expect_equal(max(seuratObj@assays$escape.ssGSEA$data[1]), 0.112, tolerance = 0.01)
    expect_equal(max(seuratObj@assays$escape.ssGSEA$scale.data[1]), -0.624, tolerance = 0.01)
    
    expect_true('pca.escape.ssGSEA' %in% names(seuratObj@reductions))
    print(names(seuratObj@reductions))
    expect_true('escape.ssGSEA.umap' %in% names(seuratObj@reductions))
    expect_true('escape.ssGSEA.nn' %in% names(seuratObj@graphs))
    print(names(seuratObj@graphs))
})
