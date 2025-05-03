context("scRNAseq")

test_that("escape works as expected", {
    seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))

    # if no genes match any feature, escape will error on its own, and also the genes missing will warn the user.
    testthat::expect_warning(testthat::expect_error(RunEscape(seuratObj, msigdbGeneSets = NULL, customGeneSets = list(
      "FakeGene" = c("FakeGene"),
      "FakeGene2" = c("FakeGene2")
    ))))
    
    testthat::expect_warning(RunEscape(seuratObj, msigdbGeneSets = NULL, customGeneSets = list(
      "CD3G" = c("CD3G", "CD3E"),
      "FakeGene2" = c("FakeGene2")
    )))
    
    #test full functionality
    seuratObj <- RunEscape(seuratObj, msigdbGeneSets = "H", customGeneSets = list("CD3" = c("CD3E", "CD3G"), "CD4" = c("CD4", "CD3E")), performDimRedux = FALSE)
    expect_equal(length(rownames(seuratObj@assays$escape.H)), 50)
    expect_true('escape.CustomGeneSet' %in% names(seuratObj@assays))

    print(seuratObj@assays$escape.CustomGeneSet)
    print(rownames(seuratObj@assays$escape.CustomGeneSet))
    expect_equal(length(rownames(seuratObj@assays$escape.CustomGeneSet)), 2)

    seuratObj <- RunEscape(seuratObj, msigdbGeneSets = "H", performDimRedux = TRUE)
    expect_equal(max(seuratObj@assays$escape.H$counts[1]), 263, tolerance = 0.5)
    expect_equal(max(seuratObj@assays$escape.H$data[1]), 0.112, tolerance = 0.01)
    expect_equal(max(seuratObj@assays$escape.H$scale.data[1]), -0.624, tolerance = 0.01)
    
    expect_true('pca.escape.H' %in% names(seuratObj@reductions))
    print(names(seuratObj@reductions))
    expect_true('escape.H.umap' %in% names(seuratObj@reductions))
    expect_true('escape.H.nn' %in% names(seuratObj@graphs))
    print(names(seuratObj@graphs))

    seuratObj <- RunEscape(seuratObj, msigdbGeneSets = c("H", "C5" = "GO:BP"), performDimRedux = TRUE)
    print(names(seuratObj@assays))
    expect_true('escape.H' %in% names(seuratObj@assays))
    expect_true('escape.C5.GO.BP' %in% names(seuratObj@assays))
})

test_that("escape works with batches", {
    seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))

    seuratObjNoBatch <- RunEscape(seuratObj, msigdbGeneSets = "H", performDimRedux = TRUE)
    expect_equal(max(seuratObjNoBatch@assays$escape.H$counts[1]), 263, tolerance = 0.5)
    expect_equal(max(seuratObjNoBatch@assays$escape.H$data[1]), 0.112, tolerance = 0.01)
    expect_equal(max(seuratObjNoBatch@assays$escape.H$scale.data[1]), -0.624, tolerance = 0.01)

    seuratObj <- RunEscape(seuratObj, msigdbGeneSets = "H", performDimRedux = TRUE, maxBatchSize = 500)
    expect_equal(rownames(seuratObjNoBatch@assays$escape.H$counts), rownames(seuratObj@assays$escape.H$counts))
    expect_equal(max(seuratObj@assays$escape.H$counts[1]), 263, tolerance = 0.5)
    expect_equal(max(seuratObj@assays$escape.H$data[1]), 0.0886, tolerance = 0.01)
    expect_equal(max(seuratObj@assays$escape.H$scale.data[1]), -0.613, tolerance = 0.01)
})