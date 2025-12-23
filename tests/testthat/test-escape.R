context("scRNAseq")

ensureEscapeCacheDir <- function() {
    # Debug saveRDS() issue:
    x <- tools::R_user_dir("escape", "cache")
    print(paste0('escape cache: ', x))
    if (! dir.exists(x)) {
        print('creating folder')
        dir.create(x, recursive = TRUE)
    }
}

test_that("escape works as expected", {
    ensureEscapeCacheDir()

    seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))

    # if no genes match any feature, escape will error on its own, and also the genes missing will warn the user.
    testthat::expect_warning(testthat::expect_error(RunEscape(seuratObj, msigdbGeneSets = NULL, customGeneSets = list(
      "FakeGene" = c("FakeGene"),
      "FakeGene2" = c("FakeGene2")
    ))))
    
    testthat::expect_warning(RunEscape(seuratObj, msigdbGeneSets = NULL, customGeneSets = list(
      "CD3G" = c("CD3G", "CD3E"),
      "CD4" = c("CD4"),
      "FakeGene2" = c("FakeGene2")
    )))
    
    #test full functionality
    seuratObj <- RunEscape(seuratObj, msigdbGeneSets = "H", customGeneSets = list("CD3" = c("CD3E", "CD3G"), "CD4" = c("CD4", "CD3E")), performDimRedux = FALSE, heatmapGroupingVars = 'ClusterNames_0.2')
    expect_equal(length(rownames(seuratObj@assays$escape.ssGSEA.H)), 50)
    expect_true('escape.ssGSEA.CustomGeneSet' %in% names(seuratObj@assays))

    print(seuratObj@assays$escape.ssGSEA.CustomGeneSet)
    print(rownames(seuratObj@assays$escape.ssGSEA.CustomGeneSet))
    expect_equal(length(rownames(seuratObj@assays$escape.ssGSEA.CustomGeneSet)), 2)

    seuratObj <- RunEscape(seuratObj, msigdbGeneSets = "H", performDimRedux = TRUE)
    expect_equal(max(seuratObj@assays$escape.ssGSEA.H$counts[1]), 263, tolerance = 0.5)
    expect_equal(max(seuratObj@assays$escape.ssGSEA.H$data[1]), 0.112, tolerance = 0.01)
    expect_equal(max(seuratObj@assays$escape.ssGSEA.H$scale.data[1]), -0.624, tolerance = 0.01)
    
    expect_true('pca.escape.ssGSEA.H' %in% names(seuratObj@reductions))
    print(names(seuratObj@reductions))
    expect_true('escape.ssGSEA.H.umap' %in% names(seuratObj@reductions))
    expect_true('escape.ssGSEA.H.nn' %in% names(seuratObj@graphs))
    print(names(seuratObj@graphs))

    seuratObj <- RunEscape(seuratObj, msigdbGeneSets = c("H", "C5" = "GO:BP"), performDimRedux = FALSE, escapeMethod = 'GSVA')
    print(names(seuratObj@assays))
    expect_true('escape.GSVA.H' %in% names(seuratObj@assays))
    expect_true('escape.GSVA.C5.GO.BP' %in% names(seuratObj@assays))

    seuratObj <- RunEscape(seuratObj, performDimRedux = FALSE, heatmapGroupingVars = 'ClusterNames_0.2')
    print(names(seuratObj@assays))
    expect_true('escape.ssGSEA.H' %in% names(seuratObj@assays))
    expect_true('escape.ssGSEA.C2.CP.KEGG' %in% names(seuratObj@assays))
    expect_true('escape.ssGSEA.C7.IMMUNESIGDB' %in% names(seuratObj@assays))
})

test_that("escape works with batches", {
    ensureEscapeCacheDir()

    seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))

    seuratObjNoBatch <- RunEscape(seuratObj, msigdbGeneSets = "H", performDimRedux = TRUE)
    print(names(seuratObjNoBatch@assays))
    expect_equal(max(seuratObjNoBatch@assays$escape.ssGSEA.H$counts[1]), 263, tolerance = 0.5)
    expect_equal(max(seuratObjNoBatch@assays$escape.ssGSEA.H$data[1]), 0.112, tolerance = 0.01)
    expect_equal(max(seuratObjNoBatch@assays$escape.ssGSEA.H$scale.data[1]), -0.624, tolerance = 0.01)

    seuratObj <- RunEscape(seuratObj, msigdbGeneSets = "H", performDimRedux = TRUE, maxBatchSize = 500)
    expect_equal(rownames(seuratObjNoBatch@assays$escape.ssGSEA.H$counts), rownames(seuratObj@assays$escape.ssGSEA.H$counts))
    expect_equal(max(seuratObj@assays$escape.ssGSEA.H$counts[1]), 263, tolerance = 0.5)
    expect_equal(max(seuratObj@assays$escape.ssGSEA.H$data[1]), 0.0886, tolerance = 0.01)
    expect_equal(max(seuratObj@assays$escape.ssGSEA.H$scale.data[1]), -0.613, tolerance = 0.01)
})