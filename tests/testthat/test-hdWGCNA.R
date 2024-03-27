context("scRNAseq")

test_that("hdWGCNA works as expected", {
    seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
    seuratObj <- RunHdWGCNA(seuratObj, wgcna_name = 'Test', groupName = 'S', groupBy = 'Phase', moduleConnectivityGroupBy = 'Phase', sampleGroupingVariables = c('Phase', 'ClusterNames_0.2'), reductionName = 'umap')
    
    expect_equal(max(seuratObj@misc$Test$module_scores), 0.9958171, tolerance = 0.005)

})
