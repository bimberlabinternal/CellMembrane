context("scRNAseq")

test_that("hdWGCNA works as expected", {
    seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
    seuratObj <- RunHdWGCNA(seuratObj, groupName = 'S', groupBy = 'Phase', sampleGroupingVariables = c('Phase', 'ClusterNames_0.2'))

})
