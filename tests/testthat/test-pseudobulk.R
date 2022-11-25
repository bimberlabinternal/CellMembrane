context("scRNAseq")

test_that("Pseudobulk works", {
	seuratObj <- readRDS('../testdata/seuratOutput.rds')

	pseudo <- PseudobulkSeurat(seuratObj, groupFields = c('ClusterNames_0.2'))
	expect_equal(length(unique(seuratObj$ClusterNames_0.2)), nrow(pseudo@meta.data))
	expect_equal(nrow(seuratObj@assays$RNA), nrow(pseudo@assays$RNA))
	
	pseudo2 <- PseudobulkSeurat(seuratObj, groupFields = c('ClusterNames_0.4'), assays = c('RNA'))
	expect_equal(length(unique(seuratObj$ClusterNames_0.4)), nrow(pseudo2@meta.data))
	expect_equal(nrow(seuratObj@assays$RNA), nrow(pseudo2@assays$RNA))
})
