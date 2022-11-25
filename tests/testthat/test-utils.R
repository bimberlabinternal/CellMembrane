context("scRNAseq")

test_that("ClrNormalizeByGroup works as expected", {
    set.seed(CellMembrane::GetSeed())

    seuratObj <- readRDS('../testdata/seuratOutput.rds')

    # This is primarily to ensure it runs w/o error:
    seuratObj <- ClrNormalizeByGroup(seuratObj, groupingVar = 'ClusterNames_0.2', assayName = 'RNA', targetAssayName = 'ADT2')
    expect_equal(7.258214, max(seuratObj@assays$ADT2@data))

})

