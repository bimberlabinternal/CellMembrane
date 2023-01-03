context("scRNAseq")
library(Seurat)

test_that("SDA works as expected", {
    seuratObj <- readRDS('../testdata/seuratOutput.rds')

    outputFolder <- paste0(tempdir(), '/sda')
    if (dir.exists(outputFolder)) {
        unlink(outputFolder, recursive = TRUE)
    }

    # NOTE: these values are designed to highly subset the matrix to make the input very small and fast to run. This is just designed to show code runs w/o errors.
    seuratObj <- subset(seuratObj, cells = colnames(seuratObj)[1:100])
    sdaResults <- RunSDA(seuratObj, outputFolder = outputFolder, numComps = 2, minFeatureCount = 1000, max_iter = 50)
    print(utils::str(sdaResults))

    expect_true('component_statistics' %in% names(sdaResults))
    expect_true('Features' %in% names(sdaResults))
    expect_true('CellBarcodes' %in% names(sdaResults))

    expect_equal(length(sdaResults$scores), 200)

    seuratObj <- SDAToSeuratMetadata(seuratObj, sdaResults)
    expect_true('SDA_1' %in% names(seuratObj@meta.data))
    expect_true('SDA_2' %in% names(seuratObj@meta.data))

    seuratObj <- SDAToSeuratReduction(seuratObj, sdaResults)
    expect_true('sda' %in% names(seuratObj@reductions))

    unlink(outputFolder)
})

