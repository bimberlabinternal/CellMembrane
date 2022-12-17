context("scRNAseq")
library(Seurat)

test_that("SDA works as expected", {
    seuratObj <- readRDS('../testdata/seuratOutput.rds')

    outputFolder <- paste0(tempdir(), '/sda')
    if (dir.exists(outputFolder)) {
        unlink(outputFolder, recursive = TRUE)
    }

    seuratObj <- subset(seuratObj, cells = colnames(seuratObj[1:100]))
    results <- RunSDA(seuratObj, outputFolder = outputFolder, numComps = 2, minAsinhThreshold = 8)
    print(results)
    
    unlink(outputFolder)
})

