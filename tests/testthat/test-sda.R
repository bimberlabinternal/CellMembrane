context("scRNAseq")

test_that("SDA works as expected", {
    seuratObj <- readRDS('../testdata/seuratOutput.rds')

    outputFolder <- tempdir()

    results <- RunSDA(seuratObj, outputFolder = outputFolder, numComps = 5)
    print(results)
    
    unlink(outputFolder)
})

