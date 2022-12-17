context("scRNAseq")

test_that("SDA works as expected", {
    seuratObj <- readRDS('../testdata/seuratOutput.rds')

    outputFolder <- paste0(tempdir(), '/sda')
    if (dir.exists(outputFolder)) {
        unlink(outputFolder, recursive = TRUE)
    }

    results <- RunSDA(seuratObj, outputFolder = outputFolder, numComps = 5)
    print(results)
    
    unlink(outputFolder)
})

