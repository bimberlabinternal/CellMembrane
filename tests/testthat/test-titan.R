context("scRNAseq")
library(Seurat)

test_that("SDA works as expected", {
    seuratObj <- readRDS('../testdata/seuratOutput.rds')

    outputFolder <- paste0(tempdir(), '/lda')
    if (dir.exists(outputFolder)) {
    	unlink(outputFolder, recursive = TRUE)
    }
    
    CellMembrane::DoLdaParameterScan(seuratObj, outputFolder = outputFolder)
    
})
