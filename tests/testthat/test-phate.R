context("scRNAseq")

test_that("PHATE works as expected", {
    if (!reticulate::py_available(initialize = TRUE)) {
        stop(paste0('Python/reticulate not configured. Run "reticulate::py_config()" to initialize python'))
    }

    if (!reticulate::py_module_available('phate')) {
        print('Phate module not found, debugging:')
        print(reticulate::py_list_packages())
        if ('phate' %in% reticulate::py_list_packages()$package) {
            reticulate::import('phate')
        }

        warning('The python phate module has not been installed!')
        return()
    }

    set.seed(CellMembrane::GetSeed())
    seuratObj <- readRDS('../testdata/seuratOutput.rds')

    seuratObj <- RunPHATE(seuratObj)
    
    Seurat::DimPlot(seuratObj, reduction = 'phate')
})

