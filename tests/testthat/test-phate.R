context("scRNAseq")

test_that("PHATE works as expected", {
    if (!reticulate::py_available(initialize = TRUE)) {
        stop(paste0('Python/reticulate not configured. Run "reticulate::py_config()" to initialize python'))
    }

    if (!reticulate::py_module_available('phate')) {
        print('Phate module not found, debugging:')
        print(reticulate::py_list_packages())
        if ('phate' %in% reticulate::py_list_packages()$package) {
            tryCatch({
                reticulate::import('phate')
            }, error = function(e){
                print("Error with reticulate::import('phate')")
                print(conditionMessage(e))
                traceback()
            })
        }

        warning('The python phate module has not been installed!')
        expect_true(reticulate::py_module_available('phate'))
    }

    set.seed(CellMembrane::GetSeed())
    seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))

    seuratObj <- RunPHATE(seuratObj)
    testthat::expect_equal(ncol(seuratObj), 1557)
    
    Seurat::DimPlot(seuratObj, reduction = 'phate')
})

