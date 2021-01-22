context("scRNAseq")

test_that("SingleR works as expected", {
    seuratObj <- readRDS('../testdata/seuratOutput.rds')

    results <- 'singleR.txt'
    singleRPrefix <- 'singleR.results'

    nGene <- nrow(seuratObj)
    nCell <- ncol(seuratObj)
    datasets <- c('hpca', 'blueprint', 'dice', 'monaco')
    seuratObj <- RunSingleR(seuratObj = seuratObj, resultTableFile = results, singlerSavePrefix = singleRPrefix, datasets = datasets)
    nGene2 <- nrow(seuratObj)
    nCell2 <- ncol(seuratObj)

    expect_equal(nGene, nGene2)
    expect_equal(nCell, nCell2)

    print(table(seuratObj$hpca.label))
    print(table(seuratObj$hpca.label.fine))

    for (dataset in datasets) {
		sr1 <- paste0(singleRPrefix, '.',dataset,'.singleR.rds')
		expect_true(file.exists(sr1))
		unlink(sr1)

		sr2 <- paste0(singleRPrefix, '.',dataset,'.singleR.fine.rds')
		expect_true(file.exists(sr2))
		unlink(sr2)
    }

    expect_equal(100, sum(seuratObj$hpca.label == 'NK_cell'))
    expect_equal(1435, sum(seuratObj$hpca.label == 'T_cells'))
    expect_equal(0, sum(seuratObj$hpca.label == 'B_cell'))
    expect_equal(0, sum(seuratObj$hpca.label == 'Neutrophils'))
    expect_equal(22, sum(seuratObj$hpca.label == 'Unknown'))

    expect_equal(153, sum(seuratObj$hpca.label.fine == 'T_cell:CD8+_Central_memory'))
    expect_equal(60, sum(seuratObj$hpca.label.fine == 'T_cell:CD8+_naive'))

    expect_equal(ncol(seuratObj), nrow(read.table(results, sep = '\t', header = T)))

    unlink(results)

    expect_equal(1252, sum(seuratObj$dice.label == 'NK cells'))
    expect_equal(228, sum(seuratObj$dice.label == 'T cells, CD4+'))
    expect_equal(3, sum(seuratObj$dice.label == 'Unknown'))

    for (dataset in datasets) {
        expect_true(!is.null(seuratObj@meta.data[[paste0(dataset,'.label')]]))
        expect_true(!is.null(seuratObj@meta.data[[paste0(dataset,'.label.fine')]]))
    }

    DimPlot_SingleR(seuratObj, datasets = datasets)

    Tabulate_SingleR(seuratObj, datasets = datasets)

    #saveRDS(seuratObj, file = '../testdata/singleR.rds')
})

